#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <queue>
#include <vector>

#include "bitset.h"
#include "common/type.h"

namespace oasis {

enum class QueryPosStatus { OUT_OF_SCOPE, EXIST, NO_IDEA };

template <typename Key>
class Models {
  using Keys = std::vector<Key>;

 public:
  Models() = default;
  virtual ~Models() = default;

  /**
   * For given query (point/range), return the (slope, bias) of the model
   * return (-1, 0): query crosses the models
   * return (0, 2): query out of scop;
   */
  virtual auto query(const Key &key, BitPos &result) -> QueryPosStatus = 0;
  /* [l_key, r_key], and return [l_pos, r_pos] */
  virtual auto query(const Key &left_key, const Key &right_key,
                     std::pair<BitPos, BitPos> &result) -> QueryPosStatus = 0;

  /* return the estimate distribution of all the key */
  virtual auto get_locations(const Keys &keys) -> BitsPos = 0;

  virtual auto size() const -> size_t = 0;

  virtual auto serialize() const -> std::pair<uint8_t *, size_t> = 0;

  /*
   * For any type of model you must implement a deserializer!
   * static auto deserialize(uint8_t *ser) -> std::unique_ptr<Models>;
   */
};

template <typename Key>
class CDFModels : public Models<Key> {
  using Keys = std::vector<Key>;
  using MinHeap = std::priority_queue<Key, Keys, std::greater<Key>>;

  static const uint64_t kCost =
      static_cast<uint64_t>(sizeof(Key) * 2 + sizeof(BitPos)) << 3;

 public:
  CDFModels(Keys &begins, Keys &ends, BitsPos &accumulate_alphas)
      : begins_(std::move(begins)),
        ends_(std::move(ends)),
        accumulate_alphas_(std::move(accumulate_alphas)) {
    assert(begins_.size() == ends_.size());
    assert(accumulate_alphas_.size() == begins_.size() + 1);
  }

  CDFModels(double bpk, size_t elem_per_block, const Keys &keys) {
    size_t nkeys = keys.size();
    double kMemBudget = bpk * nkeys;
    size_t M = static_cast<size_t>(kMemBudget / kCost);

    MinHeap min_heap;

    for (size_t i = 0; i < nkeys - 1; ++i) {
      Key diff = keys[i + 1] - keys[i];
      if (min_heap.size() >= M) {
        if (min_heap.top() > diff) {
          continue;
        }
        min_heap.pop();
      }
      min_heap.push(diff);
    }

    std::queue<Key> threshold_set;
    Key delta_sum = Key();
    while (!min_heap.empty()) {
      Key tmp = min_heap.top();
      min_heap.pop();

      threshold_set.push(tmp);
      delta_sum += tmp;
    }

    delta_sum = keys.back() - keys[0] - delta_sum;
    double remain_bpk = bpk - 2 - 64.0L / elem_per_block;

    Key threshold = get_threshold(remain_bpk, delta_sum, nkeys, threshold_set);
    build_indices(threshold, remain_bpk, keys);
  }

  ~CDFModels() = default;

  auto get_locations(const Keys &keys) -> BitsPos override {
    size_t nkeys = keys.size();

    BitsPos locations;

    size_t idx_iter = 0;
    auto params = get_params(idx_iter);
    for (size_t i = 1; i < nkeys - 1; ++i) {
      if (keys[i] >= ends_[idx_iter]) {
        params = get_params(++idx_iter);
      } else if (keys[i] > begins_[idx_iter]) {
        locations.emplace_back(get_location(keys[i], params));
      }
    }
    return locations;
  }

  auto query(const Key &key, BitPos &result) -> QueryPosStatus override {
    if (key < begins_[0] || key > ends_.back()) {
      return QueryPosStatus::OUT_OF_SCOPE;
    }

    size_t idx =
        std::distance(begins_.begin(),
                      std::upper_bound(begins_.begin(), begins_.end(), key)) -
        1;
    if (ends_[idx] < key) {
      return QueryPosStatus::OUT_OF_SCOPE;
    }
    if (begins_[idx] == key || ends_[idx] == key) {
      return QueryPosStatus::EXIST;
    }

    result = get_location(key, get_params(idx));
    return QueryPosStatus::NO_IDEA;
  }

  auto query(const Key &left_key, const Key &right_key,
             std::pair<BitPos, BitPos> &result) -> QueryPosStatus override {
    if (left_key > ends_.back() || right_key < begins_[0]) {
      return QueryPosStatus::OUT_OF_SCOPE;
    }

    size_t idx = std::distance(begins_.begin(),
                               std::upper_bound(begins_.begin(), begins_.end(),
                                                left_key)) -
                 1;

    if (right_key < begins_[idx + 1] && left_key > ends_[idx]) {
      return QueryPosStatus::OUT_OF_SCOPE;
    }

    if (!(left_key > begins_[idx] && right_key < ends_[idx])) {
      return QueryPosStatus::EXIST;
    }

    auto params = get_params(idx);

    result.first = get_location(left_key, params);
    result.second = get_location(right_key, params);
    return QueryPosStatus::NO_IDEA;
  }

  auto size() const -> size_t override {
    size_t idx_sz = begins_.size();
    return sizeof(size_t) +           /* # elements in index_ */
           sizeof(Key) * 2 * idx_sz + /* begins_ + ends_*/
           sizeof(BitPos) * idx_sz;   /* alphas */
  }

  auto serialize() const -> std::pair<uint8_t *, size_t> override {
    size_t idx_sz = begins_.size();

    size_t size = sizeof(size_t) +           /* # elements in index_ */
                  sizeof(Key) * 2 * idx_sz + /* begins_ + ends_*/
                  sizeof(BitPos) * idx_sz;   /* alphas */

    uint8_t *out = new uint8_t[size];
    uint8_t *pos = out;

    memcpy(pos, &idx_sz, sizeof(size_t));
    pos += sizeof(size_t);

    memcpy(pos, begins_.data(), sizeof(Key) * idx_sz);
    pos += sizeof(Key) * idx_sz;

    memcpy(pos, ends_.data(), sizeof(Key) * idx_sz);
    pos += sizeof(Key) * idx_sz;

    memcpy(pos, accumulate_alphas_.data() + 1, sizeof(BitPos) * idx_sz);

    return {out, size};
  }

  static auto deserialize(uint8_t *ser) -> std::unique_ptr<CDFModels> {
    size_t idx_len;
    memcpy(&idx_len, ser, sizeof(size_t));
    ser += sizeof(size_t);

    size_t index_sz = sizeof(uint64_t) * idx_len;
    std::vector<uint64_t> begins(idx_len);
    memcpy(begins.data(), ser, index_sz);
    ser += index_sz;

    std::vector<uint64_t> ends(idx_len);
    memcpy(ends.data(), ser, index_sz);
    ser += index_sz;

    std::vector<uint64_t> accumulate_nkeys(idx_len + 1, 0);
    memcpy(accumulate_nkeys.data() + 1, ser, index_sz);

    return std::make_unique<CDFModels>(begins, ends, accumulate_nkeys);
  }

 private:
  auto get_threshold(double bpk, Key delta_sum, size_t nkeys,
                     std::queue<Key> &threshold_set) -> Key {
    double param = kCost * 1.0L / nkeys;
    double min_rho = std::numeric_limits<double>::max();
    size_t m = threshold_set.size();
    Key best_threshold;

    auto get_best = [&](const Key &threshold) {
      uint64_t mp_sz = std::ceil(std::pow(2, bpk - param * (m + 1)) * nkeys);
      double rho = static_cast<double>(delta_sum) * delta_sum / mp_sz;
      if (rho <= min_rho) {
        min_rho = rho;
        best_threshold = threshold;
      }
    };

    while (!threshold_set.empty()) {
      Key threshold = threshold_set.front();

      get_best(threshold);
      while (!threshold_set.empty() && threshold_set.front() == threshold) {
        threshold_set.pop();
        delta_sum += threshold;
      }
      m = threshold_set.size();
    }

    get_best(delta_sum);

    return best_threshold;
  }

  auto get_params(size_t idx) -> std::vector<double> {
    auto begin = static_cast<double>(begins_[idx]);
    auto end = static_cast<double>(ends_[idx]);
    auto l_bound = static_cast<double>(accumulate_alphas_[idx]);
    auto r_bound = static_cast<double>(accumulate_alphas_[idx + 1]);
    return {end - begin, r_bound - l_bound, end * l_bound - begin * r_bound};
  }

  auto get_location(const double &key, const std::vector<double> &params)
      -> BitPos {
    return (params[1] * static_cast<double>(key) + params[2]) / params[0];
  }

  void build_indices(const Key threshold, const double bpk, const Keys &keys) {
    size_t nkeys = keys.size();
    // set avg_range to default value
    double avg_range = 0.0;

    begins_.clear();
    ends_.clear();

    std::vector<bool> interval_sz;
    uint32_t cnt = 0;
    /** build the indices */
    begins_.emplace_back(keys[0]);
    for (size_t i = 0; i < nkeys - 1; ++i) {
      if (keys[i + 1] - keys[i] >= threshold) {
        ends_.emplace_back(keys[i++]);

        interval_sz.emplace_back(cnt == 0);
        avg_range +=
            cnt == 0 ? 0.0 : static_cast<double>(ends_.back() - begins_.back());
        cnt = 0;
        begins_.emplace_back(keys[i]);
      } else {
        ++cnt;
      }
    }
    ends_.emplace_back(keys.back());
    interval_sz.emplace_back(cnt == 0);
    avg_range +=
        cnt == 0 ? 0.0 : static_cast<double>(ends_.back() - begins_.back());

    accumulate_alphas_.clear();
    accumulate_alphas_.emplace_back(0);

    /** Build the alpha array */
    BitPos bit_array_range =
        std::ceil(pow(2, bpk - kCost * 1.0L / nkeys * ends_.size()) * nkeys);

    for (size_t i = 0; i < begins_.size(); ++i) {
      if (interval_sz[i]) {
        accumulate_alphas_.emplace_back(accumulate_alphas_.back());
        continue;
      }
      BitPos alpha = std::ceil(static_cast<double>(ends_[i] - begins_[i]) /
                               avg_range * bit_array_range);
      if (alpha == 0) {
        alpha = 1;
      }
      accumulate_alphas_.emplace_back(accumulate_alphas_.back() + alpha);
    }
  }

 private:
  Keys begins_, ends_;
  BitsPos accumulate_alphas_;
};

}  // namespace oasis