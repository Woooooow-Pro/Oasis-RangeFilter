#pragma once

#include <cmath>
#include <cstdint>
#include <list>
#include <tuple>
#include <vector>

#include "learned_rf/bitset.h"
#include "learned_rf/learned_rf.h"
#include "proteus/proteus.h"

namespace oasis_plus {
class FilterBuilder {
 private:
  static const uint64_t kCost = 129;
  static const uint64_t kMeta_LRF = 64;
  static const uint64_t kMeta_Biset = 32;
  static constexpr long double LN2_2 = -M_LN2 * M_LN2;

 public:
  FilterBuilder(double bpk, uint32_t block_size, const size_t max_qlen = 10);
  ~FilterBuilder() = default;

  auto get_filter_types() -> std::vector<uint32_t>&& {
    return std::move(filter_types_);
  }
  auto get_begins() -> std::vector<uint64_t>&& { return std::move(begins_); }
  auto get_ends() -> std::vector<uint64_t>&& { return std::move(ends_); }

  auto get_learned_rf() -> LearnedRF* { return learned_rf_; }
  auto get_proteus() -> Proteus* { return proteus_; }

 private:
  void build_block_lists(const std::vector<uint64_t>& locations);

  auto cal_proteus_fpr(double mem_budget, size_t trie_len, size_t bf_len)
      -> long double;
  auto cal_used_bytes() const -> size_t;

  inline auto cal_mp_sz(double mem_budget, size_t key_sz, size_t m) -> uint64_t;
  inline void shink_index(const uint64_t threshold,
                          std::vector<uint64_t>& begins,
                          std::vector<uint64_t>& ends,
                          std::vector<size_t>& interval_sz);

  // experiment function
 public:
  void build(const std::vector<uint64_t>& keys);
  void segmentation(const std::vector<uint64_t>& keys);
  void build_index(const uint64_t threshold, const std::vector<uint64_t>& keys,
                   std::vector<uint64_t>& begins, std::vector<uint64_t>& ends,
                   std::vector<size_t>& interval_sz);

  auto build_seq(const std::vector<uint64_t>& begins,
                 const std::vector<uint64_t>& ends,
                 const std::vector<size_t>& interval_sz) -> std::vector<size_t>;
  inline auto build_seq(const std::vector<uint64_t>& density)
      -> std::vector<size_t>;

  void build_proteus(const std::vector<uint64_t>& keys,
                     const std::tuple<size_t, size_t, size_t>& conf,
                     double bpk);
  void get_positions(const std::vector<uint64_t>& keys, double bpk,
                     uint64_t delta_sum, size_t lrf_sz, size_t lrf_interval_num,
                     std::vector<uint64_t>& learned_pos,
                     std::vector<uint64_t>& proteus_keys);

 private:
  // metadata of the building step
  double bpk_;
  size_t max_klen_ = 64;
  size_t max_qlen_ = 10;
  size_t M_;
  std::vector<size_t> qk_dists_;
  std::vector<size_t> key_prefixes_;
  std::vector<uint64_t> threshold_set_;

  // OasisPlus metadata
  std::vector<uint32_t> filter_types_;
  std::vector<uint64_t> begins_;
  std::vector<uint64_t> ends_;

  // learned filter metadata
  uint16_t block_sz_;
  uint16_t last_block_sz_;
  size_t bitmap_sz_ = 0;
  uint8_t* bitmap_ptr_ = nullptr;
  std::vector<uint64_t> accumulate_interval_sz_;
  std::vector<BitSet> block_lists_;
  std::vector<uint64_t> block_bias_;

  // proteus
  Proteus* proteus_;
  // learned filter
  LearnedRF* learned_rf_;
};
}  // namespace oasis_plus
