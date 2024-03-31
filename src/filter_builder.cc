#include "filter_builder.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>

#include "proteus/config.h"
#include "proteus/modeling.h"
#include "proteus/prefixbf.h"
#include "util.h"

namespace oasis_plus {
FilterBuilder::FilterBuilder(double bpk, uint32_t block_size,
                             const size_t max_qlen)
    : bpk_(bpk),
      max_klen_(64),
      max_qlen_(max_qlen),
      key_prefixes_(max_klen_ + 1, 0),
      block_sz_(block_size),
      proteus_(nullptr),
      learned_rf_(nullptr) {}

void FilterBuilder::build_block_lists(const std::vector<uint64_t>& locations) {
  last_block_sz_ = locations.size() % block_sz_;

  std::vector<uint8_t> compressed_bitmap;

  block_bias_.clear();
  uint64_t low_bound = locations[0];
  block_bias_.emplace_back(low_bound);
  std::vector<uint64_t> cur_batch;
  for (size_t idx = 1; idx < locations.size();) {
    cur_batch.clear();
    cur_batch.emplace_back(0);
    for (size_t cnt = 1; cnt < block_sz_ && idx < locations.size(); ++cnt) {
      cur_batch.emplace_back(locations[idx++] - low_bound);
    }
    if (idx < locations.size()) {
      block_bias_.emplace_back(locations[idx++]);
    } else {
      block_bias_.emplace_back(locations.back());
    }
    std::vector<uint8_t> batch_block =
        BitSet::build(cur_batch, block_bias_.back() - low_bound);
    compressed_bitmap.insert(compressed_bitmap.end(), batch_block.begin(),
                             batch_block.end());
    low_bound = block_bias_.back();
  }

  bitmap_sz_ = compressed_bitmap.size();
  bitmap_ptr_ = new uint8_t[bitmap_sz_];
  memcpy(bitmap_ptr_, compressed_bitmap.data(), bitmap_sz_ * sizeof(uint8_t));

  size_t nbatches = block_bias_.size() - 1;
  uint8_t* pos = bitmap_ptr_;
  for (size_t i = 0; i < nbatches; ++i) {
    block_lists_.emplace_back(i == nbatches - 1 ? last_block_sz_ : block_sz_,
                              block_bias_[i + 1] - block_bias_[i], pos);
    pos += block_lists_.back().size();
  }
}

void FilterBuilder::shink_index(const uint64_t threshold,
                                std::vector<uint64_t>& begins,
                                std::vector<uint64_t>& ends,
                                std::vector<size_t>& interval_sz) {
  std::vector<uint64_t> tmp_begins, tmp_ends;
  std::vector<size_t> tmp_sz;
  tmp_begins.emplace_back(begins[0]);
  tmp_sz.emplace_back(interval_sz[0]);
  for (size_t i = 1; i < begins.size(); ++i) {
    if (begins[i] - ends[i - 1] < threshold) {
      tmp_sz.back() += interval_sz[i];
    } else {
      tmp_ends.emplace_back(ends[i - 1]);
      tmp_begins.emplace_back(begins[i]);
      tmp_sz.emplace_back(interval_sz[i]);
    }
  }
  tmp_ends.emplace_back(ends.back());
  begins = std::move(tmp_begins);
  ends = std::move(tmp_ends);
  interval_sz = std::move(tmp_sz);
}

auto FilterBuilder::cal_used_bytes() const -> size_t {
  // calculate metadata size
  size_t interval_num = begins_.size();
  size_t metadata_sz = sizeof(uint64_t) /* number of the interval */
                       + sizeof(uint64_t) * interval_num * 2 /* indices size*/
                       + BIT2BYTE(interval_num) /* filter types bitmap */;

  size_t learned_sz =
      sizeof(uint16_t)   /** block_sz_ */
      + sizeof(uint16_t) /** last_block_sz_ */
      + sizeof(uint32_t) /** # accumulate_interval_sz_ */
      + sizeof(uint32_t) /** # block_lists_ */
      + sizeof(size_t)   /** size of compressed bitmap */
      + sizeof(uint64_t) * (accumulate_interval_sz_.size() -
                            1) /** size of accumulate interval size array */
      + sizeof(uint64_t) * block_bias_.size() /** size of block bias array */
      + bitmap_sz_;                           /** size of all bitset */

  return learned_sz + metadata_sz;
}

auto FilterBuilder::cal_mp_sz(double bpk, size_t key_sz, size_t m) -> uint64_t {
  return pow(2, bpk - 1.0L * kMeta_LRF * m / key_sz) * key_sz;
}

auto FilterBuilder::cal_proteus_fpr(double mem_budget, size_t trie_len,
                                    size_t bf_len) -> long double {
  size_t empty_queries = key_prefixes_[max_klen_];
  size_t solved_in_trie = trie_len == 0 ? 0 : key_prefixes_[trie_len - 1];

  size_t n = key_prefixes_[bf_len - 1];
  size_t nhf = static_cast<size_t>(round(M_LN2 * mem_budget / n));
  nhf = (nhf == 0 ? 1 : nhf);
  nhf = std::min(static_cast<size_t>(MAX_PBF_HASH_FUNCS), nhf);

  long double prefix_query_fpr =
      pow((1.0L - exp(-((nhf * n * 1.0L) / mem_budget))), nhf);

  size_t resolved_in_bf = bf_len == 0 ? 0 : key_prefixes_[bf_len - 1];
  resolved_in_bf -= solved_in_trie;
  size_t base = resolved_in_bf / (max_qlen_ + 1);

  long double cumulative_fpp = empty_queries - solved_in_trie - resolved_in_bf;
  if (trie_len >= max_klen_ - max_qlen_) {
    size_t query_times = 1 << (bf_len - trie_len + 1);
    cumulative_fpp +=
        base * static_cast<long double>(trie_len - max_klen_ + max_qlen_ + 1) *
        (1.0L - pow(1.0L - prefix_query_fpr, query_times));
  }
  for (size_t q_len = 0; q_len <= max_qlen_ && max_klen_ - q_len > trie_len;
       ++q_len) {
    size_t query_times =
        max_klen_ - q_len >= bf_len ? 1 : 1 << (bf_len + q_len - max_klen_);
    cumulative_fpp += base * (1.0L - pow(1.0L - prefix_query_fpr, query_times));
  }

  return cumulative_fpp;
}

/**
 * Find best LRF then decide the filter for each interval.................
 */
void FilterBuilder::build(const std::vector<uint64_t>& keys) {
  segmentation(keys);

  uint64_t delta_sum =
      keys.back() - keys[0] -
      std::accumulate(threshold_set_.begin() + M_, threshold_set_.end(), 0ULL);

  size_t nkeys = keys.size();
  size_t m = nkeys - M_;
  double bpk = bpk_ - (2 + kMeta_Biset * 1.0L / block_sz_);

  double min_rho = std::numeric_limits<double>::max();
  size_t best_lrf_idx = M_;
  size_t best_delta_sum = delta_sum;
  for (size_t set_idx = M_; set_idx < threshold_set_.size();) {
    double bpk_slash = bpk - 1.0L * kCost * m / nkeys;
    uint64_t mp_sz = cal_mp_sz(bpk_slash, nkeys, m);
    double rho = static_cast<double>(delta_sum) * delta_sum / mp_sz;
    if (rho <= min_rho) {
      min_rho = rho;
      best_lrf_idx = set_idx;
      best_delta_sum = delta_sum;
    }

    // Update set idx
    size_t pre_idx = set_idx;
    uint64_t cur_threshold = threshold_set_[set_idx];
    while (set_idx < threshold_set_.size() &&
           threshold_set_[set_idx] == cur_threshold) {
      ++set_idx;
    }
    m = nkeys - set_idx;
    delta_sum += (set_idx - pre_idx) * cur_threshold;
  }

  std::vector<size_t> interval_sz;
  build_index(threshold_set_[best_lrf_idx], keys, begins_, ends_, interval_sz);

  std::vector<uint64_t> begins = begins_;
  std::vector<uint64_t> ends = ends_;

  /**
   * decide the intervals encoding method
   */
  ProteusModeling proteus_model(max_klen_, max_qlen_);
  proteus_model.set_key_prefixes(key_prefixes_);
  std::tuple<size_t, size_t, size_t> single_proteus_conf =
      proteus_model.modeling(keys, bpk_);

  std::tuple<size_t, size_t, size_t> best_pro_conf = single_proteus_conf;
  std::vector<bool> best_indicator(1, false);
  double min_fpr = proteus_model.get_min_fpp() * (keys.back() - keys[0]);
  size_t best_idx = best_lrf_idx;
  size_t best_lrf_sz = nkeys;

  m = ends.size();
  uint64_t batch_sum = best_delta_sum;

  double bpk_slash = bpk_ - 1.0L * kCost * m / nkeys;
  size_t pre_bpk_slash_1 = std::floor(bpk_slash);
  double pre_bpk_slash_2 = bpk_slash;
  std::tuple<size_t, size_t, size_t> cur_pro_conf =
      proteus_model.modeling(keys, bpk_slash);
  long double proteus_fpr = proteus_model.get_min_fpp();
  for (size_t set_idx = best_lrf_idx; set_idx < threshold_set_.size();) {
    bpk_slash -= 2 + kMeta_Biset * 1.0L / block_sz_;
    size_t lrf_key_sz = nkeys;
    uint64_t delta_sum = batch_sum;
    uint64_t mp_sz = cal_mp_sz(bpk_slash, lrf_key_sz, m);

    uint64_t rho = delta_sum / mp_sz;
    rho <<= 1;

    size_t rho_len = rho == 0 ? max_klen_ : __builtin_clzll(rho);
    double lrf_fpr = rho_len == 0
                         ? key_prefixes_.back()
                         : key_prefixes_.back() - key_prefixes_[rho_len];
    double cur_fpr = lrf_fpr * batch_sum;

    bool all_proteus = false;
    std::vector<bool> indicator(ends.size(), true);
    std::vector<size_t> order = build_seq(begins, ends, interval_sz);
    for (auto iter : order) {
      if (interval_sz[iter] <= 2) {
        all_proteus = true;
        break;
      }
      lrf_key_sz -= interval_sz[iter];
      delta_sum -= ends[iter] - begins[iter];
      mp_sz = cal_mp_sz(bpk_slash, lrf_key_sz, --m);

      if (mp_sz == 0 || delta_sum == 0) {
        all_proteus = true;
        break;
      } else {
        rho = delta_sum / mp_sz;
        rho <<= 1;
        rho_len = rho == 0 ? max_klen_ : __builtin_clzll(rho);
        lrf_fpr = rho_len == 0 ? key_prefixes_.back()
                               : key_prefixes_.back() - key_prefixes_[rho_len];
      }
      double fpr = delta_sum * lrf_fpr + (batch_sum - delta_sum) * proteus_fpr;
      if (fpr >= cur_fpr) {
        lrf_key_sz += interval_sz[iter];
        break;
      } else {
        indicator[iter] = false;
        cur_fpr = fpr;
      }
    }

    if (!all_proteus && cur_fpr < min_fpr) {
      min_fpr = cur_fpr;
      best_indicator = std::move(indicator);
      best_pro_conf = std::move(cur_pro_conf);
      best_idx = set_idx;
      best_lrf_sz = lrf_key_sz;
    }

    size_t pre_idx = set_idx;
    uint64_t cur_threshold = threshold_set_[set_idx];
    while (set_idx < threshold_set_.size() &&
           cur_threshold == threshold_set_[set_idx]) {
      ++set_idx;
    }

    if (set_idx >= threshold_set_.size()) {
      break;
    }

    batch_sum += (set_idx - pre_idx) * cur_threshold;

    shink_index(threshold_set_[set_idx], begins, ends, interval_sz);

    m = ends.size();
    bpk_slash = bpk_ - 1.0L * kCost * m / nkeys;
    if (bpk_slash >= pre_bpk_slash_1 + 1) {
      cur_pro_conf = proteus_model.modeling(keys, bpk_slash);
      proteus_fpr = proteus_model.get_min_fpp();
      pre_bpk_slash_1 = std::floor(bpk_slash);
      continue;
    } else if (bpk_slash - 0.1 >= pre_bpk_slash_2) {
      double mem_budget =
          std::ceil(nkeys * bpk_slash - proteus_model.get_trie_mem());
      proteus_fpr = cal_proteus_fpr(mem_budget, std::get<0>(cur_pro_conf),
                                    std::get<2>(cur_pro_conf));
      pre_bpk_slash_2 = bpk_slash;
    }
  }

  size_t lrf_interval_num = 0;
  if (double lrf_portion = 1.0 * best_lrf_sz / nkeys;
      (best_indicator.size() == 1 && best_indicator[0] == false) ||
      lrf_portion <= 0.1) {
    begins_.clear();
    ends_.clear();
    build_proteus(keys, single_proteus_conf, bpk_);
    learned_rf_ = nullptr;
    return;
  } else if (lrf_portion >= 0.9) {
    lrf_interval_num = ends_.size();
    filter_types_.resize(lrf_interval_num, 0);
    std::iota(filter_types_.begin(), filter_types_.end(), 1);
  } else {
    interval_sz.clear();
    begins.clear();
    ends.clear();
    build_index(threshold_set_[best_idx], keys, begins, ends, interval_sz);

    begins_.clear();
    ends_.clear();

    best_delta_sum = 0;
    size_t cnt = 1;
    for (size_t i = 0; i < ends.size(); ++i) {
      if (best_indicator[i]) {
        filter_types_.emplace_back(cnt++);
        best_delta_sum += ends[i] - begins[i];
      } else {
        if (i > 0 && false == best_indicator[i - 1]) {
          ends_.back() = ends[i];
          continue;
        }
        filter_types_.emplace_back(0);
      }
      begins_.emplace_back(begins[i]);
      ends_.emplace_back(ends[i]);
    }
    lrf_interval_num = cnt;
  }

  std::vector<uint64_t> learned_pos;
  std::vector<uint64_t> proteus_keys;
  get_positions(keys, bpk_, best_delta_sum, best_lrf_sz, lrf_interval_num,
                learned_pos, proteus_keys);

  double last_bits = bpk_ * nkeys;
  if (learned_pos.size() > 0) {
    build_block_lists(learned_pos);
    double used_bits = static_cast<double>(cal_used_bytes()) * 8;
    double bpk_lrf = used_bits / best_lrf_sz;
    if (std::abs(bpk_ - bpk_lrf) >= 0.2) {
      learned_pos.clear();
      block_lists_.clear();

      delete[] bitmap_ptr_;

      std::vector<uint64_t> _;
      get_positions(keys, 2 * bpk_ - bpk_lrf, best_delta_sum, best_lrf_sz,
                    lrf_interval_num, learned_pos, _);
      build_block_lists(learned_pos);
      used_bits = static_cast<double>(cal_used_bytes()) * 8;
    }

    last_bits -= used_bits;
    learned_rf_ =
        new LearnedRF(block_sz_, last_block_sz_, bitmap_sz_, bitmap_ptr_,
                      accumulate_interval_sz_, block_lists_, block_bias_);
  } else {
    learned_rf_ = nullptr;
  }

  // building proteus
  if (proteus_keys.size() > 0) {
    assert(last_bits > 0);
    double p_bpk = last_bits / proteus_keys.size();
    build_proteus(proteus_keys, best_pro_conf, p_bpk);
  } else {
    proteus_ = nullptr;
  }

  printf("\tLearned:%lf\tTradition:\t%lf\n",
         learned_pos.size() * 1.0 / keys.size(),
         proteus_keys.size() * 1.0 / keys.size());
  printf("\tinterval number:%lu\n", ends_.size());
}

void FilterBuilder::segmentation(const std::vector<uint64_t>& keys) {
  size_t nkeys = keys.size();

  std::vector<uint64_t> neigh_dists(nkeys - 1, 0);
  for (size_t i = 1; i < nkeys; ++i) {
    size_t lcp = longestCommonPrefix(keys[i], keys[i - 1], max_klen_);
    ++key_prefixes_[lcp];
    neigh_dists[i - 1] = keys[i] - keys[i - 1];
  }
  std::sort(neigh_dists.begin(), neigh_dists.end());

  key_prefixes_[0] = 1;
  uint64_t mask = UINT64_MAX;
  auto iter = std::upper_bound(neigh_dists.begin(), neigh_dists.end(), mask);
  qk_dists_.emplace_back(neigh_dists.end() - iter);
  for (size_t i = 1; i <= max_klen_; ++i) {
    key_prefixes_[i] += key_prefixes_[i - 1];
    mask >>= 1;
    iter = std::upper_bound(neigh_dists.begin(), iter, mask);
    qk_dists_.emplace_back(neigh_dists.end() - iter);
  }
  M_ = nkeys - nkeys * bpk_ / (kCost + kMeta_LRF);
  uint64_t threshold = neigh_dists[M_];

  while (M_ < neigh_dists.size() && neigh_dists[M_] == threshold) {
    ++M_;
  }

  threshold_set_ = std::move(neigh_dists);
}

void FilterBuilder::build_index(const uint64_t threshold,
                                const std::vector<uint64_t>& keys,
                                std::vector<uint64_t>& begins,
                                std::vector<uint64_t>& ends,
                                std::vector<size_t>& interval_sz) {
  begins.emplace_back(keys[0]);
  size_t back = 0;
  size_t inter_sz;
  for (size_t i = 1; i < keys.size(); ++i) {
    if (keys[i] - keys[i - 1] >= threshold) {
      ends.emplace_back(keys[i - 1]);
      inter_sz = i - back;
      interval_sz.emplace_back(inter_sz);

      back = i;
      begins.emplace_back(keys[i]);
    }
  }

  ends.emplace_back(keys.back());
  inter_sz = keys.size() - back;
  interval_sz.emplace_back(inter_sz);
}

auto FilterBuilder::build_seq(const std::vector<uint64_t>& begins,
                              const std::vector<uint64_t>& ends,
                              const std::vector<size_t>& interval_sz)
    -> std::vector<size_t> {
  std::vector<uint64_t> density(ends.size(), 0);
  for (size_t i = 0; i < ends.size(); ++i) {
    density[i] = interval_sz[i] <= 2
                     ? UINT64_MAX
                     : (ends[i] - begins[i]) / (interval_sz[i] - 2);
  }
  return build_seq(density);
}

auto FilterBuilder::build_seq(const std::vector<uint64_t>& density)
    -> std::vector<size_t> {
  std::vector<size_t> order(density.size());
  std::iota(order.begin(), order.end(), 0);

  std::stable_sort(order.begin(), order.end(), [&density](auto i1, auto i2) {
    return density[i1] < density[i2];
  });
  return order;
}

void FilterBuilder::build_proteus(
    const std::vector<uint64_t>& keys,
    const std::tuple<size_t, size_t, size_t>& conf, double bpk) {
  if (keys.empty()) {
    proteus_ = nullptr;
    return;
  }
  size_t trie_depth = std::get<0>(conf);
  size_t sparse_dense_cutoff = std::get<1>(conf);
  size_t prefix_length = std::get<2>(conf);
  if (trie_depth < 5) {
    trie_depth = 0;
    sparse_dense_cutoff = 0;
  }
  printf(
      "\tTrie Depth: %lu; Sparse-Dense Cutoff (bytes): %lu; BF Prefix "
      "Length: %lu\n",
      trie_depth, sparse_dense_cutoff, prefix_length);
  proteus_ =
      new Proteus(keys, trie_depth, sparse_dense_cutoff, prefix_length, bpk);
}

void FilterBuilder::get_positions(const std::vector<uint64_t>& keys, double bpk,
                                  uint64_t delta_sum, size_t lrf_sz,
                                  size_t lrf_interval_num,
                                  std::vector<uint64_t>& learned_pos,
                                  std::vector<uint64_t>& proteus_keys) {
  accumulate_interval_sz_.clear();
  accumulate_interval_sz_.emplace_back(0);

  size_t nkeys = keys.size();
  double bpk_slash = bpk - 1.0L * kCost * ends_.size() / nkeys -
                     (2 + kMeta_Biset * 1.0L / block_sz_);
  size_t best_mp_sz = cal_mp_sz(bpk_slash, lrf_sz, lrf_interval_num);
  size_t key_idx = 1;
  for (size_t i = 0; i < begins_.size(); ++i) {
    uint64_t begin = begins_[i];
    uint64_t end = ends_[i];
    while (key_idx < nkeys && keys[key_idx] <= begin) {
      ++key_idx;
    }

    if (filter_types_[i] > 0) {
      uint64_t param_1 = end - begin;
      uint64_t param_2 = accumulate_interval_sz_.back();
      uint64_t alpha = std::ceil(1.0L * best_mp_sz * param_1 / delta_sum);
      size_t pre_idx = key_idx;
      while (key_idx < nkeys && keys[key_idx] < end) {
        uint64_t pos = static_cast<uint64_t>(
            static_cast<double>(keys[key_idx++] - begin) * alpha / param_1 +
            param_2);
        learned_pos.emplace_back(pos);
      }
      if (key_idx - pre_idx > 0) {
        accumulate_interval_sz_.emplace_back(param_2 + alpha);
      } else {
        accumulate_interval_sz_.emplace_back(param_2);
      }
    } else {
      while (keys[key_idx] < end) {
        proteus_keys.emplace_back(keys[key_idx++]);
      }
    }
  }
}

}  // namespace oasis_plus