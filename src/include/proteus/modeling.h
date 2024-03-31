#pragma once

#include <array>
#include <cassert>
#include <cstring>
#include <vector>

#include "proteus/config.h"

namespace oasis_plus {

typedef std::array<std::pair<size_t, size_t>, 64> bin_array;

class ProteusModeling {
 public:
  ProteusModeling(const size_t max_klen, const size_t max_qlen = 10)
      : max_klen_(max_klen),
        max_qlen_(max_qlen),
        key_prefixes_(max_klen_ + 1, 0),
        bf_mem_(max_klen_ + 1, 0.0L) {}

  template <typename T>
  auto modeling(const std::vector<T>& keys,
                const std::vector<std::pair<T, T>>& sample_queries,
                const double bits_per_key,
                std::vector<size_t>* sparse_dense_cutoffs = nullptr)
      -> std::tuple<size_t, size_t, size_t>;

  template <typename T>
  auto modeling(const std::vector<T>& keys, const double bits_per_key,
                std::vector<size_t>* sparse_dense_cutoffs = nullptr)
      -> std::tuple<size_t, size_t, size_t>;

 public:
  void set_key_prefixes(const std::vector<size_t>& key_prefixes) {
    assert(key_prefixes_.size() == key_prefixes.size());
    memcpy(key_prefixes_.data(), key_prefixes.data(),
           key_prefixes.size() * sizeof(size_t));
  }

  void set_qk_dists(const std::vector<size_t>& qk_dists) {
    qk_dists_.resize(qk_dists.size());
    memcpy(qk_dists_.data(), qk_dists.data(), qk_dists.size() * sizeof(size_t));
  }

  auto get_bf_mem(size_t trie_depth) const -> long double {
    return bf_mem_[trie_depth];
  }

  auto get_trie_mem() const -> long double { return trie_mem_; }

  auto get_min_fpp() const -> long double { return min_fpp_; }

 private:
  size_t max_klen_;
  size_t max_qlen_;
  // Number of unique key prefixes for every prefix length
  std::vector<size_t> key_prefixes_;
  std::vector<size_t> qk_dists_;
  // Bloom filter memory available for every trie depth
  std::vector<long double> bf_mem_;
  long double trie_mem_;

  // set when the model find the best config
  long double min_fpp_;
};

}  // namespace oasis_plus
