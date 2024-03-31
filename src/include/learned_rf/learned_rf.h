#pragma once

#include <cstdint>
#include <vector>

#include "learned_rf/bitset.h"

namespace oasis_plus {
class LearnedRF {
 public:
  LearnedRF(uint16_t block_sz, uint16_t last_block_sz, size_t bitmap_sz,
            uint8_t *bitmap_ptr, std::vector<uint64_t> &accumulate_interval_sz,
            std::vector<BitSet> &block_lists, std::vector<uint64_t> block_bias)
      : block_sz_(block_sz),
        last_block_sz_(last_block_sz),
        bitmap_sz_(bitmap_sz),
        bitmap_ptr_(bitmap_ptr),
        accumulate_interval_sz_(std::move(accumulate_interval_sz)),
        block_lists_(std::move(block_lists)),
        block_bias_(std::move(block_bias)) {}
  ~LearnedRF();

  auto query(uint64_t key, size_t interval_idx, uint64_t low, uint64_t up)
      -> bool;
  auto query(uint64_t l_key, uint64_t r_key, size_t interval_idx, uint64_t low,
             uint64_t up) -> bool;

  auto serialize() const -> std::pair<uint8_t *, size_t>;
  static auto deserialize(uint8_t *ser) -> std::pair<LearnedRF *, size_t>;

  auto size() const -> size_t;

 private:
  /** Helping Methods */
  auto get_params(uint64_t low, uint64_t up, size_t interval_idx)
      -> std::pair<uint64_t, uint64_t> {
    return {static_cast<uint64_t>(accumulate_interval_sz_[interval_idx] -
                                  accumulate_interval_sz_[interval_idx - 1]),
            up - low};
  }

  auto get_location(double delta_key, size_t interval_idx,
                    const std::pair<uint64_t, uint64_t> &params) -> uint64_t {
    return delta_key * params.first / params.second +
           accumulate_interval_sz_[interval_idx - 1];
  }

 private:
  uint16_t block_sz_;
  uint16_t last_block_sz_;
  size_t bitmap_sz_ = 0;

  uint8_t *bitmap_ptr_ = nullptr;

  std::vector<uint64_t> accumulate_interval_sz_;

  std::vector<BitSet> block_lists_;
  std::vector<uint64_t> block_bias_;
};
}  // namespace oasis_plus
