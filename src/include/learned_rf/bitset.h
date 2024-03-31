#pragma once

#include <cstdint>
#include <vector>

namespace oasis_plus {

class BitSet {
 public:
  ~BitSet() = default;
  BitSet(size_t nkeys, uint64_t max_range, uint8_t *data);


  auto query(uint64_t query) -> bool;
  /* [left, right] */
  auto query(uint64_t left, uint64_t right) -> bool;

  auto size() const -> size_t;

  static auto size(size_t nkeys, uint64_t max_range) -> size_t;
  static auto build(const std::vector<uint64_t> &keys, uint64_t max_range)
      -> std::vector<uint8_t>;

 private:
  inline static auto get_lower_bit_length(uint64_t nkeys, uint64_t max_range)
      -> uint64_t;

  inline static auto bit_array_bsize(uint64_t nkeys, uint64_t max_range,
                                     uint64_t lower_bit_len) -> size_t;

  inline auto get_lower_bit_length(uint64_t max_range) -> uint64_t;

  inline auto bit_array_bsize(uint64_t max_range) -> size_t;

  inline auto find_next_set_bit(size_t idx) -> size_t;

 private:
  uint64_t lower_bit_len_;
  size_t b_size_;
  size_t nkeys_;

  uint8_t *data_;
};

}  // namespace oasis_plus
