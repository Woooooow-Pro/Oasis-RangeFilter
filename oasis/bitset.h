#pragma once

#include <cassert>
#include <cstdint>
#include <cstring>
#include <vector>

#include "common/type.h"
#include "common/util.h"

namespace oasis {

class BitSet {
 public:
  BitSet(size_t nkeys, BitPos max_range, uint8_t *data);

  ~BitSet() = default;

  auto query(BitPos query) -> bool;
  /* [left, right] */
  auto query(BitPos left, BitPos right) -> bool;

  auto size() const -> size_t;

  static auto size(size_t nkeys, BitPos max_range) -> size_t;
  static auto build(const BitsPos &keys, BitPos max_range)
      -> std::vector<uint8_t>;

 private:
  inline static auto get_lower_bit_length(BitPos nkeys, BitPos max_range)
      -> BitPos;

  inline static auto bit_array_bsize(BitPos nkeys, BitPos max_range,
                                     BitPos lower_bit_len) -> size_t;

  inline auto get_lower_bit_length(BitPos max_range) -> BitPos;

  inline auto bit_array_bsize(BitPos max_range) -> size_t;

  inline auto find_next_set_bit(size_t idx) -> size_t;

 private:
  BitPos lower_bit_len_;
  size_t b_size_;
  size_t nkeys_;

  // need to be serialized
  uint8_t *data_;
};

BitSet::BitSet(size_t nkeys, BitPos max_range, uint8_t *data) : nkeys_(nkeys) {
  lower_bit_len_ = get_lower_bit_length(max_range);
  b_size_ = bit_array_bsize(max_range);
  data_ = data;
}

auto BitSet::query(BitPos query) -> bool {
  size_t up_idx = lower_bit_len_ * nkeys_;
  BitPos up_query = query >> lower_bit_len_;

  // upper bit search
  size_t zero_cnt = 0;
  size_t range = 0;
  BitPos upper = 0;
  while (upper <= up_query && up_idx < b_size_) {
    if (data_[up_idx >> 3] & (1U << (up_idx & 0x7ULL))) {
      ++up_idx;
      ++upper;
    } else {
      size_t next_idx = find_next_set_bit(up_idx);
      size_t delta = next_idx - up_idx;
      if (upper == up_query) {
        range = delta;
        break;
      }
      zero_cnt += delta;
      up_idx = next_idx;
    }
  }

  if (range == 0) {
    return false;
  } else if (lower_bit_len_ == 0) {
    return true;
  }

  // lower bit search
  BitPos low_query = query & ((1 << lower_bit_len_) - 1);
  size_t low_idx = zero_cnt * lower_bit_len_;
  while (low_idx < (zero_cnt + range) * lower_bit_len_) {
    BitPos low_key = 0;
    BitPos bit_cnt = 0;
    while (bit_cnt < lower_bit_len_) {
      BitPos byte_bias = low_idx & 0x7ULL;
      BitPos delta = std::min((8 - byte_bias), lower_bit_len_ - bit_cnt);
      BitPos tmp_mask = (1 << delta) - 1;
      BitPos cur_data = (data_[low_idx >> 3] >> byte_bias) & tmp_mask;
      low_key |= cur_data << bit_cnt;
      bit_cnt += delta;
      low_idx += delta;
    }
    if (low_key >= low_query) {
      return low_key == low_query;
    }
  }
  return false;
}

auto BitSet::query(BitPos left, BitPos right) -> bool {
  size_t up_idx = lower_bit_len_ * nkeys_;
  BitPos up_left = left >> lower_bit_len_;
  BitPos up_right = right >> lower_bit_len_;

  // upper bit search
  size_t zero_cnt = 0;
  size_t range = 0;
  BitPos upper = 0;
  std::vector<std::pair<BitPos, size_t>> upper_keys;
  while (upper <= up_right && up_idx < b_size_) {
    if (data_[up_idx >> 3] & (1U << (up_idx & 0x7ULL))) {
      ++up_idx;
      ++upper;
    } else {
      size_t next_idx = find_next_set_bit(up_idx);
      size_t delta = next_idx - up_idx;
      if (upper >= up_left && upper <= up_right) {
        if (range > 0 && upper < up_right) {
          return true;
        }

        upper_keys.emplace_back(upper << lower_bit_len_, delta);
        range += delta;
      }
      zero_cnt += range > 0 ? 0 : delta;
      up_idx = next_idx;
    }
  }

  if (range == 0) {
    return false;
  } else if (lower_bit_len_ == 0) {
    return true;
  }

  size_t key_begin_idx = zero_cnt;
  size_t low_idx = key_begin_idx * lower_bit_len_;
  for (const auto &[up_key, n] : upper_keys) {
    key_begin_idx += n;
    while (low_idx < key_begin_idx * lower_bit_len_) {
      BitPos low_key = 0;
      BitPos bit_cnt = 0;
      while (bit_cnt < lower_bit_len_) {
        BitPos byte_bias = low_idx & 0x7ULL;
        BitPos delta = std::min((8 - byte_bias), lower_bit_len_ - bit_cnt);
        BitPos tmp_mask = (1 << delta) - 1;
        BitPos cur_data = (data_[low_idx >> 3] >> byte_bias) & tmp_mask;
        low_key |= cur_data << bit_cnt;
        bit_cnt += delta;
        low_idx += delta;
      }
      BitPos key = low_key | up_key;
      if (key >= left && key <= right) {
        return true;
      }
      if (key > right) {
        return false;
      }
    }
  }
  return false;
}

auto BitSet::size() const -> size_t {
  return align_bit2byte(b_size_);  // data_ size
}

auto BitSet::size(size_t nkeys, BitPos max_range) -> size_t {
  BitPos lower_bit_len = get_lower_bit_length(nkeys, max_range);
  size_t b_size_ = bit_array_bsize(nkeys, max_range, lower_bit_len);
  return align_bit2byte(b_size_);
}

/** Helping Method */
auto BitSet::get_lower_bit_length(BitPos nkeys, BitPos max_range) -> BitPos {
  return __builtin_clzl(nkeys) > __builtin_clzl(max_range)
             ? __builtin_clzl(nkeys) - __builtin_clzl(max_range)
             : 0;
}

auto BitSet::get_lower_bit_length(BitPos max_range) -> BitPos {
  return get_lower_bit_length(nkeys_, max_range);
}

auto BitSet::bit_array_bsize(BitPos nkeys, BitPos max_range,
                             BitPos lower_bit_len) -> size_t {
  BitPos total_quotient = max_range >> lower_bit_len;
  return nkeys * (1U + lower_bit_len) + total_quotient + 1;
}

auto BitSet::bit_array_bsize(BitPos max_range) -> size_t {
  return bit_array_bsize(nkeys_, max_range, lower_bit_len_);
}

auto BitSet::find_next_set_bit(size_t idx) -> size_t {
  if (!(data_[idx >> 3] >> (idx & 0x7U))) {
    idx = (idx & (~static_cast<size_t>(0x7U))) + 8U;
  }

  while (data_[idx >> 3] == 0) {
    idx += 8U;
  }

  idx += __builtin_ctz(data_[idx >> 3] >> (idx & 0x7U));
  return idx;
}

auto BitSet::build(const BitsPos &keys, BitPos max_range)
    -> std::vector<uint8_t> {
  size_t nkeys = keys.size();

  BitPos lower_bit_len = get_lower_bit_length(nkeys, max_range);
  size_t b_size_ = bit_array_bsize(nkeys, max_range, lower_bit_len);
  size_t data_size = align_bit2byte(b_size_);
  std::vector<uint8_t> data(data_size, 0);

  size_t low_idx = 0;
  size_t up_idx = lower_bit_len * nkeys;

  BitPos pre_upper = 0;
  for (size_t i = 0; i < keys.size(); ++i) {
    BitPos lower_key = keys[i] & ((1 << lower_bit_len) - 1);
    BitPos bit_cnt = 0;
    while (bit_cnt < lower_bit_len) {
      BitPos byte_bias = low_idx & 0x7ULL;
      BitPos delta = std::min((8 - byte_bias), lower_bit_len - bit_cnt);
      BitPos tmp_mask = (1 << delta) - 1;
      data[low_idx >> 3] |= (lower_key & tmp_mask) << byte_bias;
      low_idx += delta;
      bit_cnt += delta;
      lower_key >>= delta;
    }

    while (pre_upper < (keys[i] >> lower_bit_len)) {
      data[up_idx >> 3] |= 1 << (up_idx & 0x7ULL);
      ++pre_upper;
      ++up_idx;
    }
    ++up_idx;
  }
  data[up_idx >> 3] |= 1 << (up_idx & 0x7ULL);

  return data;
}

}  // namespace oasis