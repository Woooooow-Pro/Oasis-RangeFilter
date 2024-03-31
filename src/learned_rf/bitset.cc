#include "learned_rf/bitset.h"

#include <cassert>
#include <cstring>

#include "util.h"

namespace oasis_plus {
BitSet::BitSet(size_t nkeys, uint64_t max_range, uint8_t *data)
    : nkeys_(nkeys) {
  lower_bit_len_ = get_lower_bit_length(max_range);
  b_size_ = bit_array_bsize(max_range);
  data_ = data;
}

auto BitSet::query(uint64_t query) -> bool {
  size_t up_idx = lower_bit_len_ * nkeys_;
  uint64_t up_query = query >> lower_bit_len_;

  // upper bit search
  size_t zero_cnt = 0;
  size_t range = 0;
  uint64_t upper = 0;
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
  uint64_t low_query = query & ((1 << lower_bit_len_) - 1);
  size_t low_idx = zero_cnt * lower_bit_len_;
  while (low_idx < (zero_cnt + range) * lower_bit_len_) {
    uint64_t low_key = 0;
    uint64_t bit_cnt = 0;
    while (bit_cnt < lower_bit_len_) {
      uint64_t byte_bias = low_idx & 0x7ULL;
      uint64_t delta = std::min((8 - byte_bias), lower_bit_len_ - bit_cnt);
      uint64_t tmp_mask = (1 << delta) - 1;
      uint64_t cur_data = (data_[low_idx >> 3] >> byte_bias) & tmp_mask;
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

auto BitSet::query(uint64_t left, uint64_t right) -> bool {
  size_t up_idx = lower_bit_len_ * nkeys_;
  uint64_t up_left = left >> lower_bit_len_;
  uint64_t up_right = right >> lower_bit_len_;

  // upper bit search
  size_t zero_cnt = 0;
  size_t range = 0;
  uint64_t upper = 0;
  std::vector<std::pair<uint64_t, size_t>> upper_keys;
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
      uint64_t low_key = 0;
      uint64_t bit_cnt = 0;
      while (bit_cnt < lower_bit_len_) {
        uint64_t byte_bias = low_idx & 0x7ULL;
        uint64_t delta = std::min((8 - byte_bias), lower_bit_len_ - bit_cnt);
        uint64_t tmp_mask = (1 << delta) - 1;
        uint64_t cur_data = (data_[low_idx >> 3] >> byte_bias) & tmp_mask;
        low_key |= cur_data << bit_cnt;
        bit_cnt += delta;
        low_idx += delta;
      }
      uint64_t key = low_key | up_key;
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
  return BIT2BYTE(b_size_);  // data_ size
}

auto BitSet::size(size_t nkeys, uint64_t max_range) -> size_t {
  uint64_t lower_bit_len = get_lower_bit_length(nkeys, max_range);
  size_t b_size_ = bit_array_bsize(nkeys, max_range, lower_bit_len);
  return BIT2BYTE(b_size_);
}

/** Helping Method */
auto BitSet::get_lower_bit_length(uint64_t nkeys, uint64_t max_range)
    -> uint64_t {
  return __builtin_clzl(nkeys) > __builtin_clzl(max_range)
             ? __builtin_clzl(nkeys) - __builtin_clzl(max_range)
             : 0;
}

auto BitSet::get_lower_bit_length(uint64_t max_range) -> uint64_t {
  return get_lower_bit_length(nkeys_, max_range);
}

auto BitSet::bit_array_bsize(uint64_t nkeys, uint64_t max_range,
                             uint64_t lower_bit_len) -> size_t {
  uint64_t total_quotient = max_range >> lower_bit_len;
  return nkeys * (1U + lower_bit_len) + total_quotient + 1;
}

auto BitSet::bit_array_bsize(uint64_t max_range) -> size_t {
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

auto BitSet::build(const std::vector<uint64_t> &keys, uint64_t max_range)
    -> std::vector<uint8_t> {
  size_t nkeys = keys.size();

  uint64_t lower_bit_len = get_lower_bit_length(nkeys, max_range);
  size_t b_size_ = bit_array_bsize(nkeys, max_range, lower_bit_len);
  size_t data_size = BIT2BYTE(b_size_);
  std::vector<uint8_t> data(data_size, 0);

  size_t low_idx = 0;
  size_t up_idx = lower_bit_len * nkeys;

  uint64_t pre_upper = 0;
  for (size_t i = 0; i < keys.size(); ++i) {
    uint64_t lower_key = keys[i] & ((1 << lower_bit_len) - 1);
    uint64_t bit_cnt = 0;
    while (bit_cnt < lower_bit_len) {
      uint64_t byte_bias = low_idx & 0x7ULL;
      uint64_t delta = std::min((8 - byte_bias), lower_bit_len - bit_cnt);
      uint64_t tmp_mask = (1 << delta) - 1;
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
}  // namespace oasis_plus