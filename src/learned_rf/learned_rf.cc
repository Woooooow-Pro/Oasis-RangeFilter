#include "learned_rf/learned_rf.h"

#include <algorithm>
#include <cassert>
#include <cstring>

#include "util.h"

namespace oasis_plus {

LearnedRF::~LearnedRF() { delete[] bitmap_ptr_; }

auto LearnedRF::query(uint64_t key, size_t interval_idx, uint64_t low,
                      uint64_t up) -> bool {
  auto params = get_params(low, up, interval_idx);
  if (params.first == 0) {
    return false;
  }
  uint64_t location = get_location(key - low, interval_idx, params);

  if (location < block_bias_[0] || location > block_bias_.back()) {
    return false;
  }
  auto iter =
      std::upper_bound(block_bias_.begin(), block_bias_.end(), location) - 1;
  if (*iter == location) {
    return true;
  }

  size_t block_idx = iter - block_bias_.begin();
  return block_lists_[block_idx].query(location - *iter);
}

auto LearnedRF::query(uint64_t l_key, uint64_t r_key, size_t interval_idx,
                      uint64_t low, uint64_t up) -> bool {
  auto params = get_params(low, up, interval_idx);
  if (params.first == 0) {
    return false;
  }
  uint64_t l_location = get_location(l_key - low, interval_idx, params);
  uint64_t r_location = get_location(r_key - low, interval_idx, params);

  if (r_location < block_bias_[0] || l_location > block_bias_.back()) {
    return false;
  }

  auto iter =
      std::upper_bound(block_bias_.begin(), block_bias_.end(), r_location);
  if (iter == block_bias_.end() || *(--iter) == r_location ||
      l_location <= *iter) {
    return true;
  }

  size_t block_idx = iter - block_bias_.begin();
  return block_lists_[block_idx].query(l_location - *iter, r_location - *iter);
}

auto LearnedRF::serialize() const -> std::pair<uint8_t *, size_t> {
  size_t meta_sz = sizeof(uint16_t)   /** block_sz_ */
                   + sizeof(uint16_t) /** last_block_sz_ */
                   + sizeof(size_t)   /** record size of compressed bitmap */
                   + sizeof(uint32_t) /** # accumulate_interval_sz_ */
                   + sizeof(uint32_t) /** # block_lists_ entries */
                   + bitmap_sz_;      /** size of compressed bitmap  */
  size_align(meta_sz);

  uint32_t nintervals = accumulate_interval_sz_.size() - 1;
  uint32_t nbatches = block_lists_.size();
  size_t size = meta_sz +
                sizeof(uint64_t) *
                    nintervals /** size of accumulate interval size array */
                + sizeof(uint64_t) * (nbatches + 1); /** size of block_bias_ */

  uint8_t *ser = new uint8_t[size];
  uint8_t *pos = ser;

  memcpy(pos, &block_sz_, sizeof(uint16_t));
  pos += sizeof(uint16_t);

  memcpy(pos, &last_block_sz_, sizeof(uint16_t));
  pos += sizeof(uint16_t);

  memcpy(pos, &bitmap_sz_, sizeof(size_t));
  pos += sizeof(size_t);

  memcpy(pos, &nintervals, sizeof(uint32_t));
  pos += sizeof(uint32_t);

  memcpy(pos, &nbatches, sizeof(uint32_t));
  pos += sizeof(uint32_t);

  memcpy(pos, bitmap_ptr_, bitmap_sz_);
  pos += bitmap_sz_;

  align(pos);

  memcpy(pos, accumulate_interval_sz_.data() + 1,
         sizeof(uint64_t) * nintervals);
  pos += sizeof(uint64_t) * nintervals;

  memcpy(pos, block_bias_.data(), sizeof(uint64_t) * (nbatches + 1));
  pos += sizeof(uint64_t) * (nbatches + 1);

  return {ser, size};
}

auto LearnedRF::deserialize(uint8_t *ser) -> std::pair<LearnedRF *, size_t> {
  uint8_t *pos = ser;

  uint16_t block_sz;
  memcpy(&block_sz, pos, sizeof(uint16_t));
  pos += sizeof(uint16_t);

  uint16_t last_block_sz;
  memcpy(&last_block_sz, pos, sizeof(uint16_t));
  pos += sizeof(uint16_t);

  size_t bitmap_sz;
  memcpy(&bitmap_sz, pos, sizeof(size_t));
  pos += sizeof(size_t);

  uint32_t ninterval;
  memcpy(&ninterval, pos, sizeof(uint32_t));
  pos += sizeof(uint32_t);

  uint32_t nbatches;
  memcpy(&nbatches, pos, sizeof(uint32_t));
  pos += sizeof(uint32_t);

  uint8_t *bitmap_ptr = new uint8_t[bitmap_sz];
  memcpy(bitmap_ptr, pos, bitmap_sz);
  pos += bitmap_sz;

  align(pos);

  std::vector<uint64_t> accumulate_interval_sz(ninterval + 1, 0);
  memcpy(accumulate_interval_sz.data() + 1, pos, sizeof(uint64_t) * ninterval);
  pos += sizeof(uint64_t) * ninterval;

  std::vector<uint64_t> block_bias(nbatches + 1, 0);
  memcpy(block_bias.data(), pos, sizeof(uint64_t) * (nbatches + 1));
  pos += sizeof(uint64_t) * (nbatches + 1);

  std::vector<BitSet> block_lists;
  uint8_t *data = bitmap_ptr;
  for (size_t i = 0; i < nbatches - 1; ++i) {
    block_lists.emplace_back(block_sz, block_bias[i + 1] - block_bias[i], data);
    data += block_lists.back().size();
  }
  block_lists.emplace_back(
      last_block_sz, block_bias[nbatches] - block_bias[nbatches - 1], data);

  return {new LearnedRF(block_sz, last_block_sz, bitmap_sz, bitmap_ptr,
                        accumulate_interval_sz, block_lists, block_bias),
          pos - ser};
}

auto LearnedRF::size() const -> size_t {
  size_t meta_sz = sizeof(uint16_t)   /** block_sz_ */
                   + sizeof(uint16_t) /** last_block_sz_ */
                   + sizeof(size_t)   /** record size of compressed bitmap */
                   + sizeof(uint32_t) /** # accumulate_interval_sz_ */
                   + sizeof(uint32_t) /** # block_lists_ entries */
                   + bitmap_sz_;      /** size of compressed bitmap  */
  size_align(meta_sz);

  uint32_t nintervals = accumulate_interval_sz_.size() - 1;
  uint32_t nbatches = block_lists_.size();
  return meta_sz +
         sizeof(uint64_t) *
             nintervals /** size of accumulate interval size array */
         + sizeof(uint64_t) * (nbatches + 1);
}
}  // namespace oasis_plus
