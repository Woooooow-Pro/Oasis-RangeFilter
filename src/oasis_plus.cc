#include "oasis_plus.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <numeric>

#include "filter_builder.h"
#include "util.h"

namespace oasis_plus {

OasisPlus::OasisPlus(double bpk, uint32_t block_size,
                     const std::vector<uint64_t> &keys, const size_t max_qlen) {
  FilterBuilder filter_builder(bpk, block_size, max_qlen);
  filter_builder.build(keys);
  filter_types_ = filter_builder.get_filter_types();
  begins_ = filter_builder.get_begins();
  ends_ = filter_builder.get_ends();

  learned_rf_ = filter_builder.get_learned_rf();
  proteus_ = filter_builder.get_proteus();
}

auto OasisPlus::query(uint64_t key) -> bool {
  if (learned_rf_ == nullptr) {
    // single proteus
    return proteus_->Query(key);
  }

  if (key < begins_[0] || key > ends_.back()) {
    return false;
  }

  size_t idx =
      std::distance(begins_.begin(),
                    std::upper_bound(begins_.begin(), begins_.end(), key)) -
      1;

  if (ends_[idx] < key) {
    return false;
  }
  if (begins_[idx] == key || ends_[idx] == key) {
    return true;
  }

  if (size_t interval_idx = proteus_ == nullptr ? idx + 1 : filter_types_[idx];
      interval_idx != 0) {
    return learned_rf_->query(key, interval_idx, begins_[idx], ends_[idx]);
  }
  // return Proteus query
  return proteus_->Query(key);
}

auto OasisPlus::query(uint64_t l_key, uint64_t r_key) -> bool {
  assert(l_key < r_key);

  if (learned_rf_ == nullptr) {
    // single proteus
    return proteus_->Query(l_key, r_key + 1);
  }

  if (l_key > ends_.back() || r_key < begins_[0]) {
    return false;
  }

  size_t idx = std::upper_bound(begins_.begin(), begins_.end(), l_key) -
               begins_.begin() - 1;

  if (r_key < begins_[idx + 1] && l_key > ends_[idx]) {
    return false;
  }
  if (!(l_key > begins_[idx] && r_key < ends_[idx])) {
    return true;
  }

  if (size_t interval_idx = proteus_ == nullptr ? idx + 1 : filter_types_[idx];
      interval_idx != 0) {
    return learned_rf_->query(l_key, r_key, interval_idx, begins_[idx],
                              ends_[idx]);
  }

  // calling the Proteus
  return proteus_ != nullptr && proteus_->Query(l_key, r_key + 1);
}

auto OasisPlus::serialize() const -> std::pair<uint8_t *, size_t> {
  uint8_t filter_bitmap = 0;

  size_t learned_rf_sz = 0;
  uint8_t *learned_ptr = nullptr;
  if (learned_rf_ != nullptr) {
    filter_bitmap |= LEARNED_EXIST_BIT;
    auto data = learned_rf_->serialize();
    learned_ptr = data.first;
    learned_rf_sz = data.second;
  }

  size_t proteus_sz = 0;
  uint8_t *proteus_ptr = nullptr;
  if (proteus_ != nullptr) {
    filter_bitmap |= PROTEUS_EXIST_BIT;
    auto data = proteus_->serialize();
    proteus_ptr = data.first;
    proteus_sz = data.second;
  }

  size_t meta_sz = sizeof(size_t)     /* # intervals */
                   + sizeof(uint8_t); /* fitler bitmap size*/

  size_t interval_num = begins_.size();
  size_t bitmap_sz = 0;
  if (filter_bitmap == (LEARNED_EXIST_BIT | PROTEUS_EXIST_BIT)) {
    bitmap_sz = BIT2BYTE(interval_num);
    meta_sz += bitmap_sz;
  }
  size_align(meta_sz);

  size_t size = meta_sz + learned_rf_sz +
                2 * sizeof(uint64_t) * interval_num; /* indices size */
  size_align(size);
  size += proteus_sz;

  uint8_t *ser = new uint8_t[size];
  uint8_t *pos = ser;

  memcpy(pos, &interval_num, sizeof(size_t));
  pos += sizeof(size_t);

  memcpy(pos, &filter_bitmap, sizeof(uint8_t));
  pos += sizeof(uint8_t);

  if (filter_bitmap == (LEARNED_EXIST_BIT | PROTEUS_EXIST_BIT)) {
    std::vector<uint8_t> filter_type_bitmap(bitmap_sz, 0);
    for (size_t i = 0; i < interval_num; ++i) {
      if (filter_types_[i] != 0) {
        filter_type_bitmap[i >> 3] |= 1ULL << (i & 7ULL);
      }
    }
    memcpy(pos, filter_type_bitmap.data(), bitmap_sz);
    pos += bitmap_sz;
  }

  align(pos);

  memcpy(pos, begins_.data(), begins_.size() * sizeof(uint64_t));
  pos += begins_.size() * sizeof(uint64_t);

  memcpy(pos, ends_.data(), ends_.size() * sizeof(uint64_t));
  pos += ends_.size() * sizeof(uint64_t);

  if (learned_ptr != nullptr) {
    std::memcpy(pos, learned_ptr, learned_rf_sz);
    pos += learned_rf_sz;
    delete[] learned_ptr;
    align(pos);
  }

  if (proteus_ptr != nullptr) {
    memcpy(pos, proteus_ptr, proteus_sz);
    pos += proteus_sz;
    delete[] proteus_ptr;
  }

  return {ser, size};
}

auto OasisPlus::deserialize(uint8_t *ser) -> OasisPlus * {
  size_t interval_num = 0;
  memcpy(&interval_num, ser, sizeof(size_t));
  ser += sizeof(size_t);

  uint8_t bitmap = ser[0];
  ser += sizeof(uint8_t);

  std::vector<uint32_t> filter_types(interval_num, 0);
  if (bitmap == (LEARNED_EXIST_BIT | PROTEUS_EXIST_BIT)) {
    uint32_t cnt = 0;
    for (size_t i = 0; i < interval_num; ++i) {
      if (ser[i >> 3] & (1ULL << (i & 7ULL))) {
        filter_types[i] = ++cnt;
      }
    }

    ser += BIT2BYTE(interval_num);
  }

  align(ser);

  std::vector<uint64_t> begins(interval_num);
  memcpy(begins.data(), ser, interval_num * sizeof(uint64_t));
  ser += interval_num * sizeof(uint64_t);

  std::vector<uint64_t> ends(interval_num);
  memcpy(ends.data(), ser, interval_num * sizeof(uint64_t));
  ser += interval_num * sizeof(uint64_t);

  LearnedRF *learned_rf = nullptr;
  if (bitmap & LEARNED_EXIST_BIT) {
    auto data = LearnedRF::deserialize(ser);
    learned_rf = data.first;
    ser += data.second;
    align(ser);
  }

  Proteus *proteus = nullptr;
  if (bitmap & PROTEUS_EXIST_BIT) {
    auto data = Proteus::deSerialize(reinterpret_cast<char *>(ser));
    proteus = data.first;
    ser += data.second;
  }
  return {new OasisPlus(begins, ends, filter_types, learned_rf, proteus)};
}

auto OasisPlus::size() const -> size_t {
  size_t size = sizeof(size_t)                          /* # intervals */
                + 2 * sizeof(uint64_t) * begins_.size() /* indices size */
                + sizeof(uint8_t) /* fitler bitmap size*/;

  if (proteus_ != nullptr && learned_rf_ != nullptr) {
    size += BIT2BYTE(filter_types_.size());
  }

  if (learned_rf_ != nullptr) {
    size += learned_rf_->size();
  }
  if (proteus_ != nullptr) {
    auto data = proteus_->serialize();
    size += data.second;
    delete[] data.first;
  }
  return size;
}

}  // namespace oasis_plus