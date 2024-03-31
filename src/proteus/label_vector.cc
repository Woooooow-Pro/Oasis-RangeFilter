#include "proteus/label_vector.h"

#include <cstring>

namespace oasis_plus {
LabelVector::LabelVector(
    const std::vector<std::vector<label_t> >& labels_per_level,
    const level_t start_level, level_t end_level /* non-inclusive */) {
  num_bytes_ = 1;
  for (level_t level = start_level; level < end_level; level++)
    num_bytes_ += labels_per_level[level].size();

  position_t alloc_bytes = num_bytes_ * (num_bytes_ / kWordSize + 1);
  labels_ = new label_t[alloc_bytes];
  for (position_t i = 0; i < alloc_bytes; i++) labels_[i] = 0;

  position_t pos = 0;
  for (level_t level = start_level; level < end_level; level++) {
    for (position_t idx = 0; idx < labels_per_level[level].size(); idx++) {
      labels_[pos] = labels_per_level[level][idx];
      pos++;
    }
  }
}

auto LabelVector::serializedSize() const -> position_t {
  position_t size = sizeof(num_bytes_) + num_bytes_;
  sizeAlign(size);
  return size;
}

auto LabelVector::size() const -> position_t {
  return (sizeof(LabelVector) + num_bytes_);
}

auto LabelVector::read(const position_t pos) const -> label_t {
  return labels_[pos];
}

auto LabelVector::operator[](const position_t pos) const -> label_t {
  return labels_[pos];
}

auto LabelVector::search(const label_t target, position_t& pos,
                         position_t search_len) const -> bool {
  if (search_len < 3) return linearSearch(target, pos, search_len);
  if (search_len < 12)
    return binarySearch(target, pos, search_len);
  else
    return simdSearch(target, pos, search_len);
}

auto LabelVector::searchGreaterThan(const label_t target, position_t& pos,
                                    position_t search_len) const -> bool {
  if (search_len < 3)
    return linearSearchGreaterThan(target, pos, search_len);
  else
    return binarySearchGreaterThan(target, pos, search_len);
}

auto LabelVector::binarySearch(const label_t target, position_t& pos,
                               const position_t search_len) const -> bool {
  position_t l = pos;
  position_t r = pos + search_len;
  while (l < r) {
    position_t m = (l + r) >> 1;
    if (target < labels_[m]) {
      r = m;
    } else if (target == labels_[m]) {
      pos = m;
      return true;
    } else {
      l = m + 1;
    }
  }
  return false;
}

auto LabelVector::simdSearch(const label_t target, position_t& pos,
                             const position_t search_len) const -> bool {
  position_t num_labels_searched = 0;
  position_t num_labels_left = search_len;
  while ((num_labels_left >> 4) > 0) {
    label_t* start_ptr = labels_ + pos + num_labels_searched;
    __m128i cmp =
        _mm_cmpeq_epi8(_mm_set1_epi8(target),
                       _mm_loadu_si128(reinterpret_cast<__m128i*>(start_ptr)));
    unsigned check_bits = _mm_movemask_epi8(cmp);
    if (check_bits) {
      pos += (num_labels_searched + __builtin_ctz(check_bits));
      return true;
    }
    num_labels_searched += 16;
    num_labels_left -= 16;
  }

  if (num_labels_left > 0) {
    label_t* start_ptr = labels_ + pos + num_labels_searched;
    __m128i cmp =
        _mm_cmpeq_epi8(_mm_set1_epi8(target),
                       _mm_loadu_si128(reinterpret_cast<__m128i*>(start_ptr)));
    unsigned leftover_bits_mask = (1 << num_labels_left) - 1;
    unsigned check_bits = _mm_movemask_epi8(cmp) & leftover_bits_mask;
    if (check_bits) {
      pos += (num_labels_searched + __builtin_ctz(check_bits));
      return true;
    }
  }

  return false;
}

auto LabelVector::linearSearch(const label_t target, position_t& pos,
                               const position_t search_len) const -> bool {
  for (position_t i = 0; i < search_len; i++) {
    if (target == labels_[pos + i]) {
      pos += i;
      return true;
    }
  }
  return false;
}

auto LabelVector::binarySearchGreaterThan(const label_t target, position_t& pos,
                                          const position_t search_len) const
    -> bool {
  position_t l = pos;
  position_t r = pos + search_len;
  while (l < r) {
    position_t m = (l + r) >> 1;
    if (target < labels_[m]) {
      r = m;
    } else if (target == labels_[m]) {
      if (m < pos + search_len - 1) {
        pos = m + 1;
        return true;
      }
      return false;
    } else {
      l = m + 1;
    }
  }

  if (l < pos + search_len) {
    pos = l;
    return true;
  }
  return false;
}

auto LabelVector::linearSearchGreaterThan(const label_t target, position_t& pos,
                                          const position_t search_len) const
    -> bool {
  for (position_t i = 0; i < search_len; i++) {
    if (labels_[pos + i] > target) {
      pos += i;
      return true;
    }
  }
  return false;
}

void LabelVector::serialize(char*& dst) const {
  memcpy(dst, &num_bytes_, sizeof(num_bytes_));
  dst += sizeof(num_bytes_);
  memcpy(dst, labels_, num_bytes_);
  dst += num_bytes_;
  align(dst);
}

auto LabelVector::deSerialize(char*& src) -> LabelVector* {
  LabelVector* lv = new LabelVector();
  memcpy(&(lv->num_bytes_), src, sizeof(lv->num_bytes_));
  src += sizeof(lv->num_bytes_);

  lv->labels_ = new label_t[lv->num_bytes_];
  memcpy(lv->labels_, src, lv->num_bytes_);
  src += lv->num_bytes_;

  align(src);
  return lv;
}

}  // namespace oasis_plus