#pragma once

#include <emmintrin.h>

#include <vector>

#include "proteus/config.h"

namespace oasis_plus {

class LabelVector {
 public:
  LabelVector() : num_bytes_(0), labels_(nullptr){};

  LabelVector(const std::vector<std::vector<label_t> >& labels_per_level,
              const level_t start_level, level_t end_level /* non-inclusive */);

  ~LabelVector() {}

  auto getNumBytes() const -> position_t { return num_bytes_; }

  auto serializedSize() const -> position_t;

  auto size() const -> position_t;

  auto read(const position_t pos) const -> label_t;

  auto operator[](const position_t pos) const -> label_t;

  auto search(const label_t target, position_t& pos,
              const position_t search_len) const -> bool;
  auto searchGreaterThan(const label_t target, position_t& pos,
                         const position_t search_len) const -> bool;

  auto binarySearch(const label_t target, position_t& pos,
                    const position_t search_len) const -> bool;
  auto simdSearch(const label_t target, position_t& pos,
                  const position_t search_len) const -> bool;
  auto linearSearch(const label_t target, position_t& pos,
                    const position_t search_len) const -> bool;

  auto binarySearchGreaterThan(const label_t target, position_t& pos,
                               const position_t search_len) const -> bool;
  auto linearSearchGreaterThan(const label_t target, position_t& pos,
                               const position_t search_len) const -> bool;

  void serialize(char*& dst) const;

  static auto deSerialize(char*& src) -> LabelVector*;

  void destroy() { delete[] labels_; }

 private:
  position_t num_bytes_;
  label_t* labels_;
};

}  // namespace oasis_plus