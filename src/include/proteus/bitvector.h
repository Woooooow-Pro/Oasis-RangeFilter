#pragma once

#include <cassert>
#include <cstring>
#include <vector>

#include "proteus/config.h"

namespace oasis_plus {

class Bitvector {
 public:
  Bitvector() : num_bits_(0), bits_(nullptr){};

  Bitvector(const std::vector<std::vector<word_t> >& bitvector_per_level,
            const std::vector<position_t>& num_bits_per_level,
            const level_t start_level, level_t end_level /* non-inclusive */);
  ~Bitvector() {}

  auto numBits() const -> position_t;

  auto numWords() const -> position_t;

  // in bytes
  auto bitsSize() const -> position_t;

  // in bytes
  auto size() const -> position_t;

  auto readBit(const position_t pos) const -> bool;

  auto distanceToNextSetBit(const position_t pos) const -> position_t;
  auto distanceToPrevSetBit(const position_t pos) const -> position_t;

 private:
  auto totalNumBits(const std::vector<position_t>& num_bits_per_level,
                    const level_t start_level,
                    const level_t end_level /* non-inclusive */) -> position_t;

  void concatenateBitvectors(
      const std::vector<std::vector<word_t> >& bitvector_per_level,
      const std::vector<position_t>& num_bits_per_level,
      const level_t start_level, const level_t end_level /* non-inclusive */);

 protected:
  position_t num_bits_;
  word_t* bits_;
};

}  // namespace oasis_plus