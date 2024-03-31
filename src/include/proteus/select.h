#pragma once

#include <vector>

#include "proteus/bitvector.h"
#include "proteus/config.h"

namespace oasis_plus {

class BitvectorSelect : public Bitvector {
 public:
  BitvectorSelect() : sample_interval_(0), num_ones_(0), select_lut_(nullptr){};

  BitvectorSelect(const position_t sample_interval,
                  const std::vector<std::vector<word_t> >& bitvector_per_level,
                  const std::vector<position_t>& num_bits_per_level,
                  const level_t start_level,
                  const level_t end_level /* non-inclusive */);

  ~BitvectorSelect() {}

  // Returns the postion of the rank-th 1 bit.
  // posistion is zero-based; rank is one-based.
  // E.g., for bitvector: 100101000, select(3) = 5
  auto select(position_t rank) const -> position_t;

  auto selectLutSize() const -> position_t;

  auto serializedSize() const -> position_t;

  auto size() const -> position_t;

  auto numOnes() const -> position_t { return num_ones_; }

  void serialize(char*& dst) const;

  static auto deSerialize(char*& src) -> BitvectorSelect*;

  void destroy() {
    delete[] bits_;
    delete[] select_lut_;
  }

 private:
  // This function currently assumes that the first bit in the
  // bitvector is one.
  void initSelectLut();

 private:
  position_t sample_interval_;
  position_t num_ones_;
  position_t* select_lut_;  // select look-up table
};

}  // namespace oasis_plus
