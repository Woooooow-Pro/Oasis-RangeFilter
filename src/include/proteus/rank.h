#pragma once

#include <cstring>
#include <vector>

#include "proteus/bitvector.h"
#include "proteus/config.h"

namespace oasis_plus {

class BitvectorRank : public Bitvector {
 public:
  BitvectorRank() : basic_block_size_(0), rank_lut_(nullptr){};

  BitvectorRank(const position_t basic_block_size,
                const std::vector<std::vector<word_t> >& bitvector_per_level,
                const std::vector<position_t>& num_bits_per_level,
                const level_t start_level,
                const level_t end_level /* non-inclusive */);

  ~BitvectorRank() = default;

  // Counts the number of 1's in the bitvector up to position pos.
  // pos is zero-based; count is one-based.
  // E.g., for bitvector: 100101000, rank(3) = 2
  auto rank(position_t pos) const -> position_t;

  auto rankLutSize() const -> position_t;

  auto serializedSize() const -> position_t;

  auto size() const -> position_t;

  void prefetch(position_t pos) const;

  void serialize(char*& dst) const;

  static auto deSerialize(char*& src) -> BitvectorRank*;

  void destroy() {
    delete[] bits_;
    delete[] rank_lut_;
  }

 private:
  void initRankLut();

  position_t basic_block_size_;
  position_t* rank_lut_;  // rank look-up table
};

}  // namespace oasis_plus
