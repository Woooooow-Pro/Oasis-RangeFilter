#pragma once

#include <vector>

#include "proteus/bitvector.h"
#include "proteus/config.h"

namespace oasis_plus {

// Max suffix_len_ = 64 bits
// For kReal suffixes, if the stored key is not long enough to provide
// suffix_len_ suffix bits, its suffix field is cleared (i.e., all 0's)
// to indicate that there is no suffix info associated with the key.
class BitvectorSuffix : public Bitvector {
 private:
  level_t start_level_;
  std::vector<position_t> num_suffixes_per_level_;

 public:
  BitvectorSuffix() {}

  BitvectorSuffix(const std::vector<std::vector<word_t> >& bitvector_per_level,
                  const std::vector<position_t>& num_bits_per_level,
                  const std::vector<position_t> num_suffixes_per_level,
                  const level_t start_level,
                  level_t end_level /* non-inclusive */)
      : Bitvector(bitvector_per_level, num_bits_per_level, start_level,
                  end_level),
        start_level_(start_level),
        num_suffixes_per_level_(num_suffixes_per_level) {}

  static auto constructSuffix(const std::string& key, const level_t level,
                              const level_t len) -> word_t;

  // Determines the start position of the desired suffix within the suffix
  // bitvector
  auto calcBitPos(const position_t idx, const level_t level,
                  const uint32_t trie_depth) const -> position_t;

  auto getSuffixLen(const level_t level, const uint32_t trie_depth) const
      -> level_t;

  auto serializedSize() const -> position_t;

  auto size() const -> position_t;

  auto read(const position_t bit_pos, const level_t level,
            const uint32_t trie_depth) const -> word_t;
  auto checkEquality(const position_t idx, const std::string& key,
                     const level_t level, const uint32_t trie_depth) const
      -> bool;

  // Compare stored suffix to querying suffix.
  auto compare(const position_t idx, const std::string& key,
               const level_t level, const uint32_t trie_depth) const -> int;

  void serialize(char*& dst) const;

  static auto deSerialize(char*& src) -> BitvectorSuffix*;

  void destroy();
};

}  // namespace oasis_plus