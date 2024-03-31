#pragma once

#include <string>
#include <vector>

#include "proteus/config.h"

namespace oasis_plus {

const uint32_t MAX_PBF_HASH_FUNCS = 32;

class PrefixBF {
 public:
  /*
      Hash key prefixes of specified length into the Bloom filter
      Based on the number of unique key prefixes hashed, the optimal
      number of hash functions is calculated and the appropriate number
      of seeds is generated. We set a maximum of 32 hash functions to bound
      filter latency when the number of filter elements is small. This can
      happen if a shorter prefix length is chosen.
  */
  PrefixBF(uint32_t prefix_len, uint64_t nbits,
           const std::vector<uint64_t>& keys);

  PrefixBF(uint32_t prefix_len, uint64_t nbits,
           const std::vector<std::string>& keys);

  PrefixBF(uint32_t prefix_len, uint8_t* data, std::vector<uint32_t> seeds32,
           std::vector<std::pair<uint64_t, uint64_t>> seeds64, uint64_t nmod);

  ~PrefixBF() { delete[] data_; }

  uint32_t getPrefixLen() const { return prefix_len_; }

  auto Query(const uint64_t key, bool shift = true) -> bool;
  auto Query(const uint64_t from, const uint64_t to) -> bool;

  auto Query(const std::string& key) -> bool;
  auto Query(const std::string& from, const std::string& to) -> bool;

  auto serialize() const -> std::pair<uint8_t*, uint64_t>;
  static auto deserialize(uint8_t* ser) -> std::pair<PrefixBF*, uint64_t>;

  auto get(uint64_t i) const -> bool;
  void set(uint64_t i, bool v);

  auto hash(const uint64_t edited_key, const uint32_t& seed) -> uint64_t;

 private:
  uint32_t prefix_len_;
  uint8_t* data_;
  // Number of hash functions is given by the size of the seed vector
  // (seeds32_ for uint64 keys, and seeds128_ for string keys)
  std::vector<uint32_t> seeds32_;
  std::vector<std::pair<uint64_t, uint64_t>> seeds64_;
  std::vector<void*> seeds128_;
  uint64_t nmod_;
};

}  // namespace oasis_plus
