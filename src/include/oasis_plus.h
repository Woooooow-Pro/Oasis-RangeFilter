#pragma once

#include <cstdint>
#include <vector>

#include "learned_rf/learned_rf.h"
#include "proteus/proteus.h"

namespace oasis_plus {
class OasisPlus {
 private:
  static const uint8_t LEARNED_EXIST_BIT = 1U;
  static const uint8_t PROTEUS_EXIST_BIT = 2U;

 public:
  OasisPlus(double bpk, uint32_t block_size, const std::vector<uint64_t> &keys,
            const size_t max_qlen = 10);

  OasisPlus(std::vector<uint64_t> &begins, std::vector<uint64_t> &ends,
            std::vector<uint32_t> &filter_types, LearnedRF *learned_rf,
            Proteus *&proteus)
      : begins_(std::move(begins)),
        ends_(std::move(ends)),
        filter_types_(std::move(filter_types)),
        learned_rf_(learned_rf),
        proteus_(proteus) {}

  ~OasisPlus() {
    delete learned_rf_;
    delete proteus_;
  }

  auto query(uint64_t key) -> bool;
  auto query(uint64_t l_key, uint64_t r_key) -> bool;

  auto serialize() const -> std::pair<uint8_t *, size_t>;
  static auto deserialize(uint8_t *ser) -> OasisPlus *;

  auto size() const -> size_t;

 private:
  /* Indices of the Intervals */
  std::vector<uint64_t> begins_;
  std::vector<uint64_t> ends_;
  /* 0 for Proteus, others indicate leanred filter's indice */
  std::vector<uint32_t> filter_types_;

  LearnedRF *learned_rf_;
  Proteus *proteus_;
};

}  // namespace oasis_plus