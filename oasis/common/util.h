#pragma once

#include <cstdint>
#include <vector>

namespace oasis {

inline auto align_bit2byte(uint64_t b_size) -> size_t {
  return (b_size + 7ULL) >> 3;
}

inline void align(uint8_t *&ptr) {
  ptr = reinterpret_cast<uint8_t *>((reinterpret_cast<uint64_t>(ptr) + 7ULL) &
                                    ~(7ULL));
}

inline void sizeAlign(uint32_t &size) { size = (size + 7UL) & ~(7UL); }
inline void sizeAlign(uint64_t &size) { size = (size + 7ULL) & ~(7ULL); }

}  // namespace oasis
