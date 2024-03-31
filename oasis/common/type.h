#pragma once

#include <cstdint>
#include <vector>

namespace oasis {
using bitmap_t = uint64_t;

using BitPos = bitmap_t;
using BitsPos = std::vector<BitPos>;
}  // namespace oasis