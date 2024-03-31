#pragma once

#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <mutex>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace benchmark {
inline bool isPointQuery(const uint64_t& a, const uint64_t& b) {
  return b == (a + 1);
}


void intLoadKeys(std::string keyFilePath, std::vector<uint64_t>& keys,
                 std::set<uint64_t>& keyset);

void intLoadQueries(std::string lQueryFilePath, std::string uQueryFilePath,
                    std::vector<std::pair<uint64_t, uint64_t>>& range_queries);


void intLoadKeys(std::string keyFilePath, std::vector<uint64_t>& keys,
                 std::set<uint64_t>& keyset) {
  std::ifstream keyFile;
  uint64_t key;

  keyFile.open(keyFilePath);
  while (keyFile >> key) {
    keyset.insert(key);
    keys.push_back(key);
  }
  keyFile.close();

  std::sort(keys.begin(), keys.end());
}

void intLoadQueries(std::string lQueryFilePath, std::string uQueryFilePath,
                    std::vector<std::pair<uint64_t, uint64_t>>& range_queries) {
  std::ifstream lQueryFile, uQueryFile;
  uint64_t lq, uq;

  lQueryFile.open(lQueryFilePath);
  uQueryFile.open(uQueryFilePath);
  while ((lQueryFile >> lq) && (uQueryFile >> uq)) {
    assert(lq <= uq);
    range_queries.push_back(std::make_pair(lq, uq));
  }
  lQueryFile.close();
  uQueryFile.close();

  std::sort(range_queries.begin(), range_queries.end());
}

}  // namespace benchmark
