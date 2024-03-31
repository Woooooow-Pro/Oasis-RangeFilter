#pragma once

#include <string>
#include <vector>

#include "oasis_plus.h"
#include "test_wrapper.hpp"
#include "util.hpp"

namespace benchmark {
class OasisPlusWrapper : public TestWrapper {
 public:
  OasisPlusWrapper() = default;
  OasisPlusWrapper(bool is_int_bench) : TestWrapper() {}

  ~OasisPlusWrapper() { delete filter_; }

  void init(const std::vector<std::string> &argv) override {
    bpk_ = strtod(argv[0].c_str(), nullptr);
    block_sz_ = strtoul(argv[1].c_str(), nullptr, 10);

    if (argv.size() > 2) {
      size_t max_qrange = strtoul(argv[2].c_str(), nullptr, 10);
      max_qlen_ = 64ULL - __builtin_clzll(max_qrange) - 1;
    }
  }
  auto construct(std::vector<std::uint64_t> &keys) -> clock_t override;

  /* [left, right) */
  auto query(const uint64_t &left, const uint64_t &right) -> bool override {
    if (isPointQuery(left, right)) {
      return filter_->query(left);
    }
    return filter_->query(left, right - 1);
  }

  auto size() const -> size_t override;

 private:
  double bpk_;
  uint32_t block_sz_;
  size_t max_qlen_ = 10;

  oasis_plus::OasisPlus *filter_ = nullptr;
};

auto OasisPlusWrapper::construct(std::vector<std::uint64_t> &keys) -> clock_t {
  clock_t begin_time;
  begin_time = clock();
  filter_ = new oasis_plus::OasisPlus(bpk_, block_sz_, keys, max_qlen_);
  return clock() - begin_time;
}

auto OasisPlusWrapper::size() const -> size_t { return filter_->size(); }

}  // namespace benchmark