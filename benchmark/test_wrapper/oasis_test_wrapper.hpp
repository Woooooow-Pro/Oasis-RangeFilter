#pragma once

#include <string>
#include <vector>

#include "oasis/oasis.h"
#include "test_wrapper.hpp"
#include "util.hpp"

namespace benchmark {
class OasisWrapper : public TestWrapper {
 public:
  OasisWrapper() = default;
  OasisWrapper(bool is_int_bench) : TestWrapper() {}

  ~OasisWrapper() { delete filter_; }

  void init(const std::vector<std::string> &argv) override {
    bpk_ = strtod(argv[0].c_str(), nullptr);
    block_sz_ = strtoul(argv[1].c_str(), nullptr, 10);
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

  oasis::Oasis<uint64_t> *filter_ = nullptr;
};

auto OasisWrapper::construct(std::vector<std::uint64_t> &keys) -> clock_t {
  clock_t begin_time;

  begin_time = clock();
  filter_ = new oasis::Oasis<uint64_t>(bpk_, block_sz_, keys);
  return clock() - begin_time;
}

auto OasisWrapper::size() const -> size_t { return filter_->size(); }
}  // namespace benchmark