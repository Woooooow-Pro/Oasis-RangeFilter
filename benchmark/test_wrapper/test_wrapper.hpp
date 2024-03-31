#pragma once

#include <cstdint>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace benchmark {
class TestWrapper {
 public:
  /** Set the int bench as default option */
  TestWrapper() = default;

  virtual ~TestWrapper() = default;

  virtual void init(const std::vector<std::string> &argv) {}
  virtual void init(const std::vector<std::pair<uint64_t, uint64_t>> &queries,
                    const std::vector<std::string> &argv) {}

  virtual auto construct(std::vector<std::uint64_t> &keys) -> clock_t = 0;

  /* [left, right) */
  virtual auto query(const uint64_t &left, const uint64_t &right) -> bool = 0;

  auto query_time(const std::vector<std::pair<uint64_t, uint64_t>> &queries)
      -> clock_t;

  auto fpr_counter(const std::vector<std::pair<uint64_t, uint64_t>> &queries,
                   std::set<uint64_t> &keyset)
      -> std::tuple<size_t, size_t, size_t>;

  virtual auto size() const -> size_t = 0;

 public:
};

auto TestWrapper::query_time(
    const std::vector<std::pair<uint64_t, uint64_t>> &queries) -> clock_t {
  clock_t begin_time = clock();
  for (const std::pair<uint64_t, uint64_t> &q : queries) {
    query(q.first, q.second);
  }

  return clock() - begin_time;
}

auto TestWrapper::fpr_counter(
    const std::vector<std::pair<uint64_t, uint64_t>> &queries,
    std::set<uint64_t> &keyset) -> std::tuple<size_t, size_t, size_t> {
  size_t neg = 0;
  size_t fp = 0;
  size_t fn = 0;
  for (const std::pair<uint64_t, uint64_t> &q : queries) {
    bool filter_ans = query(q.first, q.second);

    auto it = keyset.lower_bound(q.first);
    bool full = it != keyset.end() && ((*it) < q.second);
    if (!full) {
      ++neg;
      fp += filter_ans ? 1 : 0;
    } else if (!filter_ans) {
      ++fn;
      printf("False negative!\t[%lu, %lu)\n", q.first, q.second);
    }
  }
  return std::make_tuple(neg, fp, fn);
}

}  // namespace benchmark