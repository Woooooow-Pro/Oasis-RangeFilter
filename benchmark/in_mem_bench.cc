#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

#include "test_wrapper/oasis_plus_test_wrapper.hpp"
#include "test_wrapper/oasis_test_wrapper.hpp"
#include "test_wrapper/util.hpp"

namespace benchmark {
enum FilterType { Oasis, OasisPlus };

static std::unordered_map<std::string, FilterType> kName2Type{
    {"Oasis", FilterType::Oasis}, {"OasisPlus", FilterType::OasisPlus}};

const std::string dataPath = "./my_data/";

const std::string keyFilePath = dataPath + "data0.txt";
const std::string lQueryFilePath = dataPath + "txn0.txt";
const std::string uQueryFilePath = dataPath + "upper_bound0.txt";

class Experiment {
 public:
  Experiment(std::vector<uint64_t> &keys, std::set<uint64_t> &keyset,
             std::vector<std::pair<uint64_t, uint64_t>> &queries,
             std::vector<std::string> &argv)
      : keys_(std::move(keys)),
        keyset_(std::move(keyset)),
        queries_(std::move(queries)) {
    csv_file_ = argv[0];

    filter_type_ = kName2Type[argv[1]];
    argv.erase(argv.begin(), argv.begin() + 2);
    switch (filter_type_) {
      case FilterType::OasisPlus:
        assert(argv.size() >= 3);
        test_wrapper_ = new OasisPlusWrapper();
        test_wrapper_->init(argv);
        break;
      case FilterType::Oasis:
        assert(argv.size() >= 2);
        test_wrapper_ = new OasisWrapper();
        test_wrapper_->init(argv);
        break;
    }
    rescsv.open(csv_file_.c_str(), std::ios::app);
  }

  ~Experiment() {
    delete test_wrapper_;
    rescsv.close();
  }

  void test() {
    construct_time_ = test_wrapper_->construct(keys_);
    query_time_ = test_wrapper_->query_time(queries_);
    auto fpr_result = test_wrapper_->fpr_counter(queries_, keyset_);
    size_ = test_wrapper_->size();
    neg_num_ = std::get<0>(fpr_result);
    fp_num_ = std::get<1>(fpr_result);
    fn_num_ = std::get<2>(fpr_result);
  }

  void print_result_human() {
    printf("s/Constructing Time:\t%lf\n",
           (double)construct_time_ / CLOCKS_PER_SEC);
    printf("us/Query:\t%lf\n",
           (double)1000000 * query_time_ / CLOCKS_PER_SEC / queries_.size());
    printf("FPR:\t%lf\n", (double)fp_num_ / neg_num_);
    printf("\tNeg: %zu, fn: %lu, fp: %lu\n", neg_num_, fn_num_, fp_num_);
    printf("BPK:\t%lf\n", (double)size_ * 8 / keys_.size());
  }

  void print_csv() {
    rescsv << (double)size_ * 8 / keys_.size() /* Acctual bpk of the filter */
           << ","
           << (double)construct_time_ /
                  CLOCKS_PER_SEC /* Filter construction time */
           << ","
           << (double)1000000 * query_time_ / CLOCKS_PER_SEC /
                  queries_.size()               /* avg query time */
           << "," << (double)fp_num_ / neg_num_ /* FPR */
           << "," << std::endl;
  }

 private:
  /* exp params */
  FilterType filter_type_;
  std::vector<uint64_t> keys_;
  std::set<uint64_t> keyset_;
  std::vector<std::pair<uint64_t, uint64_t>> queries_;

  std::string csv_file_;

  /* exp data */
  clock_t construct_time_;
  clock_t query_time_;
  size_t neg_num_;
  size_t fp_num_;
  size_t fn_num_;
  size_t size_;

  TestWrapper *test_wrapper_ = nullptr;
  std::ofstream rescsv;
};

}  // namespace benchmark

int main(int argc, char *argv[]) {
  using namespace benchmark;
  std::vector<std::string> v_argv(argc - 1);
  for (int i = 1; i < argc; ++i) {
    v_argv[i - 1] = argv[i];
  }

  std::vector<uint64_t> keys;
  std::set<uint64_t> keyset;
  std::vector<std::pair<uint64_t, uint64_t>> queries;

  intLoadKeys(keyFilePath, keys, keyset);
  intLoadQueries(lQueryFilePath, uQueryFilePath, queries);
  Experiment exp(keys, keyset, queries, v_argv);
  exp.test();
  exp.print_csv();
  exp.print_result_human();
  return 0;
}
