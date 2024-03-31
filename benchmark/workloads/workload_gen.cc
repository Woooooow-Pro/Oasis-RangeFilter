#include <assert.h>
#include <sys/stat.h>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>  // C++17 compatible
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <sstream>
#include <unordered_map>
#include <vector>

// using namespace std;

const double readWriteProportion = 0.5;
static size_t SOSD_IDX = 0;
std::vector<uint64_t> SOSD_DATA;

enum kdist_type {
  kuniform,
  knormal,
  ksosd_books,
  ksosd_fb,
};

enum qdist_type { quniform, qsosd_books, qsosd_fb };

static std::vector<std::string> kdist_names = {"kuniform", "knormal",
                                               "ksosd_books", "ksosd_fb"};
static std::vector<std::string> qdist_names = {"quniform", "qsosd_books",
                                               "qsosd_fb"};

static std::unordered_map<std::string, kdist_type> str_to_kdist = {
    {"kuniform", kuniform},
    {"knormal", knormal},
    {"ksosd_books", ksosd_books},
    {"ksosd_fb", ksosd_fb}};
static std::unordered_map<kdist_type, std::string> ksosd_to_file_name = {
    {ksosd_books, "books_800M_uint64"}, {ksosd_fb, "fb_200M_uint64"}};
static std::unordered_map<std::string, qdist_type> str_to_qdist = {
    {"quniform", quniform},
    {"qsosd_books", qsosd_books},
    {"qsosd_fb", qsosd_fb}};

auto generateKeysUniform(uint64_t nkeys, uint64_t kmax)
    -> std::vector<uint64_t> {
  std::vector<uint64_t> keys;
  keys.reserve(nkeys);

  // Will be used to obtain a seed for the random number engine
  std::random_device rd;
  // Standard mersenne_twister_engine seeded with rd()
  std::mt19937_64 gen(rd());
  std::uniform_int_distribution<uint64_t> uni_dist(0, kmax);

  while (keys.size() < nkeys) {
    uint64_t number = uni_dist(gen);
    keys.push_back(number);
  }

  return keys;
}

auto generateKeysNormal(size_t nkeys, uint64_t kmax,
                        long double standard_deviation)
    -> std::vector<uint64_t> {
  std::vector<uint64_t> keys;
  keys.reserve(nkeys);

  // Mean is in the middle of the key space
  std::normal_distribution<long double> nor_dist =
      std::normal_distribution<long double>(1ULL << (64 - 1),
                                            standard_deviation);
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);

  while (keys.size() < nkeys) {
    uint64_t number = (uint64_t)nor_dist(generator);
    if (number <= kmax) {
      keys.push_back(number);
    }
  }

  return keys;
}

auto generateUniqueUnsortedKeys(size_t nkeys, uint64_t kmax)
    -> std::vector<uint64_t> {
  // Will be used to obtain a seed for the random number engine
  std::random_device rd;
  // Standard mersenne_twister_engine seeded with rd()
  std::mt19937_64 gen(rd());
  std::uniform_int_distribution<unsigned long long> uni_dist(0, kmax);
  std::vector<uint64_t> unsort_keys;
  unsort_keys.reserve(nkeys);
  // if the nkeys is very close to key space size, then we need to rethink
  // the following code;
  std::unordered_map<uint64_t, bool> hash;
  while (unsort_keys.size() < nkeys) {
    uint64_t number = uni_dist(gen);
    if (hash.find(number) == hash.end()) {
      unsort_keys.push_back(number);
      hash[number] = true;
    }
  }
  return unsort_keys;
}

auto getSOSDKeys(size_t nkeys) -> std::vector<uint64_t> {
  std::vector<uint64_t> keys(SOSD_DATA.begin() + SOSD_IDX,
                             SOSD_DATA.begin() + SOSD_IDX + nkeys);
  SOSD_IDX += nkeys;
  return keys;
}

auto generateRangeQueries(size_t nqueries, uint64_t min_range,
                          uint64_t max_range, double pqratio,
                          const std::vector<uint64_t>& keys, qdist_type qdist,
                          double pnratio)
    -> std::vector<std::pair<uint64_t, uint64_t>> {
  std::vector<std::pair<uint64_t, uint64_t>> txn_keys;
  std::vector<uint64_t> range_lefts;

  switch (qdist) {
    case qdist_type::quniform:
      range_lefts = generateKeysUniform(
          nqueries, std::numeric_limits<uint64_t>::max() - max_range);
      break;
    case qdist_type::qsosd_books:
    case qdist_type::qsosd_fb:
      range_lefts = getSOSDKeys(nqueries);
    default:
      break;
  }

  uint64_t range_size, left;
  uint64_t index = range_lefts.size() - 1;
  srand(time(NULL));

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<uint64_t> index_dist(0, index);

  for (size_t i = 1; i <= nqueries; i++) {
    // Generate range query size
    double randdbl = index_dist(gen) * 1.0 / index;
    if (randdbl < pqratio) {
      range_size = 1;
    } else if (min_range == max_range) {
      range_size = 2;
    } else {
      range_size =
          (rand() % (max_range - min_range)) + std::max(2UL, min_range);
    }

    // Generate left query bound
    randdbl = index_dist(gen) * 1.0 / index;
    if (randdbl < pnratio) {
      left = range_size > 1 ? keys[index_dist(gen)] - 1 : keys[index_dist(gen)];
    } else {
      left = range_lefts[index_dist(gen)];
    }

    if (std::numeric_limits<uint64_t>::max() - left > range_size) {
      txn_keys.push_back(
          std::pair<uint64_t, uint64_t>(left, left + range_size));
    }
  }

  if (txn_keys.size() != nqueries) {
    std::cout << "WARNING: The number of queries generated is less than the "
                 "specified amount."
              << std::endl;
    std::cout
        << "This is because some queries generated go past the max uint64_t."
        << std::endl;
  }

  return txn_keys;
}

auto generateKeys(size_t nkeys, kdist_type kdist) -> std::vector<uint64_t> {
  switch (kdist) {
    case kdist_type::kuniform:
      return generateKeysUniform(nkeys, std::numeric_limits<uint64_t>::max());
    case kdist_type::knormal:
      return generateKeysNormal(nkeys, std::numeric_limits<uint64_t>::max(),
                                std::numeric_limits<uint64_t>::max() * 0.01L);
    case kdist_type::ksosd_books:
    case kdist_type::ksosd_fb:
      return getSOSDKeys(nkeys);
    default:
      assert(false);
  }
}

// True = Read (Range Query), False = Write (Put)
auto generateTrace(const double read_write_proportion, size_t nkeys,
                   size_t nqueries) -> std::vector<bool> {
  std::default_random_engine generator;
  std::bernoulli_distribution distribution(read_write_proportion);

  std::vector<bool> trace;
  size_t nreads = 0;
  size_t nwrites = 0;

  for (size_t i = 0; i < nkeys + nqueries; i++) {
    if (distribution(generator)) {
      trace.push_back(true);
      nreads++;
    } else {
      trace.push_back(false);
      nwrites++;
    }

    if (nreads == nqueries) {
      for (size_t j = i + 1; j < nkeys + nqueries; j++) {
        trace.push_back(false);
        nwrites++;
      }
      assert(nwrites == nkeys);
      break;
    } else if (nwrites == nkeys) {
      for (size_t j = i + 1; j < nkeys + nqueries; j++) {
        trace.push_back(true);
        nreads++;
      }
      assert(nreads == nqueries);
      break;
    }
  }

  return trace;
}

template <typename T>
void interleaveWorkload(std::vector<T>& v1, std::vector<T>& v2) {
  std::default_random_engine generator;
  std::vector<T> merged;
  size_t nv1 = 0, nv2 = 0, total = v1.size() + v2.size();
  for (size_t i = 0; i < total; i++) {
    std::bernoulli_distribution distribution(((total - i) * 1.0L) / total);
    if (distribution(generator)) {
      merged.push_back(v1[nv1]);
      nv1++;
    } else {
      merged.push_back(v2[nv2]);
      nv2++;
    }

    if (nv1 == v1.size()) {
      for (size_t j = i + 1; j < total; j++) {
        merged.push_back(v2[nv2]);
        nv2++;
      }
      assert(nv2 == v2.size());
      break;
    } else if (nv2 == v2.size()) {
      for (size_t j = i + 1; j < total; j++) {
        merged.push_back(v1[nv1]);
        nv1++;
      }
      assert(nv1 == v1.size());
      break;
    }
  }

  assert(merged.size() == total);
  std::copy(merged.begin(), merged.begin() + v1.size(), v1.begin());
  std::copy(merged.begin() + v1.size(), merged.end(), v2.begin());
}

template <typename T>
void writeValuesToFile(const std::vector<T>& v, const std::string& f) {
  std::stringstream ss;
  ss << "my_data/" << f << ".txt";
  std::ofstream output_file(ss.str());
  std::ostream_iterator<T> output_iterator(output_file, "\n");
  std::copy(v.begin(), v.end(), output_iterator);
  output_file.close();
}

template <typename T>
void writePairsToFile(const std::vector<std::pair<T, T>>& v,
                      const std::string& f1, const std::string& f2) {
  std::stringstream ss1;
  std::stringstream ss2;
  ss1 << "my_data/" << f1 << ".txt";
  ss2 << "my_data/" << f2 << ".txt";
  std::ofstream output_file1(ss1.str());
  std::ofstream output_file2(ss2.str());
  for (const auto& p : v) {
    output_file1 << p.first << "\n";
    output_file2 << p.second << "\n";
  }
  output_file1.close();
  output_file2.close();
}

template <typename T>
void shuffleVector(std::vector<T>& v) {
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::shuffle(std::begin(v), std::end(v), std::default_random_engine(seed));
}

template <typename T>
std::vector<T> parseArg(const char* arg, std::function<T(std::string)> parse) {
  std::string to_parse = std::string(arg);
  std::stringstream ss(to_parse);
  std::vector<T> parsed;
  std::string w;
  while (ss >> w) {
    parsed.push_back(parse(w));
  }
  return parsed;
}

int main(int argc, char* argv[]) {
  assert(argc == 10);

  std::function<uint64_t(std::string)> parse_ull = [](std::string s) {
    return strtoull(s.c_str(), NULL, 0);
  };
  std::function<kdist_type(std::string)> parse_kdist = [](std::string s) {
    assert(str_to_kdist.count(s) != 0);
    return str_to_kdist[s];
  };
  std::function<qdist_type(std::string)> parse_qdist = [](std::string s) {
    assert(str_to_qdist.count(s) != 0);
    return str_to_qdist[s];
  };
  std::function<double(std::string)> parse_double = [](std::string s) {
    return atof(s.c_str());
  };

  const std::string SOSD_DATA_DIR = argv[1];
  std::vector<uint64_t> nkeys = parseArg<uint64_t>(argv[2], parse_ull);
  std::vector<uint64_t> nqueries = parseArg<uint64_t>(argv[3], parse_ull);
  std::vector<uint64_t> min_range = parseArg<uint64_t>(argv[4], parse_ull);
  std::vector<uint64_t> max_range = parseArg<uint64_t>(argv[5], parse_ull);
  std::vector<kdist_type> kdist = parseArg<kdist_type>(argv[6], parse_kdist);
  std::vector<qdist_type> qdist = parseArg<qdist_type>(argv[7], parse_qdist);
  std::vector<double> pqratio = parseArg<double>(argv[8], parse_double);
  std::vector<double> pnratio = parseArg<double>(argv[9], parse_double);

  // Safety Checks
  assert(std::filesystem::exists(SOSD_DATA_DIR));
  assert(nqueries.size() == nkeys.size());
  assert(min_range.size() == nkeys.size());
  assert(max_range.size() == nkeys.size());
  assert(kdist.size() == nkeys.size());
  assert(qdist.size() == nkeys.size());
  assert(pqratio.size() == nkeys.size());
  assert(pnratio.size() == nkeys.size());

  assert(nkeys[0] > 0);
  assert(nqueries[0] > 0);
  for (size_t i = 0; i < nkeys.size(); i++) {
    assert(min_range[i] >= 2);
    assert(min_range[i] <= max_range[i]);
    assert(pqratio[i] >= 0 && pqratio[i] <= 1);
    assert(pnratio[i] >= 0 && pnratio[i] <= 1);

    // Print out to double check
    printf(
        "%zu - SOSD_DIR: %s; KNum: %lu; QNum: %lu; MinR: %lu; MaxR: %lu\n\t"
        "KDist: %s; QDist: %s; PQR: %f; PNR: %f\n",
        i, SOSD_DATA_DIR.c_str(), nkeys[i], nqueries[i], min_range[i],
        max_range[i], kdist_names[kdist[i]].c_str(),
        qdist_names[qdist[i]].c_str(), pqratio[i], pnratio[i]);
  }

  // Shuffle SOSD dataset before using them
  auto it = std::find_if(kdist.begin(), kdist.end(), [](kdist_type kd) {
    return kd == ksosd_books || kd == ksosd_fb;
  });
  if (it != kdist.end()) {
    std::string filename = SOSD_DATA_DIR + ksosd_to_file_name[*it];
    assert(std::filesystem::exists(filename));
    std::fstream f;
    f.open(filename.c_str(), std::ios::out | std::ios::in | std::ios::binary);
    if (!f.is_open()) {
      std::cerr << "Unable to open " << filename << std::endl;
      exit(EXIT_FAILURE);
    }

    // Read size
    uint64_t size;
    f.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
    SOSD_DATA.resize(size);

    // Read values
    f.read(reinterpret_cast<char*>(SOSD_DATA.data()), size * sizeof(uint64_t));

    // Shuffle data
    std::default_random_engine rng(rand());
    std::shuffle(std::begin(SOSD_DATA), std::end(SOSD_DATA), rng);
  }

  // Create data directory
  std::string DIR = "my_data";
  mkdir(DIR.c_str(), 0777);

  std::vector<uint64_t> keys;
  keys.reserve(std::accumulate(nkeys.begin(), nkeys.end(), 0ULL));

  std::vector<uint64_t> prev_keys = generateKeys(nkeys[0], kdist[0]);
  shuffleVector(prev_keys);
  keys.insert(keys.end(), prev_keys.begin(), prev_keys.end());

  std::vector<std::pair<uint64_t, uint64_t>> prev_queries =
      generateRangeQueries(nqueries[0], min_range[0], max_range[0], pqratio[0],
                           keys, qdist[0], pnratio[0]);
  shuffleVector(prev_queries);

  for (size_t i = 1; i < nkeys.size(); i++) {
    std::vector<uint64_t> gen_keys = generateKeys(nkeys[i], kdist[i]);
    shuffleVector(gen_keys);
    keys.insert(keys.end(), gen_keys.begin(), gen_keys.end());

    std::vector<std::pair<uint64_t, uint64_t>> gen_queries =
        generateRangeQueries(nqueries[i], min_range[i], max_range[i],
                             pqratio[i], keys, qdist[i], pnratio[i]);
    shuffleVector(gen_queries);

    // Interleave successive workloads starting from the 2nd one
    if (i > 1) {
      interleaveWorkload(prev_keys, gen_keys);
      interleaveWorkload(prev_queries, gen_queries);
    }

    // Write keys and queries from previous iteration to file
    std::string kfilename = "data" + std::to_string(i - 1);
    std::string qfilename1 = "txn" + std::to_string(i - 1);
    std::string qfilename2 = "upper_bound" + std::to_string(i - 1);
    writeValuesToFile(prev_keys, kfilename);
    writePairsToFile(prev_queries, qfilename1, qfilename2);

    // Generate read-write interleave trace and write it to file
    // True = Read (Range Query), False = Write (Put)
    if (gen_keys.size() > 0 && gen_queries.size() > 0) {
      std::vector<bool> trace = generateTrace(
          readWriteProportion, gen_keys.size(), gen_queries.size());
      std::string filename = "trace" + std::to_string(i);
      writeValuesToFile(trace, filename);
    }

    prev_keys = gen_keys;
    prev_queries = gen_queries;
  }

  std::string kfilename = "data" + std::to_string(nkeys.size() - 1);
  std::string qfilename1 = "txn" + std::to_string(nkeys.size() - 1);
  std::string qfilename2 = "upper_bound" + std::to_string(nkeys.size() - 1);
  writeValuesToFile(prev_keys, kfilename);
  writePairsToFile(prev_queries, qfilename1, qfilename2);

  return 0;
}
