#pragma once

#include <string>
#include <vector>

#include "proteus/config.h"
#include "proteus/louds_dense.h"
#include "proteus/louds_sparse.h"
#include "proteus/prefixbf.h"
#include "proteus/surf_builder.h"

namespace oasis_plus {

class Proteus {
 public:
  class Iter {
   public:
    Iter(){};
    Iter(const Proteus* filter) {
      dense_iter_ = filter->validLoudsDense()
                        ? LoudsDense::Iter(filter->louds_dense_)
                        : LoudsDense::Iter();
      sparse_iter_ = filter->validLoudsSparse()
                         ? LoudsSparse::Iter(filter->louds_sparse_)
                         : LoudsSparse::Iter();
      could_be_fp_ = false;
    }

    auto prefixFilterTrue(bool valid_dense, bool valid_sparse) const -> bool;
    void clear(bool valid_dense, bool valid_sparse);
    auto isValid(bool valid_dense, bool valid_sparse) const -> bool;
    auto getFpFlag() const -> bool;

    // PROTEUS
    template <typename T>
    auto compare(const T& key, bool valid_dense, bool valid_sparse,
                 PrefixBF* prefix_filter) const -> int;

    auto getKey() const -> std::string;

   private:
    void passToSparse();
    auto incrementDenseIter() -> bool;

   private:
    // true implies that dense_iter_ is valid
    LoudsDense::Iter dense_iter_;
    LoudsSparse::Iter sparse_iter_;
    bool could_be_fp_;

    friend class Proteus;
  };

 public:
  Proteus(){};

  //------------------------------------------------------------------
  // Input keys must be SORTED
  //------------------------------------------------------------------

  template <typename T>
  Proteus(const std::vector<T>& keys, const size_t trie_depth,
          const size_t sparse_dense_cutoff, const size_t prefix_length,
          const double bpk);

  ~Proteus();

  // PROTEUS
  template <typename T>
  auto Query(const T& key) const -> bool;
  template <typename T>
  auto Query(const T& left_key, const T& right_key) -> bool;

  auto trieSerializedSize() const -> uint64_t;
  auto getMemoryUsage() const -> uint64_t;
  auto getHeight() const -> level_t;
  auto getSparseStartLevel() const -> level_t;

  auto serialize() const -> std::pair<uint8_t*, size_t>;

  static auto deSerialize(char* src) -> std::pair<Proteus*, size_t>;

  void destroy();

  auto validLoudsDense() const -> bool { return sparse_dense_cutoff_ > 0; }

  auto validLoudsSparse() const -> bool {
    return sparse_dense_cutoff_ < ((trie_depth_ + 7) / 8);
  }

 private:
  LoudsDense* louds_dense_;
  LoudsSparse* louds_sparse_;
  SuRFBuilder* builder_;
  Proteus::Iter iter_;
  PrefixBF* prefix_filter;
  uint32_t trie_depth_;
  uint32_t sparse_dense_cutoff_;
};

}  // namespace oasis_plus
