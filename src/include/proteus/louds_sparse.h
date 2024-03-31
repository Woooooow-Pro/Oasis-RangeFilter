#pragma once

#include <climits>
#include <string>

#include "proteus/config.h"
#include "proteus/label_vector.h"
#include "proteus/prefixbf.h"
#include "proteus/rank.h"
#include "proteus/select.h"
#include "proteus/suffix.h"
#include "proteus/surf_builder.h"

namespace oasis_plus {

class LoudsSparse {
 public:
  class Iter {
   public:
    Iter()
        : is_valid_(false),
          is_done_(false),
          trie_(nullptr),
          start_level_(0),
          start_node_num_(0),
          key_len_(0){};

    Iter(LoudsSparse* trie);

    void clear();
    auto isValid() const -> bool { return is_valid_; };
    auto isDone() const -> bool { return is_done_; };

    // PROTEUS
    template <typename T>
    auto compare(const T& key, PrefixBF* prefix_filter,
                 const std::string& dense_prefix) const -> int;

    auto getKey() const -> std::string;

    auto getStartNodeNum() const -> position_t { return start_node_num_; };
    void setStartNodeNum(position_t node_num) { start_node_num_ = node_num; };

    void setToFirstLabelInRoot();
    void setToLastLabelInRoot();
    void moveToLeftMostKey();
    void moveToRightMostKey();
    void operator++(int);
    void operator--(int);

   private:
    void append(const position_t pos);
    void append(const label_t label, const position_t pos);
    void set(const level_t level, const position_t pos);

   private:
    bool is_valid_;  // True means the iter currently points to a valid key

    // PROTEUS
    bool is_done_;  // True means range query is done and is true overall

    LoudsSparse* trie_;
    level_t start_level_;
    position_t start_node_num_;  // Passed in by the dense iterator; default = 0
    level_t
        key_len_;  // Start counting from start_level_; does NOT include suffix

    std::vector<label_t> key_;
    std::vector<position_t> pos_in_trie_;

    friend class LoudsSparse;
  };

 public:
  LoudsSparse(){};
  LoudsSparse(const SuRFBuilder* builder);

  ~LoudsSparse() {}

  // point query: trie walk starts at node "in_node_num" instead of root
  // in_node_num is provided by louds-dense's lookupKey function
  template <typename T>
  auto lookupKey(const T& key, PrefixBF* prefix_filter,
                 const position_t in_node_num) const -> bool;

  // return value indicates potential false positive
  template <typename T>
  auto moveToKeyGreaterThan(const T& lq, const T& rq, LoudsSparse::Iter& iter,
                            PrefixBF* prefix_filter) const -> bool;

  auto getHeight() const -> level_t { return height_; };
  auto getStartLevel() const -> level_t { return start_level_; };
  auto getTrieDepth() const -> uint32_t { return trie_depth_; };
  auto serializedSize() const -> uint64_t;
  auto getMemoryUsage() const -> uint64_t;

  void serialize(char*& dst) const;
  static auto deSerialize(char*& src, uint32_t trie_depth) -> LoudsSparse*;

  void destroy();

 private:
  auto getChildNodeNum(const position_t pos) const -> position_t;
  auto getFirstLabelPos(const position_t node_num) const -> position_t;
  auto getLastLabelPos(const position_t node_num) const -> position_t;
  auto getSuffixPos(const position_t pos) const -> position_t;
  auto nodeSize(const position_t pos) const -> position_t;

  void moveToLeftInNextSubtrie(position_t pos, const position_t node_size,
                               const label_t label,
                               LoudsSparse::Iter& iter) const;

  // return value indicates potential false positive
  auto compareSuffixGreaterThan(const position_t pos, const level_t level,
                                const uint64_t lq, const uint64_t rq,
                                std::string& edited_lq, LoudsSparse::Iter& iter,
                                PrefixBF* prefix_filter) const -> bool;
  auto compareSuffixGreaterThan(const position_t pos, const level_t level,
                                const std::string& lq, const std::string& rq,
                                std::string& edited_lq, LoudsSparse::Iter& iter,
                                PrefixBF* prefix_filter) const -> bool;

 private:
  static const position_t kRankBasicBlockSize = 512;
  static const position_t kSelectSampleInterval = 64;

  level_t height_;       // trie height
  level_t start_level_;  // louds-sparse encoding starts at this level
  // number of nodes in louds-dense encoding
  position_t node_count_dense_;
  // number of children(1's in child indicator bitmap) in louds-dense encoding
  position_t child_count_dense_;

  uint32_t trie_depth_;

  LabelVector* labels_;
  BitvectorRank* child_indicator_bits_;
  BitvectorSelect* louds_bits_;
  BitvectorSuffix* suffixes_;
};

}  // namespace oasis_plus