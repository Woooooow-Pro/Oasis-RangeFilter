#pragma once

#include <string>

#include "proteus/config.h"
#include "proteus/prefixbf.h"
#include "proteus/rank.h"
#include "proteus/suffix.h"
#include "proteus/surf_builder.h"

namespace oasis_plus {

class LoudsDense {
 public:
  class Iter {
   public:
    Iter()
        : is_valid_(false),
          is_search_complete_(false),
          is_move_left_complete_(false),
          is_move_right_complete_(false),
          prefix_filter_true_(false),
          trie_(nullptr),
          send_out_node_num_(0),
          key_len_(0){};

    Iter(LoudsDense* trie);

    void clear();
    auto isValid() const -> bool { return is_valid_; };
    auto isSearchComplete() const -> bool { return is_search_complete_; };
    auto isMoveLeftComplete() const -> bool { return is_move_left_complete_; };
    auto isMoveRightComplete() const -> bool {
      return is_move_right_complete_;
    };
    auto isComplete() const -> bool {
      return (is_search_complete_ &&
              (is_move_left_complete_ && is_move_right_complete_));
    }

    // PROTEUS
    auto prefixFilterTrue() const -> bool { return prefix_filter_true_; };
    template <typename T>
    auto compare(const T& key, PrefixBF* prefix_filter) const -> int;

    auto getKey() const -> std::string;
    auto getSendOutNodeNum() const -> position_t { return send_out_node_num_; };

    void setToFirstLabelInRoot();
    void setToLastLabelInRoot();
    void moveToLeftMostKey();
    void moveToRightMostKey();
    void operator++(int);
    void operator--(int);

   private:
    inline void append(position_t pos);
    inline void set(level_t level, position_t pos);
    inline void setSendOutNodeNum(position_t node_num) {
      send_out_node_num_ = node_num;
    };
    inline void setFlags(const bool is_valid, const bool is_search_complete,
                         const bool is_move_left_complete,
                         const bool is_move_right_complete,
                         bool prefix_filter_true = false);

   private:
    // True means the iter either points to a valid key
    // or to a prefix with length trie_->getHeight()
    bool is_valid_;
    // If false, call moveToKeyGreaterThan in LoudsSparse to complete
    bool is_search_complete_;
    // If false, call moveToLeftMostKey in LoudsSparse to complete
    bool is_move_left_complete_;
    // If false, call moveToRightMostKey in LoudsSparse to complete
    bool is_move_right_complete_;

    // PROTEUS
    // If true, immediately return true overall
    bool prefix_filter_true_;

    LoudsDense* trie_;
    position_t send_out_node_num_;
    level_t key_len_;  // Does NOT include suffix

    std::vector<label_t> key_;
    std::vector<position_t> pos_in_trie_;

    friend class LoudsDense;
    friend class Proteus;  // Allow LoudsSparse to get access to key prefix in
                           // LoudsDense Iter
  };

 public:
  LoudsDense(){};
  LoudsDense(const SuRFBuilder* builder);

  ~LoudsDense() {}

  // Returns whether key exists in the trie so far
  // out_node_num == 0 means search terminates in louds-dense.
  template <typename T>
  auto lookupKey(const T& key, PrefixBF* prefix_filter,
                 position_t& out_node_num) const -> bool;

  // return value indicates potential false positive
  template <typename T>
  auto moveToKeyGreaterThan(const T& lq, const T& rq, LoudsDense::Iter& iter,
                            PrefixBF* prefix_filter) const -> bool;

  auto getHeight() const -> uint64_t { return height_; };
  auto getTrieDepth() const -> uint32_t { return trie_depth_; };
  auto serializedSize() const -> uint64_t;
  auto getMemoryUsage() const -> uint64_t;

  void serialize(char*& dst) const;

  static auto deSerialize(char*& src, uint32_t trie_depth) -> LoudsDense*;

  void destroy();

 private:
  auto getChildNodeNum(const position_t pos) const -> position_t;
  auto getSuffixPos(const position_t pos) const -> position_t;
  auto getNextPos(const position_t pos) const -> position_t;
  auto getPrevPos(const position_t pos, bool* is_out_of_bound) const
      -> position_t;

  auto compareSuffixGreaterThan(const position_t pos, const level_t level,
                                const uint64_t lq, const uint64_t rq,
                                std::string& edited_lq, LoudsDense::Iter& iter,
                                PrefixBF* prefix_filter) const -> bool;
  auto compareSuffixGreaterThan(const position_t pos, const level_t level,
                                const std::string& lq, const std::string& rq,
                                std::string& edited_lq, LoudsDense::Iter& iter,
                                PrefixBF* prefix_filter) const -> bool;

 private:
  static const position_t kNodeFanout = 256;
  static const position_t kRankBasicBlockSize = 512;

  level_t height_;
  uint32_t trie_depth_;

  BitvectorRank* label_bitmaps_;
  BitvectorRank* child_indicator_bitmaps_;
  BitvectorSuffix* suffixes_;

  /*
          PROTEUS
          prefix_key_indicator_bits is unnecessary as keys
          will be padded out to the same length (the trie depth).
  */
};

}  // namespace oasis_plus