#pragma once

#include <string>
#include <vector>

#include "proteus/config.h"

namespace oasis_plus {

class SuRFBuilder {
 public:
  /*
      Proteus only uses Real suffixes to create a uniform depth trie.
      Instead of using a sparse-dense ratio like SuRF, Proteus' modeling
      is able to determine the optimal sparse-dense cutoff for the given
      dataset and specified trie depth.
  */
  explicit SuRFBuilder(size_t sparse_dense_cutoff, size_t trie_depth)
      : sparse_dense_cutoff_(sparse_dense_cutoff), trie_depth_(trie_depth){};

  ~SuRFBuilder(){};

  // Fills in the LOUDS-dense and sparse vectors (members of this class)
  // through a single scan of the sorted key list.
  // After build, the member vectors are used in SuRF constructor.
  // REQUIRED: provided key list must be sorted.
  template <typename T>
  void build(const std::vector<T>& keys);

  static auto readBit(const std::vector<word_t>& bits, const position_t pos)
      -> bool;

  static void setBit(std::vector<word_t>& bits, const position_t pos);

  level_t getTreeHeight() const { return labels_.size(); }

  // const accessors
  auto getBitmapLabels() const -> const std::vector<std::vector<word_t> >& {
    return bitmap_labels_;
  }
  auto getBitmapChildIndicatorBits() const
      -> const std::vector<std::vector<word_t> >& {
    return bitmap_child_indicator_bits_;
  }
  auto getLabels() const -> const std::vector<std::vector<label_t> >& {
    return labels_;
  }
  auto getChildIndicatorBits() const
      -> const std::vector<std::vector<word_t> >& {
    return child_indicator_bits_;
  }
  auto getLoudsBits() const -> const std::vector<std::vector<word_t> >& {
    return louds_bits_;
  }
  auto getSuffixes() const -> const std::vector<std::vector<word_t> >& {
    return suffixes_;
  }
  auto getSuffixCounts() const -> const std::vector<position_t>& {
    return suffix_counts_;
  }
  auto getNodeCounts() const -> const std::vector<position_t>& {
    return node_counts_;
  }

  // PROTEUS
  auto getSuffixLen(level_t level) const -> level_t {
    // Returns the suffix length for a node terminating at the supplied level
    // (which is 0-indexed). A node that terminates at level 1 has 1 key byte,
    // etc.
    return (level * 8 < trie_depth_) ? trie_depth_ - (level * 8) : 0;
  }
  auto getTrieDepth() const -> level_t { return trie_depth_; }
  auto getSparseDenseCutoff() const -> level_t { return sparse_dense_cutoff_; }
  auto isSameEditedKey(const uint64_t a, const uint64_t b) -> bool {
    return (a >> (64 - trie_depth_)) == (b >> (64 - trie_depth_));
  }
  auto isSameEditedKey(const std::string& a, const std::string& b) -> bool {
    return compare(a, b, trie_depth_) == 0;
  }

 private:
  static auto isSameKey(const std::string& a, const std::string& b) -> bool {
    return a.compare(b) == 0;
  }

  // Fill in the LOUDS-Sparse vectors through a single scan
  // of the sorted key list.
  template <typename T>
  void buildSparse(const std::vector<T>& keys);

  // Walks down the current partially-filled trie by comparing key to
  // its previous key in the list until their prefixes do not match.
  // The previous key is stored as the last items in the per-level
  // label vector.
  // For each matching prefix byte(label), it sets the corresponding
  // child indicator bit to 1 for that label.
  auto skipCommonPrefix(const std::string& key) -> level_t;

  // Starting at the start_level of the trie, the function inserts
  // key bytes to the trie vectors until the first byte/label where
  // key and next_key do not match.
  // This function is called after skipCommonPrefix. Therefore, it
  // guarantees that the stored prefix of key is unique in the trie.
  auto insertKeyBytesToTrieUntilUnique(const std::string& key,
                                       const std::string& next_key,
                                       const level_t start_level) -> level_t;

  // PROTEUS - Fills in suffix bytes for key up to trie_depth
  inline void insertSuffix(const std::string& key, const level_t level);

  inline auto isCharCommonPrefix(const label_t c, const level_t level) const
      -> bool;
  inline auto isLevelEmpty(const level_t level) const -> bool;
  inline void moveToNextItemSlot(const level_t level);
  void insertKeyByte(const char c, const level_t level,
                     const bool is_start_of_node);
  inline void storeSuffix(const level_t level, const word_t suffix);

  // Fill in the LOUDS-Dense vectors based on the built sparse vectors.
  void buildDense();

  void initDenseVectors(const level_t level);
  void setLabelAndChildIndicatorBitmap(const level_t level,
                                       const position_t node_num,
                                       const position_t pos);

  auto getNumItems(const level_t level) const -> position_t;
  void addLevel();
  auto isStartOfNode(const level_t level, const position_t pos) const -> bool;

 private:
  // trie level < sparse_dense_cutoff_: LOUDS-Dense
  // trie level >= sparse_dense_cutoff_: LOUDS-Sparse
  size_t sparse_dense_cutoff_;
  size_t trie_depth_;

  // LOUDS-Sparse bit/byte vectors
  std::vector<std::vector<label_t> > labels_;
  std::vector<std::vector<word_t> > child_indicator_bits_;
  std::vector<std::vector<word_t> > louds_bits_;

  // LOUDS-Dense bit vectors
  std::vector<std::vector<word_t> > bitmap_labels_;
  std::vector<std::vector<word_t> > bitmap_child_indicator_bits_;

  std::vector<std::vector<word_t> > suffixes_;
  std::vector<position_t> suffix_counts_;

  // auxiliary per level bookkeeping vectors
  std::vector<position_t> node_counts_;
};

}  // namespace oasis_plus