#pragma once

#include <memory>
#include <vector>

#include "bitset.h"
#include "common/type.h"
#include "common/util.h"

namespace oasis {

class BitMap {
 public:
  BitMap() = default;
  virtual ~BitMap() = default;

  virtual void build(const BitsPos &bits_pos) = 0;

  virtual void destroy() {}

  virtual auto query(BitPos key) -> bool = 0;
  virtual auto query(BitPos left_key, BitPos right_key) -> bool = 0;

  virtual auto size() const -> size_t = 0;

  virtual auto serialize() const -> std::pair<uint8_t *, size_t> = 0;

  /*
   * For any type of bitmap you must implement a deserializer!
   *  static auto deserialize(uint8_t *ser) -> std::unique_ptr<BitMap>;
   */
};

class BlockEliasFanoBitMap : public BitMap {
  using Blocks = std::vector<BitSet>;

 public:
  BlockEliasFanoBitMap(uint16_t block_sz) : block_sz_(block_sz) {}

  BlockEliasFanoBitMap(uint16_t block_sz, size_t bitmap_sz,
                       uint16_t last_block_sz, uint8_t *bitmap_ptr,
                       Blocks &blocks, BitsPos &blocks_bias)
      : block_sz_(block_sz),
        bitmap_sz_(bitmap_sz),
        last_block_sz_(last_block_sz),
        bitmap_ptr_(bitmap_ptr),
        blocks_(std::move(blocks)),
        blocks_bias_(std::move(blocks_bias)) {}

  ~BlockEliasFanoBitMap() { delete[] bitmap_ptr_; }

  void build(const BitsPos &bits_pos) override {
    last_block_sz_ = bits_pos.size() % block_sz_;
    std::vector<uint8_t> compressed_bitmap;

    BitPos low_bound = bits_pos[0];
    blocks_bias_.emplace_back(low_bound);
    for (size_t idx = 1; idx < bits_pos.size();) {
      BitsPos cur_batch;
      cur_batch.emplace_back(0);
      for (size_t cnt = 1; cnt < block_sz_ && idx < bits_pos.size(); ++cnt) {
        cur_batch.emplace_back(bits_pos[idx++] - low_bound);
      }
      if (idx < bits_pos.size()) {
        blocks_bias_.emplace_back(bits_pos[idx++]);
      } else {
        blocks_bias_.emplace_back(bits_pos.back());
      }
      std::vector<uint8_t> batch_block =
          BitSet::build(cur_batch, blocks_bias_.back() - low_bound);
      compressed_bitmap.insert(compressed_bitmap.end(), batch_block.begin(),
                               batch_block.end());
      low_bound = blocks_bias_.back();
    }

    bitmap_sz_ = compressed_bitmap.size();
    bitmap_ptr_ = new uint8_t[bitmap_sz_];
    memcpy(bitmap_ptr_, compressed_bitmap.data(), bitmap_sz_ * sizeof(uint8_t));

    size_t nbatches = blocks_bias_.size() - 1;
    uint8_t *pos = bitmap_ptr_;
    for (size_t i = 0; i < nbatches; ++i) {
      blocks_.emplace_back(i == nbatches - 1 ? last_block_sz_ : block_sz_,
                           blocks_bias_[i + 1] - blocks_bias_[i], pos);
      pos += blocks_.back().size();
    }
  }

  void destroy() override {
    delete[] bitmap_ptr_;
    bitmap_ptr_ = nullptr;
    blocks_.clear();
    blocks_bias_.clear();
  }

  auto query(BitPos key) -> bool override {
    if (key < blocks_bias_[0] || key > blocks_bias_.back()) {
      return false;
    }
    auto iter =
        std::upper_bound(blocks_bias_.begin(), blocks_bias_.end(), key) - 1;
    if (*iter == key) {
      return true;
    }

    size_t block_idx = iter - blocks_bias_.begin();

    return blocks_[block_idx].query(key - *iter);
  }

  auto query(BitPos left_key, BitPos right_key) -> bool override {
    if (right_key < blocks_bias_[0] || left_key > blocks_bias_.back()) {
      return false;
    }

    auto iter =
        std::upper_bound(blocks_bias_.begin(), blocks_bias_.end(), right_key);
    if (iter == blocks_bias_.end() || *(--iter) == right_key ||
        left_key <= *iter) {
      return true;
    }

    size_t block_idx = iter - blocks_bias_.begin();
    return blocks_[block_idx].query(left_key - *iter, right_key - *iter);
  }

  auto size() const -> size_t override {
    return sizeof(size_t) +                       /* nbatches */
           sizeof(size_t) +                       /* bitmap_sz_ */
           sizeof(uint16_t) * 2 +                 /* block_sz */
           blocks_bias_.size() * sizeof(BitPos) + /* blocks_bias_ size */
           bitmap_sz_;                            /* block_list_ */
  }

  auto serialize() const -> std::pair<uint8_t *, size_t> override {
    size_t size = sizeof(size_t) +      /* nbatches */
                  sizeof(size_t) +      /* bitmap_sz_ */
                  sizeof(uint16_t) * 2; /* block_sz */

    size_t nbatches = blocks_.size();
    size_t bias_list_sz = (nbatches + 1) * sizeof(BitPos);
    size += bitmap_sz_ + bias_list_sz;

    uint8_t *ser = new uint8_t[size];
    uint8_t *pos = ser;

    memcpy(pos, &nbatches, sizeof(size_t));
    pos += sizeof(size_t);

    memcpy(pos, &bitmap_sz_, sizeof(size_t));
    pos += sizeof(size_t);

    memcpy(pos, &block_sz_, sizeof(uint16_t));
    pos += sizeof(uint16_t);

    memcpy(pos, &last_block_sz_, sizeof(uint16_t));
    pos += sizeof(uint16_t);

    memcpy(pos, blocks_bias_.data(), bias_list_sz);
    pos += bias_list_sz;

    memcpy(pos, bitmap_ptr_, bitmap_sz_);

    return {ser, size};
  }

  static auto deserialize(uint8_t *ser) -> std::unique_ptr<BitMap> {
    size_t nbatches;
    memcpy(&nbatches, ser, sizeof(size_t));
    ser += sizeof(size_t);

    size_t bitmap_sz;
    memcpy(&bitmap_sz, ser, sizeof(size_t));
    ser += sizeof(size_t);

    uint16_t block_sz;
    memcpy(&block_sz, ser, sizeof(uint16_t));
    ser += sizeof(uint16_t);

    uint16_t last_block_sz;
    memcpy(&last_block_sz, ser, sizeof(uint16_t));
    ser += sizeof(uint16_t);

    BitsPos blocks_bias(nbatches + 1);
    size_t bias_sz = blocks_bias.size() * sizeof(BitPos);
    memcpy(blocks_bias.data(), ser, bias_sz);
    ser += bias_sz;

    uint8_t *bitmap_ptr = new uint8_t[bitmap_sz];
    memcpy(bitmap_ptr, ser, bitmap_sz);

    Blocks blocks;
    blocks.reserve(nbatches);

    uint8_t *pos = bitmap_ptr;
    for (size_t i = 0; i < nbatches; ++i) {
      blocks.emplace_back(i == nbatches - 1 ? last_block_sz : block_sz,
                          blocks_bias[i + 1] - blocks_bias[i], pos);
      pos += blocks.back().size();
    }

    return std::make_unique<BlockEliasFanoBitMap>(
        block_sz, bitmap_sz, last_block_sz, bitmap_ptr, blocks, blocks_bias);
  }

 private:
  uint16_t block_sz_;
  size_t bitmap_sz_ = 0;
  uint16_t last_block_sz_;

  uint8_t *bitmap_ptr_ = nullptr;
  Blocks blocks_;
  BitsPos blocks_bias_;
};

}  // namespace oasis