#pragma once

#include <memory>
#include <vector>

#include "bitmap.h"
#include "bitset.h"
#include "common/type.h"
#include "models.h"

namespace oasis {

template <typename Key>
class Oasis {
  using Keys = std::vector<Key>;

  using OasisModels = CDFModels<Key>;
  using OasisBitMap = BlockEliasFanoBitMap;

  using ModelsPtr = std::unique_ptr<Models<Key>>;
  using BitMapPtr = std::unique_ptr<BitMap>;

 public:
  Oasis(double bpk, uint16_t block_sz, const Keys &keys)
      : models_(std::make_unique<OasisModels>(bpk, block_sz, keys)),
        bitmap_(std::make_unique<OasisBitMap>(block_sz)) {
    bitmap_->build(models_->get_locations(keys));

    if (double bpk_strict = bpk - size() * 8.0 / keys.size();
        bpk_strict < 0.2) {
      return;
    } else {
      bpk += bpk_strict;
      // delete the pointers
      ModelsPtr tmp_1 = std::move(models_);
    }
    models_ = std::make_unique<OasisModels>(bpk, block_sz, keys);

    bitmap_->destroy();
    bitmap_->build(models_->get_locations(keys));
  }

  Oasis(ModelsPtr &models, BitMapPtr &bitmap)
      : models_(std::move(models)), bitmap_(std::move(bitmap)) {}

  ~Oasis() = default;

  auto query(const Key &key) -> bool {
    BitPos bit_pos;
    auto status = models_->query(key, bit_pos);

    switch (status) {
      case QueryPosStatus::EXIST:
        return true;
      case QueryPosStatus::OUT_OF_SCOPE:
        return false;
      case QueryPosStatus::NO_IDEA:
        return bitmap_->query(bit_pos);
      default:
        return false;
    }
  }

  auto query(const Key &left_key, const Key &right_key) -> bool {
    std::pair<BitPos, BitPos> bit_pos;
    QueryPosStatus status = models_->query(left_key, right_key, bit_pos);
    switch (status) {
      case QueryPosStatus::EXIST:
        return true;
      case QueryPosStatus::OUT_OF_SCOPE:
        return false;
      default:
        return bitmap_->query(bit_pos.first, bit_pos.second);
    }
  }

  auto size() const -> size_t { return models_->size() + bitmap_->size(); }

  auto serialize() const -> std::pair<uint8_t *, size_t> {
    size_t oasis_sz = size();
    uint8_t *ser = new uint8_t[oasis_sz];
    uint8_t *pos = ser;

    auto ser_model = models_->serialize().first;
    memcpy(pos, ser_model.first, sizeof(size_t));
    delete[] ser_model.first;
    pos += ser_model.second;

    uint8_t *ser_bitmap = bitmap_->serialize().first;
    memcpy(pos, ser_bitmap, sizeof(size_t));
    delete[] ser_bitmap;

    return {ser, oasis_sz};
  }

  static auto deserialize(uint8_t *ser) -> std::unique_ptr<Oasis> {
    ModelsPtr models = OasisModels::deserialize(ser);
    ser += models->size();
    BitMapPtr bitmap = OasisBitMap::deserialize(ser);
    return std::make_unique<Oasis>(std::move(models), std::move(bitmap));
  }

 private:
  ModelsPtr models_ = nullptr;
  BitMapPtr bitmap_ = nullptr;
};

}  // namespace oasis