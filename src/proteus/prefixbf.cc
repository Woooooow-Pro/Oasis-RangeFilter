#include "proteus/prefixbf.h"

#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstring>
#include <functional>
#include <map>
#include <memory>
#include <random>

#include "proteus/clhash.h"
#include "proteus/MurmurHash3.h"

namespace oasis_plus {
PrefixBF::PrefixBF(uint32_t prefix_len, uint64_t nbits,
                   const std::vector<uint64_t>& keys)
    : prefix_len_(prefix_len),
      nmod_(std::min(static_cast<uint64_t>(UINT32_MAX), (nbits + 7) / 8 * 8)) {
  // Filter uses 32-bit MurmurHash function so the number of filter bits
  // cannot be more than UINT32_MAX. Such situations are unlikely to happen
  // in RocksDB as it would entail an extremely high BPK or number of keys.

  assert(nbits > 0);

  // Create zeroed out bit array - parentheses value initializes array to 0
  data_ = new uint8_t[div8(nmod_)]();
  assert(data_);

  std::vector<size_t> uniq_idxs;
  uniq_idxs.emplace_back(0);
  uint64_t prev_key = keys[0] >> (64 - prefix_len_);
  uint64_t current_key = 0;

  for (size_t i = 0; i < keys.size(); i++) {
    current_key = keys[i] >> (64 - prefix_len_);
    if (current_key != prev_key) {
      uniq_idxs.emplace_back(i);
      prev_key = current_key;
    }
  }

  uint32_t nhf = static_cast<uint32_t>(round(M_LN2 * nmod_ / uniq_idxs.size()));
  nhf = (nhf == 0 ? 1 : nhf);
  nhf = std::min(MAX_PBF_HASH_FUNCS, nhf);

  seeds32_.resize(nhf);
  std::mt19937 gen(1337);
  for (uint32_t i = 0; i < nhf; ++i) {
    seeds32_[i] = gen();
  }

  uint64_t edited_key = 0ULL;
  for (auto const& idx : uniq_idxs) {
    edited_key = keys[idx] >> (64 - prefix_len_);
    for (uint32_t i = 0; i < nhf; ++i) {
      set(hash(edited_key, seeds32_[i]), 1);
    }
  }
}

PrefixBF::PrefixBF(uint32_t prefix_len, uint64_t nbits,
                   const std::vector<std::string>& keys)
    : prefix_len_(prefix_len), nmod_((nbits + 7) / 8 * 8) {
  assert(nbits > 0);

  // Create zeroed out bit array - parentheses value initializes array to 0
  data_ = new uint8_t[(nmod_ / 8)]();
  assert(data_);

  std::vector<size_t> uniq_prefixes;
  uniq_prefixes.emplace_back(0);
  std::string prev_key = keys[0];
  std::string current_key = std::string();

  for (size_t i = 0; i < keys.size(); i++) {
    current_key = keys[i];
    if (compare(current_key, prev_key, prefix_len_) != 0) {
      uniq_prefixes.emplace_back(i);
      prev_key = current_key;
    }
  }

  uint32_t nhf =
      static_cast<uint32_t>(round(M_LN2 * nmod_ / uniq_prefixes.size()));
  nhf = (nhf == 0 ? 1 : nhf);
  nhf = std::min(MAX_PBF_HASH_FUNCS, nhf);

  seeds64_.resize(nhf);
  seeds128_.resize(nhf);

  std::mt19937_64 gen(1337);
  for (uint32_t i = 0; i < nhf; ++i) {
    seeds64_[i] = std::make_pair(gen(), gen());
    seeds128_[i] =
        get_random_key_for_clhash(seeds64_[i].first, seeds64_[i].second);
  }

  uint32_t prefix_byte_len = div8(prefix_len_ + 7);
  if (mod8(prefix_len_) == 0) {
    for (auto const& idx : uniq_prefixes) {
      for (uint32_t i = 0; i < nhf; ++i) {
        set(clhash(seeds128_[i], keys[idx].data(), prefix_byte_len) % nmod_, 1);
      }
    }
  } else {
    std::string edited_key;
    for (auto const& idx : uniq_prefixes) {
      edited_key = editKey(keys[idx], prefix_len_, true);
      for (uint32_t i = 0; i < nhf; ++i) {
        set(clhash(seeds128_[i], edited_key.data(), prefix_byte_len) % nmod_,
            1);
      }
    }
  }
}

PrefixBF::PrefixBF(uint32_t prefix_len, uint8_t* data,
                   std::vector<uint32_t> seeds32,
                   std::vector<std::pair<uint64_t, uint64_t>> seeds64,
                   uint64_t nmod)
    : prefix_len_(prefix_len),
      seeds32_(seeds32),
      seeds64_(seeds64),
      nmod_(nmod) {
  // The string hash function (CLHash) generates 128-bit seeds from a pair of
  // 64-bit integers.
  seeds128_.resize(seeds64_.size());
  for (size_t i = 0; i < seeds64_.size(); ++i) {
    seeds128_[i] =
        get_random_key_for_clhash(seeds64_[i].first, seeds64_[i].second);
  }

  // Copy over bit array
  data_ = new uint8_t[(nmod_ / 8)];
  assert(data_);
  memcpy(reinterpret_cast<void*>(data_), data, (nmod_ / 8));
}

auto PrefixBF::hash(const uint64_t edited_key, const uint32_t& seed)
    -> uint64_t {
  uint32_t h;
  murmur3::MurmurHash3_x86_32(&edited_key, 8, seed, &h);
  return static_cast<uint64_t>(h) % nmod_;
}

auto PrefixBF::get(uint64_t i) const -> bool {
  return (data_[div8(i)] >> (7 - mod8(i))) & 1;
}

void PrefixBF::set(uint64_t i, bool v) {
  if (get(i) != v) {
    data_[div8(i)] ^= (1 << (7 - mod8(i)));
  }
}

auto PrefixBF::Query(const uint64_t key, bool shift) -> bool {
  bool out = true;
  uint64_t k = shift ? key >> (64 - prefix_len_) : key;
  for (size_t i = 0; i < seeds32_.size() && out; ++i) {
    out &= get(hash(k, seeds32_[i]));
  }
  return out;
}

/*
    To execute a range query, we shift the query bounds to the
    specified prefix length and do point queries for all the
    values in the shifted query range.
*/
auto PrefixBF::Query(const uint64_t from, const uint64_t to) -> bool {
  uint64_t iterator = from >> (64 - prefix_len_);
  uint64_t upper_bound = (to - 1) >> (64 - prefix_len_);
  while (iterator <= upper_bound) {
    if (Query(iterator, false)) {
      return true;
    }
    iterator++;
  }
  return false;
}

auto PrefixBF::Query(const std::string& key) -> bool {
  bool out = true;
  uint32_t prefix_byte_len = div8(prefix_len_ + 7);

  std::string padded_key = key;
  if (key.length() < prefix_byte_len) {
    padded_key.resize(prefix_byte_len, '\0');
  }

  if (mod8(prefix_len_) == 0) {
    for (size_t i = 0; i < seeds128_.size() /* # of hash functions */ && out;
         ++i) {
      out &=
          get(clhash(seeds128_[i], padded_key.data(), prefix_byte_len) % nmod_);
    }
  } else {
    std::string edited_key = editKey(padded_key, prefix_len_, true);
    for (size_t i = 0; i < seeds128_.size() /* # of hash functions */ && out;
         ++i) {
      out &=
          get(clhash(seeds128_[i], edited_key.data(), prefix_byte_len) % nmod_);
    }
  }

  return out;
}

/*
    For string range queries, we choose to make both `from` and `to` inclusive
   in order to avoid generating (to - 1). When integrated in RocksDB, this may
   incur an additional prefix query in some (rare) situations. We also
   pre-calculate the number of prefix queries to avoid multiple expensive string
   comparisons. Similar to the integer range query, we create a query iterator
   and increment it successively.
*/
auto PrefixBF::Query(const std::string& from, const std::string& to) -> bool {
  uint32_t prefix_byte_len = div8(prefix_len_ + 7);
  uint32_t shift_bits = mod8(8 - mod8(prefix_len_));

  // Pad `to` with zeroes because a shorter string is
  // lexicographically smaller than longer strings
  // that share the same prefix
  std::string padded_from = editKey(from, prefix_len_, true);
  std::string padded_to = editKey(to, prefix_len_, true);

  // count_prefixes returns 0 if the result is bigger than a uint64_t - treat as
  // guaranteed false positive
  uint64_t total_queries = count_prefixes(padded_from, padded_to, prefix_len_);
  if (total_queries == 0) {
    return true;
  }

  bool carry = false;
  uint32_t idx = 0;
  uint8_t shifted_last_char = 0;
  std::string& iterator = padded_from;

  for (uint64_t i = 0; i < total_queries; i++) {
    if (Query(iterator)) {
      return true;
    }

    idx = prefix_byte_len - 1;

    // prefix length may fall within a byte so we shift the byte accordingly
    shifted_last_char = static_cast<uint8_t>(iterator[idx]) >> shift_bits;
    if (shifted_last_char == ((MAX_UINT8) >> shift_bits)) {
      carry = true;
      iterator[idx] = 0;
    } else {
      carry = false;
      iterator[idx] = (shifted_last_char + 1) << shift_bits;
    }

    // increment prior bytes if there is a carry
    while (carry && idx > 0) {
      idx -= 1;
      if (static_cast<uint8_t>(iterator[idx]) == MAX_UINT8) {
        iterator[idx] = 0;
      } else {
        iterator[idx] = static_cast<uint8_t>(iterator[idx]) + 1;
        carry = false;
      }
    }
  }

  return false;
}

auto PrefixBF::serialize() const -> std::pair<uint8_t*, uint64_t> {
  // Size (bytes) of the serialized Prefix BF.
  uint64_t serlen =
      sizeof(uint32_t)   /* prefix_len_ */
      + sizeof(uint64_t) /* nmod_ */
      + sizeof(size_t) /* size of seeds32_ */ +
      seeds32_.size() * sizeof(uint32_t) /* actual seeds32_ vector */
      + sizeof(size_t) /* size of seeds64_ */ +
      seeds64_.size() *
          sizeof(std::pair<uint64_t, uint64_t>) /* actual seeds64_ vector */
      + div8(nmod_);

  uint8_t* out = new uint8_t[serlen];
  uint8_t* pos = out;

  memcpy(pos, &prefix_len_, sizeof(uint32_t));
  pos += sizeof(uint32_t);

  memcpy(pos, &nmod_, sizeof(uint64_t));
  pos += sizeof(uint64_t);

  // seeds array for integer prefixbf
  size_t seeds_sz = seeds32_.size();
  memcpy(pos, &seeds_sz, sizeof(size_t));
  pos += sizeof(size_t);

  memcpy(pos, seeds32_.data(), seeds32_.size() * sizeof(uint32_t));
  pos += seeds32_.size() * sizeof(uint32_t);

  // seeds array for string prefixbf
  size_t seeds64_sz = seeds64_.size();
  memcpy(pos, &seeds64_sz, sizeof(size_t));
  pos += sizeof(size_t);

  memcpy(pos, seeds64_.data(),
         seeds64_.size() * sizeof(std::pair<uint64_t, uint64_t>));
  pos += seeds64_.size() * sizeof(std::pair<uint64_t, uint64_t>);

  memcpy(pos, data_, div8(nmod_));
  return {out, serlen};
}

auto PrefixBF::deserialize(uint8_t* ser) -> std::pair<PrefixBF*, uint64_t> {
  uint8_t* pos = ser;

  uint32_t prefix_len;
  memcpy(&prefix_len, pos, sizeof(uint32_t));
  pos += sizeof(uint32_t);

  uint64_t nmod;
  memcpy(&nmod, pos, sizeof(uint64_t));
  pos += sizeof(uint64_t);

  size_t seeds32_sz;
  memcpy(&seeds32_sz, pos, sizeof(size_t));
  pos += sizeof(size_t);

  uint32_t* pos32 = reinterpret_cast<uint32_t*>(pos);
  std::vector<uint32_t> seeds32(pos32, pos32 + seeds32_sz);
  pos += seeds32_sz * sizeof(uint32_t);

  size_t seeds64_sz;
  memcpy(&seeds64_sz, pos, sizeof(size_t));
  pos += sizeof(size_t);

  std::pair<uint64_t, uint64_t>* pos64 =
      reinterpret_cast<std::pair<uint64_t, uint64_t>*>(pos);
  std::vector<std::pair<uint64_t, uint64_t>> seeds64(pos64, pos64 + seeds64_sz);
  pos += seeds64_sz * sizeof(std::pair<uint64_t, uint64_t>);

  size_t len = pos - ser + div8(nmod);
  return {new PrefixBF(prefix_len, pos, seeds32, seeds64, nmod), len};
}

}  // namespace oasis_plus