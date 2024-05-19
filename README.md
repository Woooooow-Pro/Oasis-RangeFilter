# Oasis: An Optimal Disjoint Segmented Learned Range Filter [[paper](https://www.vldb.org/pvldb/vol17/p1911-luo.pdf)]

We introduce Oasis (<u>**O**</u>ptim<u>**a**</u>l Di<u>**s**</u>jo<u>**i**</u>nt <u>**S**</u>egmented Learned Range Filter), a novel learned range filter that divides the key space into disjointed intervals by excluding large empty ranges explicitly and optimally maps those unpruned intervals into a compressed bitmap. Besides, Oasis can identify its optimal configuration through a theoretical analysis. To enhance the versatility of Oasis, we further propose Oasis+, which integrates the design space of both learned and non-learned filters, delivering robust performance across a wide range of workloads. We evaluate the performance of both Oasis and Oasis+, when integrated into the key-value system RocksDB, using a diverse set of real-world and synthetic datasets and workloads. In RocksDB, Oasis and Oasis+ improve the performance by up to 1.4× and 6.2× when compared to state-of-the-art learned and non-learned range filters.

Guanduo Chen, Meng Li, Zhenying He, Siqiang Luo.


Note: Oasis+ requires the processor to support SIMD instructions.

## Directory Layout

- `oasis/` - contains all the files necessary to the implementation of the Oasis filter.

- `src/` - contains all the files necessary to the implementation of the OasisPlus filter.

- `benchmark/` - contains all the files used for standalone filter benchmarks.
  - `workloads/` - contains all the files used for generating synthetic workloads and downloading real-world datasets.`

## Dependency

- `jq` (to download Internet domains dataset)
- `CMake`, `gcc9`, `g++9`

```
# Linux
sudo apt-get install jq cmake
```

## Dataset

```
cd benchmark/workloads
./setup.sh
```

`books_800M_uint64` and `fb_200M_uint64` from https://github.com/learnedsystems/SOSD


## Standalone Filter Benchmarks

```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make all

cd ../benchmark
./in_mem_bench.sh test > test.txt&
```
