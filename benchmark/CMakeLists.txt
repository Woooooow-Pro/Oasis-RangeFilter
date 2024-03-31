add_executable(in_mem_bench in_mem_bench.cc)

set(IN_MEM_BENCH_LIBS
    OasisPlus
)
target_link_libraries(in_mem_bench ${IN_MEM_BENCH_LIBS})

add_subdirectory(workloads)