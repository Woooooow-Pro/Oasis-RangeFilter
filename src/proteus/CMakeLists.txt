add_library(
  oasis_plus_proteus
  OBJECT
  proteus.cc
  bitvector.cc
  config.cc
  label_vector.cc
  louds_dense.cc
  louds_sparse.cc
  modeling.cc
  prefixbf.cc
  rank.cc
  select.cc
  suffix.cc
  surf_builder.cc
  clhash.cc
  MurmurHash3.cc
)

target_link_libraries(oasis_plus_proteus clhash murmur3)

set(ALL_OBJECT_FILES
    ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:oasis_plus_proteus>
    PARENT_SCOPE)