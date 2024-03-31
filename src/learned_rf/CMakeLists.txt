add_library(
  oasis_plus_learned_rf
  OBJECT
  learned_rf.cc
  bitset.cc
)

set(ALL_OBJECT_FILES
    ${ALL_OBJECT_FILES} $<TARGET_OBJECTS:oasis_plus_learned_rf>
    PARENT_SCOPE)