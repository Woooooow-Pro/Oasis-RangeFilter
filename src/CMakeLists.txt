add_subdirectory(learned_rf)
add_subdirectory(proteus)

add_library(
  OasisPlus
  STATIC
  oasis_plus.cc
  filter_builder.cc
  ${ALL_OBJECT_FILES}
  )

set(OasisPlus_LIBS
    oasis_plus_learned_rf
    oasis_plus_proteus
    )

target_link_libraries(
  OasisPlus
  PUBLIC
  ${HybridRF_LIBS}
)

target_include_directories(
    OasisPlus PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)