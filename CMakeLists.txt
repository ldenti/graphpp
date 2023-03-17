cmake_minimum_required(VERSION 3.14)

include(ExternalProject)
include(FetchContent)

#set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(GRAPHPP VERSION 1.0.0)
add_executable(graphpp main.cpp)
set_target_properties(graphpp PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")

# sdsl-lite
############
message(STATUS "sdsl-lite will be built from source")
ExternalProject_Add(sdsl
  GIT_REPOSITORY https://github.com/vgteam/sdsl-lite.git
  GIT_TAG caa2e60f7dcc2c85578eb85a89b88b485462eee4
  BUILD_IN_SOURCE 1
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND bash install.sh .
  )
ExternalProject_Get_Property(sdsl SOURCE_DIR)
SET(SDSL_SOURCE_DIR ${SOURCE_DIR})
SET(SDSL_INCLUDE_DIR ${SDSL_SOURCE_DIR}/include)
add_library(SDSL STATIC IMPORTED)
set_target_properties(SDSL PROPERTIES IMPORTED_LOCATION ${SDSL_SOURCE_DIR}/lib/libsdsl.a)

# gcsa2
########
message(STATUS "gcsa2 will be built from source")
ExternalProject_Add(gcsa
  GIT_REPOSITORY https://github.com/jltsiren/gcsa2.git
  GIT_TAG 6b4e1dd9af2c6f89e1164b3e20d3f4f36767126c
  UPDATE_COMMAND ""
  BUILD_IN_SOURCE 1
  PATCH_COMMAND sed -i "s@../sdsl-lite@${SDSL_SOURCE_DIR}@g" Makefile
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make
  INSTALL_COMMAND ""
  )
ExternalProject_Get_Property(gcsa SOURCE_DIR)
SET(GCSA_SOURCE_DIR ${SOURCE_DIR})
SET(GCSA_INCLUDE_DIR ${GCSA_SOURCE_DIR}/include)
add_library(GCSA STATIC IMPORTED)
set_target_properties(GCSA PROPERTIES IMPORTED_LOCATION ${GCSA_SOURCE_DIR}/lib/libgcsa2.a)

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

############################################
### final setup (includes and libraries) ###
############################################

add_dependencies(gcsa sdsl)
add_dependencies(graphpp gcsa sdsl)

target_include_directories(graphpp
  PRIVATE ${SDSL_INCLUDE_DIR}
  PRIVATE ${GCSA_INCLUDE_DIR}
)

target_link_libraries(graphpp
  PUBLIC GCSA
  PUBLIC SDSL
  PUBLIC z
  PUBLIC pthread
  PUBLIC OpenMP::OpenMP_CXX
)