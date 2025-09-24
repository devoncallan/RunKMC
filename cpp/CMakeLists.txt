cmake_minimum_required(VERSION 3.15)
project(RunKMC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Version configuration - can be overridden by build system
if(NOT DEFINED RUNKMC_VERSION)
    set(RUNKMC_VERSION "0.1.0-dev")
    message(WARNING "RUNKMC_VERSION not provided, using dev version")
endif()

# Parse version components
string(REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)" _ ${RUNKMC_VERSION})
set(RUNKMC_VERSION_MAJOR ${CMAKE_MATCH_1})
set(RUNKMC_VERSION_MINOR ${CMAKE_MATCH_2})
set(RUNKMC_VERSION_PATCH ${CMAKE_MATCH_3})

# Generate version header
configure_file(
    "${CMAKE_SOURCE_DIR}/include/runkmc/version.h.in"
    "${CMAKE_BINARY_DIR}/include/runkmc/version.h"
    @ONLY
)

# Disable tests and examples for dependencies
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)

# Fetch dependencies
include(FetchContent)

# Fetch Eigen
FetchContent_Declare(
  Eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
  GIT_SHALLOW TRUE
)

# Fetch yaml-cpp
FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG 0.8.0  # Latest stable version
  GIT_SHALLOW TRUE
)

FetchContent_MakeAvailable(Eigen3 yaml-cpp)

# Include directories
include_directories(include/runkmc)
include_directories(${CMAKE_BINARY_DIR}/include/runkmc)  # For generated version.h

add_executable(RunKMC src/RunKMC.cpp)

# Link libraries
target_link_libraries(RunKMC Eigen3::Eigen yaml-cpp::yaml-cpp)
target_compile_options(RunKMC PRIVATE -O3)
