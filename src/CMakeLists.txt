cmake_minimum_required (VERSION 2.6)
project (regply)

add_subdirectory(tinyply)

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

set (regply_VERSION_MAJOR 1)
set (regply_VERSION_MINOR 0)
set (regply_VERSION_PATCH 1)
set (regply_GIT_BRANCH "\"${GIT_BRANCH}\"")
set (regply_GIT_COMMIT "\"${GIT_COMMIT}\"")

configure_file (
  "${PROJECT_SOURCE_DIR}/regply.h.in"
  "${PROJECT_BINARY_DIR}/regply.h"
  )

include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/tinyply/source")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(PCL 1.8 REQUIRED COMPONENTS common io)
include_directories(${PCL_INCLUDE_DIRS} ../cloudcompare/CC/include ../cloudcompare/libs/qCC_db)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
link_directories(/usr/lib/cloudcompare)

add_executable(regply regply.cpp ply_io.cpp ../cloudcompare/CC/src/BoundingBox.cpp ../cloudcompare/CC/src/ScalarField.cpp)

# ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES} adds too much unused symbols
target_link_libraries(regply tinyply CC_CORE_LIB)

install(TARGETS regply DESTINATION /usr/bin)
