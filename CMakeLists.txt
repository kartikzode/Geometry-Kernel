cmake_minimum_required(VERSION 3.20)
project(cadkernel VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(CAD_SANITIZE "Build with ASan + UBSan" OFF)
option(CAD_WERROR "Treat warnings as errors" OFF)

add_library(cad_warnings INTERFACE)
if(MSVC)
  target_compile_options(cad_warnings INTERFACE /W4 $<$<BOOL:${CAD_WERROR}>:/WX>)
else()
  target_compile_options(cad_warnings INTERFACE -Wall -Wextra -Wpedantic -Wconversion -Wshadow
    $<$<BOOL:${CAD_WERROR}>:-Werror>)
  if(CAD_SANITIZE)
    target_compile_options(cad_warnings INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(cad_warnings INTERFACE -fsanitize=address,undefined)
  endif()
endif()

add_library(geom src/Matrix4.cpp)
target_include_directories(geom PUBLIC include)
target_link_libraries(geom PRIVATE cad_warnings)

include(CTest)
if(BUILD_TESTING)
  include(FetchContent)
  FetchContent_Declare(googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  enable_testing()
  add_executable(geom_tests
    tests/test_tolerance.cpp
    tests/test_vec3.cpp
    tests/test_matrix4.cpp
    tests/test_quaternion.cpp
    tests/test_bbox.cpp)
  target_link_libraries(geom_tests PRIVATE geom cad_warnings GTest::gtest_main)
  include(GoogleTest)
  gtest_discover_tests(geom_tests)
endif()
