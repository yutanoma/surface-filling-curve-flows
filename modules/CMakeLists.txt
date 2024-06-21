cmake_minimum_required(VERSION 3.10.0)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

# Collect all .cpp files in the subdirectories
file(GLOB_RECURSE MODULES_SRCS src/*.cpp src/*.cc)

# Create a single library for the project
add_library(modules ${MODULES_SRCS})

target_link_libraries(modules geometry-central TinyAD igl_copyleft::tetgen igl_restricted::triangle igl::predicates)
target_include_directories(modules PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../deps/libigl/include")

target_include_directories(modules PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# tinyad needs C++17 (not sure if that means this project does too)
target_compile_features(modules PRIVATE cxx_std_17)
