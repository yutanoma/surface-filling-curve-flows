cmake_minimum_required(VERSION 3.10.0)

project(gc_polyscope_example)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Print the build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()
message(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")

### Configure the compiler
# This is a basic, decent setup that should do something sane on most compilers

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++17 -Wall -Wextra")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-unused-private-field")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3 -fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Setting gcc-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()
  
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)
else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()


# == Deps
add_subdirectory(deps/geometry-central)
add_subdirectory(deps/polyscope)

# Add libIGL
# With these lines commented out, we use libIGL in "header only" mode.  Uncomment to recurse and get the full package.
# (and see the libIGL build documentation)
message("\n\n == CMAKE recursively building libIGL\n")
# set(LIBIGL_WITH_EMBREE            OFF CACHE BOOL "Use Embree")
# set(LIBIGL_WITH_OPENGL            OFF CACHE BOOL "Use OpenGL")
# set(LIBIGL_WITH_OPENGL            OFF CACHE BOOL "Use OpenGL")
# set(LIBIGL_WITH_OPENGL_GLFW       OFF CACHE BOOL "Use GLFW"  )
# set(LIBIGL_WITH_OPENGL_GLFW_IMGUI OFF CACHE BOOL "Use ImGui" )
set(LIBIGL_RESTRICTED_TRIANGLE      ON CACHE BOOL "Use Triangle")
set(LIBIGL_COPYLEFT_TETGEN          ON CACHE BOOL "Use TetGen")
# note: the following line should come after setting triangle and tetgen
set(LIBIGL_PREDICATES               ON CACHE BOOL "Use exact predicates")
add_subdirectory(deps/libigl)

# read the modules
add_subdirectory(modules)

# tinyad
include(FetchContent)
FetchContent_Declare(
  tinyad
  GIT_REPOSITORY https://github.com/patr-schm/tinyad.git
  GIT_TAG 75093e14ef0d7bb39657c5f3b2aba1251afaa38c
)

FetchContent_GetProperties(tinyad)
if(NOT tinyad_POPULATED)
  # Fetch the content using previously declared details
  FetchContent_Populate(tinyad)
  message(STATUS "tinyad_SOURCE_DIR: ${tinyad_SOURCE_DIR}")
  message(STATUS "tinyad_BINARY_DIR: ${tinyad_BINARY_DIR}")
  add_subdirectory(${tinyad_SOURCE_DIR} ${tinyad_BINARY_DIR})
endif()

# == Build our project stuff

# 1. curve on surface
set(SRCS_curve_on_surface
  src/curve_on_surface.cpp
)

add_executable(curve_on_surface "${SRCS_curve_on_surface}")
target_include_directories(curve_on_surface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/modules/include)
target_include_directories(curve_on_surface PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# add the args.hxx project which we use for command line args
target_include_directories(curve_on_surface PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args")
target_include_directories(curve_on_surface PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl/include")
target_link_libraries(curve_on_surface PUBLIC geometry-central polyscope TinyAD modules igl::core igl_copyleft::tetgen igl_restricted::triangle igl::predicates)

# tinyad needs C++17 (not sure if that means this project does too)
target_compile_features(curve_on_surface PRIVATE cxx_std_17)

# 2. mesh deformation
set(SRCS_TEST_3
  src/mesh_deformation.cpp
)

add_executable(mesh_deformation "${SRCS_TEST_3}")
target_include_directories(mesh_deformation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/modules/include)
target_include_directories(mesh_deformation PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# add the args.hxx project which we use for command line args
target_include_directories(mesh_deformation PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/polyscope/deps/args")
target_include_directories(mesh_deformation PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl/include")
target_link_libraries(mesh_deformation PUBLIC geometry-central polyscope TinyAD modules igl::core igl_copyleft::tetgen igl_restricted::triangle igl::predicates)

# tinyad needs C++17 (not sure if that means this project does too)
target_compile_features(mesh_deformation PRIVATE cxx_std_17)
