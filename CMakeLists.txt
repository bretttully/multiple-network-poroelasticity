cmake_minimum_required(VERSION 2.8)

# ====================================================================
# use clang instead of gcc? (this must come before declaring the project)
option(USE_CLANG "Build application with clang" ON)
if(USE_CLANG)
    SET (CMAKE_C_COMPILER   "/usr/bin/clang")
    SET (CMAKE_CXX_COMPILER "/usr/bin/clang++")
    SET (CMAKE_AR      "/usr/bin/llvm-ar")
    SET (CMAKE_LINKER  "/usr/bin/llvm-ld")
    SET (CMAKE_NM      "/usr/bin/llvm-nm")
    SET (CMAKE_OBJDUMP "/usr/bin/llvm-objdump")
    SET (CMAKE_RANLIB  "/usr/bin/llvm-ranlib")
else(USE_CLANG)
    SET (CMAKE_C_COMPILER   "/usr/bin/gcc")
    SET (CMAKE_CXX_COMPILER "/usr/bin/g++")
    SET (CMAKE_AR      "/usr/bin/ar")
    SET (CMAKE_LINKER  "/usr/bin/ld")
    SET (CMAKE_NM      "/usr/bin/nm")
    SET (CMAKE_OBJDUMP "/usr/bin/objdump")
    SET (CMAKE_RANLIB  "/usr/bin/ranlib")
endif(USE_CLANG)

# ====================================================================
# The project must be defined after selecting the compiler
project(mpet)

# ====================================================================
# use these macros to remove any flags
macro(remove_cc_flag flag)
    string(REPLACE "${flag}" "" CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS}")
endmacro()
macro(remove_cxx_flag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endmacro()
macro(remove_c_flag flag)
    string(REPLACE "${flag}" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
endmacro()

# ====================================================================
# default build type
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to RelWithDebInfo")
    set(CMAKE_BUILD_TYPE "Release")
#    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

SET (CMAKE_CXX_FLAGS                "-Wall -std=c++11")
SET (CMAKE_CXX_FLAGS_DEBUG          "-g")
SET (CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
# -----------------------------------------------------

# ====================================================================
# find the external packages
find_package(PkgConfig REQUIRED)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

find_package(PythonLibs 2.7 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
link_directories(${PYTHON_LIBRARY_DIRS})

# ----
# for the threaded testing
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)

# ====================================================================
# set up the default directories
set(base_dir ${CMAKE_SOURCE_DIR})
set(CMAKE_INSTALL_PREFIX ${base_dir})

set(src_dir ${base_dir}/src)
include_directories(${src_dir})

include_directories(${src_dir}/tnt) # used for the JAMA library

set(lib_dir ${base_dir}/lib)
link_directories(${lib_dir})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${lib_dir})

set(bin_dir ${base_dir}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${bin_dir})

# ====================================================================
# set up the libraries
set(GTEST_LIBS gtest ${CMAKE_THREAD_LIBS_INIT})
set(BOOST_LIBS boost_filesystem boost_system boost_thread boost_iostreams)

# ====================================================================
# set up the git version
# see: http://stackoverflow.com/questions/1435953/how-can-i-pass-git-sha1-to-compiler-as-definition-using-cmake
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# ====================================================================
# add the src and testing directories
add_subdirectory(src)

message(STATUS "-- -- -- User set options -- -- -- --")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Use Clang: ${USE_CLANG}")
message(STATUS "Report stack traces: ${REPORT_EXCEPTION_STACKTRACE}")

