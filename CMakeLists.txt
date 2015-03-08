cmake_minimum_required(VERSION 2.8.4)
project(Matasano)

set(CMAKE_CXX_COMPILER clang++)
set(WARNINGS "-Wall -Wpedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNINGS} -O3 -std=c++1y -lc++abi -g -stdlib=libc++ -nostdinc++ -I/Users/ineol/Code/llvm/projects/libcxx/include -L/Users/ineol/Code/llvm/build/lib")

include_directories("/usr/local/include" "/usr/include" "/Users/ineol/Code/llvm/projects/libcxx/include")

set(SOURCE_FILES main.cc)
add_executable(Matasano ${SOURCE_FILES})
