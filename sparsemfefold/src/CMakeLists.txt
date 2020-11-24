# Mateo Gray 2020

cmake_minimum_required(VERSION 3.8)
project(src)

# enable c++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)



include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${CMAKE_SOURCE_DIR}/src")

set(sparsemfefold_SOURCE
    cmdline.cc
    trace_arrow.hh trace_arrow.cc 
    SparseMFEFold_1.cc
)



# create executables
add_executable(SparseMFEFold ${sparsemfefold_SOURCE})

# enable C++11
#target_compile_features(SparseMFEFold PRIVATE cxx_nullptr)

# link to simfold
target_link_libraries(SparseMFEFold LINK_PUBLIC RNA)


