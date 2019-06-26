cmake_minimum_required (VERSION 2.6)
set (CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

project (connectivity)
add_executable(connectivity cnct.cpp geometry.cpp tools/kdtree.cpp)