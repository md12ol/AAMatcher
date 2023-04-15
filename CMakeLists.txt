cmake_minimum_required(VERSION 3.16.5)
project(AAMatcher)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE "Release")

add_executable(AAMatcher
        main.cpp
        main.h
        Bitsprayer.h
        Bitsprayer.cpp)
