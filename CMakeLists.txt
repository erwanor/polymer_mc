cmake_minimum_required(VERSION 2.8)
set(CMAKE_VERBOSE_MAKEFILE 1)

file(GLOB c_srcs ./src/*.c)
include_directories(./include)

add_definitions(-DSIMULATION_3D)

# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")

add_executable(polymer_mc ${c_srcs})
target_link_libraries(polymer_mc m)
