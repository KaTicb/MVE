cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

add_executable(lab1)
add_executable(lab2)

target_sources(lab1 PRIVATE src/main.cpp)
target_sources(lab2 PRIVATE src/main_1.cpp)