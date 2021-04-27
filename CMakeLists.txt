cmake_minimum_required(VERSION 3.17)
project(glpkpp)

add_subdirectory(src)

link_libraries(glpkpp)
add_executable(glpktest main.cpp include/glpkpp.h)
