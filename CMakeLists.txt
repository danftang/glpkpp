cmake_minimum_required(VERSION 3.17)
project(glpkpp)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src)

link_libraries(glpkpp)
add_executable(glpktest main.cpp include/glpkpp.h)
