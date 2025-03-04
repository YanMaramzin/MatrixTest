cmake_minimum_required(VERSION 3.5)
project(MatrixTest)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(lib)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME}
        matrix)
