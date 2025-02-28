cmake_minimum_required(VERSION 3.30)
project(MatrixTest)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(lib)
add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME}
        matrix)