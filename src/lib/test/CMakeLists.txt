cmake_minimum_required(VERSION 3.22)
project(testi)

find_package(Qt5 REQUIRED COMPONENTS Test)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOMOC ON)

enable_testing(true)

add_executable(testi test_matrix.h test_matrix.cpp)
add_test(NAME testi COMMAND testi)

target_link_libraries(${PROJECT_NAME} PUBLIC
    matrix
    PRIVATE
    Qt5::Test)

