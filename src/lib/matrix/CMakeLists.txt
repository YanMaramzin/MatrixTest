cmake_minimum_required(VERSION 3.22)
project(matrix)

add_library(${PROJECT_NAME}
        matrix.h
        matrix.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
        "${PROJECT_BINARY_DIR}"
        "${PROJECT_SOURCE_DIR}"
)
