cmake_minimum_required(VERSION 3.20)

project(step_detection LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

add_library(step_algorithms STATIC
    algorithms/peak_detector.c
    algorithms/threshold_detector.c
)

target_include_directories(step_algorithms PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/algorithms
)

add_executable(step_test
    host/main.c
    host/csv_reader.c
)

target_include_directories(step_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/host
)

target_link_libraries(step_test PRIVATE step_algorithms)

if(NOT MSVC)
    target_link_libraries(step_algorithms PUBLIC m)
endif()

if(MSVC)
    target_compile_options(step_algorithms PRIVATE /W4)
    target_compile_options(step_test PRIVATE /W4)
else()
    target_compile_options(step_algorithms PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(step_test PRIVATE -Wall -Wextra -Wpedantic)
endif()
