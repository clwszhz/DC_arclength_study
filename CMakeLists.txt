cmake_minimum_required(VERSION 3.16)

project(dcsolve LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(dcsolve_core
    src/Circuit.cpp
    src/Homotopy.cpp
    src/Lup.cpp
    src/Mna.cpp
    src/Nonlinear.cpp
    src/Parser.cpp
)

target_include_directories(dcsolve_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(MSVC)
    target_compile_options(dcsolve_core PRIVATE /W4)
else()
    target_compile_options(dcsolve_core PRIVATE -Wall -Wextra -pedantic)
endif()

add_executable(dcsolve src/main.cpp)
target_link_libraries(dcsolve PRIVATE dcsolve_core)

add_executable(homotopy_solve src/homotopy_main.cpp)
target_link_libraries(homotopy_solve PRIVATE dcsolve_core)

add_executable(dcsolve_tests tests/test_lup.cpp)
target_link_libraries(dcsolve_tests PRIVATE dcsolve_core)

enable_testing()
add_test(NAME dcsolve_tests COMMAND dcsolve_tests)
