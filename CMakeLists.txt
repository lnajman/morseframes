cmake_minimum_required(VERSION 3.16)

project(morseframes LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(MORSEFRAMES_BUILD_PYTHON "Build nanobind Python extension" OFF)
option(MORSEFRAMES_BUILD_GUDHI_TOOLS
    "Build optional GUDHI Simplex_tree integration tests and benchmark when headers are available"
    ON)
set(MORSEFRAMES_GUDHI_INCLUDE_DIR "" CACHE PATH
    "Optional path to a GUDHI include directory containing gudhi/Simplex_tree.h")
set(MORSEFRAMES_BOOST_INCLUDE_DIR "" CACHE PATH
    "Optional path to a Boost include directory needed by GUDHI headers")

if(MORSEFRAMES_BUILD_GUDHI_TOOLS)
    if(NOT MORSEFRAMES_GUDHI_INCLUDE_DIR)
        find_path(MORSEFRAMES_GUDHI_INCLUDE_DIR
            NAMES gudhi/Simplex_tree.h
            HINTS
                ENV GUDHI_ROOT
                ENV GUDHI_DIR
            PATH_SUFFIXES include
            DOC "Path to a GUDHI include directory containing gudhi/Simplex_tree.h")
    endif()

    if(NOT MORSEFRAMES_BOOST_INCLUDE_DIR)
        find_path(MORSEFRAMES_BOOST_INCLUDE_DIR
            NAMES boost/core/empty_value.hpp
            HINTS
                ENV BOOST_ROOT
                ENV Boost_ROOT
                ENV CONDA_PREFIX
            PATH_SUFFIXES include
            DOC "Path to a Boost include directory needed by GUDHI headers")
    endif()
endif()

add_library(morseframes INTERFACE)
target_include_directories(morseframes INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_executable(morseframes_tests
    tests/test_morse_reference.cpp
)
target_link_libraries(morseframes_tests PRIVATE morseframes)

add_executable(morseframes_example_barcode
    examples/example_barcode.cpp
)
target_link_libraries(morseframes_example_barcode PRIVATE morseframes)

if(MORSEFRAMES_BUILD_GUDHI_TOOLS
   AND MORSEFRAMES_GUDHI_INCLUDE_DIR
   AND EXISTS "${MORSEFRAMES_GUDHI_INCLUDE_DIR}/gudhi/Simplex_tree.h"
   AND MORSEFRAMES_BOOST_INCLUDE_DIR
   AND EXISTS "${MORSEFRAMES_BOOST_INCLUDE_DIR}/boost/core/empty_value.hpp")
    add_executable(morseframes_gudhi_simplex_tree_view_tests
        tests/test_gudhi_simplex_tree_view.cpp
    )
    target_link_libraries(morseframes_gudhi_simplex_tree_view_tests PRIVATE morseframes)
    target_include_directories(morseframes_gudhi_simplex_tree_view_tests PRIVATE
        "${MORSEFRAMES_GUDHI_INCLUDE_DIR}"
        "${MORSEFRAMES_BOOST_INCLUDE_DIR}"
    )

    add_executable(morseframes_benchmark_gudhi_view
        benchmarks/benchmark_gudhi_view.cpp
    )
    target_link_libraries(morseframes_benchmark_gudhi_view PRIVATE morseframes)
    target_include_directories(morseframes_benchmark_gudhi_view PRIVATE
        "${MORSEFRAMES_GUDHI_INCLUDE_DIR}"
        "${MORSEFRAMES_BOOST_INCLUDE_DIR}"
    )

    add_executable(morseframes_gudhi_simplex_tree_example
        examples/gudhi_simplex_tree_morse.cpp
    )
    target_link_libraries(morseframes_gudhi_simplex_tree_example PRIVATE morseframes)
    target_include_directories(morseframes_gudhi_simplex_tree_example PRIVATE
        "${MORSEFRAMES_GUDHI_INCLUDE_DIR}"
        "${MORSEFRAMES_BOOST_INCLUDE_DIR}"
    )

    add_executable(morseframes_gudhi_style_simplex_tree_example
        examples/example_morse_persistence_from_simplex_tree.cpp
    )
    target_link_libraries(morseframes_gudhi_style_simplex_tree_example PRIVATE morseframes)
    target_include_directories(morseframes_gudhi_style_simplex_tree_example PRIVATE
        "${MORSEFRAMES_GUDHI_INCLUDE_DIR}"
        "${MORSEFRAMES_BOOST_INCLUDE_DIR}"
    )
elseif(MORSEFRAMES_BUILD_GUDHI_TOOLS
       AND MORSEFRAMES_GUDHI_INCLUDE_DIR
       AND EXISTS "${MORSEFRAMES_GUDHI_INCLUDE_DIR}/gudhi/Simplex_tree.h")
    message(STATUS
        "Skipping GUDHI Simplex_tree tools because Boost headers were not found. "
        "Set MORSEFRAMES_BOOST_INCLUDE_DIR or BOOST_ROOT to enable them.")
elseif(MORSEFRAMES_BUILD_GUDHI_TOOLS)
    message(STATUS
        "Skipping GUDHI Simplex_tree tools because GUDHI headers were not found. "
        "Set MORSEFRAMES_GUDHI_INCLUDE_DIR or GUDHI_ROOT to enable them.")
endif()

if(MORSEFRAMES_BUILD_PYTHON)
    find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE
        OUTPUT_VARIABLE nanobind_ROOT
    )
    find_package(nanobind CONFIG REQUIRED)

    nanobind_add_module(_morse_core
        python/morseframes/_core_bindings.cpp
    )
    target_link_libraries(_morse_core PRIVATE morseframes)
    set_target_properties(_morse_core PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python/morseframes"
    )
endif()

enable_testing()
add_test(NAME morseframes_tests COMMAND morseframes_tests)
if(TARGET morseframes_gudhi_simplex_tree_view_tests)
    add_test(NAME morseframes_gudhi_simplex_tree_view_tests
        COMMAND morseframes_gudhi_simplex_tree_view_tests)
endif()
