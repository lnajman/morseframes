cmake_minimum_required(VERSION 3.16)

if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin" AND NOT CMAKE_OSX_ARCHITECTURES)
    set(_MORSEFRAMES_PYTHON_FOR_ARCH "")
    if(DEFINED Python_EXECUTABLE)
        set(_MORSEFRAMES_PYTHON_FOR_ARCH "${Python_EXECUTABLE}")
    elseif(DEFINED PYTHON_EXECUTABLE)
        set(_MORSEFRAMES_PYTHON_FOR_ARCH "${PYTHON_EXECUTABLE}")
    endif()

    if(_MORSEFRAMES_PYTHON_FOR_ARCH)
        execute_process(
            COMMAND "${_MORSEFRAMES_PYTHON_FOR_ARCH}" -c
                    "import platform; print(platform.machine())"
            OUTPUT_VARIABLE _MORSEFRAMES_PYTHON_ARCH
            RESULT_VARIABLE _MORSEFRAMES_PYTHON_ARCH_RESULT
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(_MORSEFRAMES_PYTHON_ARCH_RESULT EQUAL 0
           AND (_MORSEFRAMES_PYTHON_ARCH STREQUAL "arm64"
                OR _MORSEFRAMES_PYTHON_ARCH STREQUAL "x86_64"))
            set(CMAKE_OSX_ARCHITECTURES "${_MORSEFRAMES_PYTHON_ARCH}"
                CACHE STRING "Build architectures for macOS" FORCE)
        endif()
    endif()
endif()

project(morseframes LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(DEFINED SKBUILD)
    set(MORSEFRAMES_DEFAULT_BUILD_PYTHON ON)
    set(MORSEFRAMES_DEFAULT_BUILD_TESTS OFF)
    set(MORSEFRAMES_DEFAULT_BUILD_EXAMPLES OFF)
    set(MORSEFRAMES_DEFAULT_BUILD_GUDHI_TOOLS OFF)
else()
    set(MORSEFRAMES_DEFAULT_BUILD_PYTHON OFF)
    set(MORSEFRAMES_DEFAULT_BUILD_TESTS ON)
    set(MORSEFRAMES_DEFAULT_BUILD_EXAMPLES ON)
    set(MORSEFRAMES_DEFAULT_BUILD_GUDHI_TOOLS ON)
endif()

option(MORSEFRAMES_BUILD_PYTHON
    "Build nanobind Python extension"
    ${MORSEFRAMES_DEFAULT_BUILD_PYTHON})
option(MORSEFRAMES_BUILD_TESTS
    "Build MorseFrames C++ tests"
    ${MORSEFRAMES_DEFAULT_BUILD_TESTS})
option(MORSEFRAMES_BUILD_EXAMPLES
    "Build MorseFrames C++ examples"
    ${MORSEFRAMES_DEFAULT_BUILD_EXAMPLES})
option(MORSEFRAMES_BUILD_GUDHI_TOOLS
    "Build optional GUDHI Simplex_tree integration tests and benchmark when headers are available"
    ${MORSEFRAMES_DEFAULT_BUILD_GUDHI_TOOLS})
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

if(MORSEFRAMES_BUILD_TESTS)
    add_executable(morseframes_tests
        tests/test_morse_reference.cpp
    )
    target_link_libraries(morseframes_tests PRIVATE morseframes)
endif()

if(MORSEFRAMES_BUILD_EXAMPLES)
    add_executable(morseframes_example_barcode
        examples/example_barcode.cpp
    )
    target_link_libraries(morseframes_example_barcode PRIVATE morseframes)
endif()

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

    add_executable(morseframes_gudhi_upstream_simplex_tree_tests
        upstream/gudhi/test/Morse_persistence/test_morse_persistence_simplex_tree.cpp
    )
    target_link_libraries(morseframes_gudhi_upstream_simplex_tree_tests PRIVATE morseframes)
    target_include_directories(morseframes_gudhi_upstream_simplex_tree_tests PRIVATE
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
    if(DEFINED SKBUILD)
        install(TARGETS _morse_core
            LIBRARY DESTINATION morseframes
            RUNTIME DESTINATION morseframes
        )
    endif()
endif()

enable_testing()
if(TARGET morseframes_tests)
    add_test(NAME morseframes_tests COMMAND morseframes_tests)
endif()
if(TARGET morseframes_gudhi_simplex_tree_view_tests)
    add_test(NAME morseframes_gudhi_simplex_tree_view_tests
        COMMAND morseframes_gudhi_simplex_tree_view_tests)
endif()
if(TARGET morseframes_gudhi_upstream_simplex_tree_tests)
    add_test(NAME morseframes_gudhi_upstream_simplex_tree_tests
        COMMAND morseframes_gudhi_upstream_simplex_tree_tests)
endif()
