cmake_minimum_required(VERSION 3.21)

project(morseframes_ttk_gradient_benchmark LANGUAGES C CXX)

find_package(Threads REQUIRED)

if(NOT DEFINED TTK_SOURCE_DIR)
  message(FATAL_ERROR "Set TTK_SOURCE_DIR to the pinned TTK source checkout")
endif()

set(TTK_BUILD_VTK_WRAPPERS OFF CACHE BOOL "" FORCE)
set(TTK_BUILD_PARAVIEW_PLUGINS OFF CACHE BOOL "" FORCE)
set(TTK_BUILD_STANDALONE_APPS OFF CACHE BOOL "" FORCE)
set(TTK_BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)
set(TTK_ENABLE_MPI OFF CACHE BOOL "" FORCE)
set(TTK_ENABLE_TORCH OFF CACHE BOOL "" FORCE)
set(TTK_ENABLE_CGAL OFF CACHE BOOL "" FORCE)
set(TTK_ENABLE_EMBREE OFF CACHE BOOL "" FORCE)
set(TTK_ENABLE_GRAPHVIZ OFF CACHE BOOL "" FORCE)
set(TTK_ENABLE_SQLITE3 OFF CACHE BOOL "" FORCE)
set(TTK_ENABLE_TBB OFF CACHE BOOL "" FORCE)
set(TTK_ENABLE_ZFP OFF CACHE BOOL "" FORCE)
set(TTK_ENABLE_EIGEN OFF CACHE BOOL "" FORCE)
set(TTK_ENABLE_SCIKIT_LEARN OFF CACHE BOOL "" FORCE)
set(TTK_ENABLE_QHULL OFF CACHE BOOL "" FORCE)
set(TTK_CELL_ARRAY_LAYOUT SingleArray CACHE STRING "" FORCE)

add_subdirectory("${TTK_SOURCE_DIR}" ttk EXCLUDE_FROM_ALL)

add_executable(morseframes_ttk_gradient_benchmark
  benchmark_ttk_process_lower_stars.cpp
)
target_compile_features(morseframes_ttk_gradient_benchmark PRIVATE cxx_std_17)
target_include_directories(morseframes_ttk_gradient_benchmark PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}/../../include"
)
target_link_libraries(morseframes_ttk_gradient_benchmark PRIVATE
  discreteGradient
  Threads::Threads
)
