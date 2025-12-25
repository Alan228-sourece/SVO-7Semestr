cmake_minimum_required(VERSION 2.0)
project(inmost_solver_project)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_subdirectory(INMOST-master)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/INMOST-master/Source/Headers)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_executable(main main.cpp src/solver.cpp )
add_executable(visualize visualize.cpp src/solver.cpp )
target_link_libraries(visualize inmost)

target_link_libraries(main inmost)
set_target_properties(main PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
)
set_target_properties(visualize PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
)