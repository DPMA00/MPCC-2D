cmake_minimum_required(VERSION 3.20)

project(mpcc_core LANGUAGES CXX)

find_package(Eigen3 REQUIRED)
find_package(qpOASES)

# Add mpcc as a library
add_library(mpcc_core_lib
    src/Spline1D.cpp
    src/ParametricSpline.cpp
    src/MPCC.cpp)

target_include_directories(mpcc_core_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Link the libraries with Eigen    
target_link_libraries(mpcc_core_lib 
    Eigen3::Eigen
    qpOASES)

option(BUILD_SPLINE_TEST "Build spline test executable" ON)
if(BUILD_SPLINE_TEST)
    add_executable(SplineTest test/main.cpp)
    target_link_libraries(SplineTest mpcc_core_lib)
endif()


set(CMAKE_EXPORT_COMPILE_COMMANDS ON)