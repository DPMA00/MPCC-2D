cmake_minimum_required(VERSION 3.20)

project(mpcc_core LANGUAGES CXX)

find_package(Eigen3 REQUIRED)
find_package(qpOASES)

# Add mpcc as a library
add_library(mpcc_core_lib
    src/Spline1D.cpp
    src/ParametricSpline.cpp
    src/Integrator.cpp
    src/ForwardEuler.cpp
    src/Model.cpp
    src/DiffDriveModel.cpp
    src/SecondOrderModel.cpp
    src/RunningCost.cpp
    src/DiffDriveCost.cpp
    src/SecondOrderModelCost.cpp
    src/ExtendedBicycleModel.cpp
    src/ExtendedBicycleCost.cpp
    src/RunningConstraintModel.cpp
    src/NoConstraints.cpp
    src/SecondOrderModelConstraints.cpp
    src/MPCC.cpp)

target_include_directories(mpcc_core_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Link the libraries with Eigen    
target_link_libraries(mpcc_core_lib 
    Eigen3::Eigen
    qpOASES)

option(BUILD_TESTS "Build all test executables" ON)
if(BUILD_TESTS)
    add_executable(Test1 test/test1.cpp)
    target_link_libraries(Test1 mpcc_core_lib)

    add_executable(Test2 test/test2.cpp)
    target_link_libraries(Test2 mpcc_core_lib)

    #add_executable(Test3 test/test3.cpp)
    #target_link_libraries(Test3 mpcc_core_lib)

endif()


set(CMAKE_EXPORT_COMPILE_COMMANDS ON)