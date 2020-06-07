cmake_minimum_required(VERSION 3.1)
project(DDP-SQP-ADMM)
SET (CMAKE_CXX_COMPILER "/usr/bin/g++")


if(NOT CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(DDP_SYSTEM_EIGEN "Use system-installed version of Eigen" OFF)
option(DDP_BUILD_EXAMPLES "Build examples" ON)

set(EIGEN_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/eigen)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(DDP_SYSTEM_EIGEN)
    find_package(Eigen)
    if(EIGEN_FOUND)
        message(STATUS "Eigen found: ${EIGEN_INCLUDE_DIR}")
    else()
        message(FATAL_ERROR "Eigen NOT found!")
    endif()
endif()


find_package(OROCOS_KDL)
if (OROCOS_KDL)
endif()

# set(OROCOS_KDL_INCLUDE_DIR ${OROCOS_KDL_INCLUDE_DIR})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-march=native -mtune=native -mavx2 -mfma")
set(CMAKE_CXX_FLAGS "-std=c++11")


include_directories(SYSTEM ${EIGEN_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${OROCOS_KDL_INCLUDE_DIR} /usr/local/include)

message(STATUS ${OROCOS_KDL_INCLUDE_DIR})

add_library(kuka-models STATIC)
target_include_directories(kuka-models PUBLIC include)
set(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/KukaModel.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/models.cc)
target_sources(kuka-models PUBLIC ${SOURCES})

add_library(ddp-solver STATIC)
target_include_directories(ddp-solver PUBLIC include)
set(SOURCES_DDP ${CMAKE_CURRENT_SOURCE_DIR}/src/ilqrsolver.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/cost_function_kuka_arm.cc ${CMAKE_CURRENT_SOURCE_DIR}/src/kuka_arm.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/ddp.cpp)
target_link_libraries(ddp-solver PUBLIC kuka-models)
target_sources(ddp-solver PUBLIC ${SOURCES_DDP})

#add_library(ddp-admm-solver STATIC)
#target_include_directories(ddp-admm-solver PUBLIC include)
#set(SOURCES_DDP_ADMM ${CMAKE_CURRENT_SOURCE_DIR}/src/ilqrsolver_track.cc ${CMAKE_CURRENT_SOURCE_DIR}/src/cost_function_kuka_arm_track.cc #${CMAKE_CURRENT_SOURCE_DIR}/src/kuka_arm_track.cc)
#target_link_libraries(ddp-admm-solver PUBLIC kuka-models)
#target_sources(ddp-admm-solver PUBLIC ${SOURCES_DDP_ADMM})

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")


#######################################################################################


#add_executable(ddp-contact src/runddp_softcontact.cc)
#target_include_directories(ddp-contact PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
#target_link_libraries(ddp-contact orocos-kdl ddp-solver kuka-models)


#add_executable(ddp-admm-contact src/runADMM_softcontact.cc)
#target_include_directories(ddp-admm-contact PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
#target_link_libraries(ddp-admm-contact orocos-kdl ddp-admm-solver kuka-models)


add_executable(ddp-optim src/ddp_main.cpp)
target_include_directories(ddp-optim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ddp-optim orocos-kdl ddp-solver kuka-models)


add_executable(mpc-optim src/mpc_main.cpp)
target_include_directories(mpc-optim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(mpc-optim orocos-kdl ddp-solver kuka-models)
