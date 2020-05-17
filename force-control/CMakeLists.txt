cmake_minimum_required(VERSION 3.1)
project(DDP-SQP-ADMM)

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

#find_package(MKL)
if(MKL_FOUND)
    message(STATUS "Intel MKL Found")
    message(STATUS "MKL root directory: ${MKL_ROOT_DIR}")
    message(STATUS "MKL include directory: ${MKL_INCLUDE_DIR}")
    message(STATUS "MKL libraries: ${MKL_LIBRARIES}")
    add_definitions(-DEIGEN_USE_MKL_ALL)
    include_directories(SYSTEM ${MKL_INCLUDE_DIR})
    link_libraries(${MKL_LIBRARIES})
endif()

find_package(OROCOS_KDL)
if (OROCOS_KDL)
endif()


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")



include_directories(SYSTEM ${EIGEN_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${OROCOS_KDL_INCLUDE_DIR})

#add_library( INTERFACE)
#target_include_directories(ddp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
#file(GLOB DDP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
#target_sources(ddp INTERFACE ${DDP_SOURCES})


add_library(force-control STATIC)
target_include_directories(force-control PUBLIC include)
target_sources(force-control PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/force_control/AdmittanceForceController.cpp)

add_library(kuka-models STATIC)
target_include_directories(kuka-models PUBLIC include)
set(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src_test/KukaModel.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src_test/models.cpp ${CMAKE_CURRENT_SOURCE_DIR}/analytical/RobotAnalytical.cpp ${CMAKE_CURRENT_SOURCE_DIR}/analytical/KUKAAnalyticalSolutions.cpp)
target_sources(kuka-models PUBLIC ${SOURCES})


#set(SOURCES force_control/AdmittanceForceController.cpp)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")


add_executable(KukaRobot src_test/KukaBot.cpp)

target_include_directories(KukaRobot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(KukaRobot orocos-kdl force-control kuka-models)



