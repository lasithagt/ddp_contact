cmake_minimum_required(VERSION 3.1)
project(force-control)
include(GenerateExportHeader)

set (CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/install)

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/install/cmake)

option(DDP_SYSTEM_EIGEN "Use system-installed version of Eigen" OFF)
option(DDP_BUILD_EXAMPLES "Build examples" ON)

set(EIGEN_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/eigen)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# find system eigen
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


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


# include_directories(SYSTEM ${EIGEN_INCLUDE_DIR})
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/contact_model ${OROCOS_KDL_INCLUDE_DIR})


add_library(force-control INTERFACE)
# target_include_directories(force-control INTERFACE include/AdmittanceForceController.h include/SoftContactModel.h)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")


#add_library( INTERFACE)
#target_include_directories(ddp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
#file(GLOB DDP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
#target_sources(ddp INTERFACE ${DDP_SOURCES})


# set(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src_test/KukaModel.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src_test/models.cpp ${CMAKE_CURRENT_SOURCE_DIR}/analytical/RobotAnalytical.cpp ${CMAKE_CURRENT_SOURCE_DIR}/analytical/KUKAAnalyticalSolutions.cpp)
# add_library(kuka-model-analytical STATIC ${SOURCES})
# target_include_directories(kuka-model-analytical PUBLIC include)


# add_executable(KukaRobot src_test/KukaBot.cpp)
# target_include_directories(KukaRobot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_link_libraries(KukaRobot orocos-kdl force-control kuka-model-analytical)




###########################################INSTALL AND EXPORT##############################################
include(CMakePackageConfigHelpers)
install (TARGETS force-control
        EXPORT ForceControlTargets
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
)
# add include directories
target_include_directories(force-control
                           INTERFACE
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                           $<INSTALL_INTERFACE:include>
)

# install header file
install(FILES include/AdmittanceForceController.h include/SoftContactModel.h  DESTINATION include)

# # set version
set(version 3.4.1)

set_property(TARGET force-control PROPERTY VERSION ${version})
set_property(TARGET force-control PROPERTY SOVERSION 3)
set_property(TARGET force-control PROPERTY
  INTERFACE_force-control_MAJOR_VERSION 3)
set_property(TARGET force-control APPEND PROPERTY
  COMPATIBLE_INTERFACE_STRING force-control_MAJOR_VERSION
)

# generate the export targets for the build tree
export(EXPORT ForceControlTargets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/ForceControlTargets.cmake"
)

# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/ForceControlConfigVersion.cmake"
  VERSION "${version}"
  COMPATIBILITY AnyNewerVersion
)

# create config file
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/ForceControlConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/ForceControlConfig.cmake"
  INSTALL_DESTINATION cmake
)

# install config files
install(FILES
          "${CMAKE_CURRENT_BINARY_DIR}/ForceControlConfig.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/ForceControlConfigVersion.cmake"
        DESTINATION cmake
)


# generate and install export file
install(EXPORT ForceControlTargets
        FILE ForceControlTargets.cmake
        DESTINATION cmake
)

