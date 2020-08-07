cmake_minimum_required(VERSION 3.5)
project(Matplotlib)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include/python2.7")
set(CMAKE_CXX_STANDARD 11)


# find_package(Python2 COMPONENTS Development NumPy)




# add_library(ModernRoboticsCpp SHARED src/modern_robotics.cpp include/modern_robotics.h)

# set(LIB_TEST_SOURCES
#   ${CMAKE_CURRENT_SOURCE_DIR}/src/lib_test.cpp
# )

add_executable(myproject  examples/modern.cpp)
target_link_libraries(myproject python2.7)


# install (TARGETS ModernRoboticsCpp
#         ARCHIVE DESTINATION lib
#         LIBRARY DESTINATION lib
#         RUNTIME DESTINATION bin)
