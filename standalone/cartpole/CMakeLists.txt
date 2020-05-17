find_package(Qt5Widgets)
find_package(Qt5PrintSupport)

if(Qt5Widgets_FOUND AND Qt5PrintSupport_FOUND)
    message(STATUS "Found Qt5::Widgets and Qt5::PrintSupport")
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")
    add_definitions(-DUSE_GUI ${Qt5Widgets_DEFINITIONS})
    qt5_wrap_ui(CARTPOLE_UI_HDRS cartpolewindow.ui)
    add_executable(CartPole cartpole.cpp
            ${CMAKE_SOURCE_DIR}/examples/qcustomplot/qcustomplot.cpp
            cartpolewindow.h
            cartpolewindow.cpp
            ${CARTPOLE_UI_HDRS})
    target_include_directories(CartPole PUBLIC ${Qt5Widgets_INCLUDE_DIRS}
            ${Qt5PrintSupport_INCLUDE_DIRS}
            ${CMAKE_SOURCE_DIR}/examples)
    target_link_libraries(CartPole Qt5::Widgets Qt5::PrintSupport)
else()
    message(WARNING "Qt5 not found -- CartPole simulation not built")
endif()
