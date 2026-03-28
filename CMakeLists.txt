cmake_minimum_required(VERSION 3.16)

project(cytozoic VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Boost 1.80 REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
endif()

find_package(Eigen3 3.3 REQUIRED)

find_package(Qt6 REQUIRED COMPONENTS Widgets)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

add_executable(cytozoic
    src/main.cpp
    src/main_window.cpp
    src/cytozoic_widget.cpp
    src/voronoi.cpp 
    src/util.cpp
)

target_link_libraries(
    cytozoic PRIVATE 
    Qt6::Widgets  
    Eigen3::Eigen
)

set_target_properties(cytozoic PROPERTIES
    WIN32_EXECUTABLE ON
    MACOSX_BUNDLE ON
)