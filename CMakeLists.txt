cmake_minimum_required(VERSION 3.16)

project(cytozoic VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use Vcpkg toolchain if available
if(DEFINED ENV{VCPKG_ROOT})
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
endif()

find_package(CGAL REQUIRED)

find_package(Boost 1.80 REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
endif()


find_package(Qt6 REQUIRED COMPONENTS Widgets)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

add_executable(cytozoic
    src/main.cpp
    src/main_window.cpp
    src/cytozoic_widget.cpp
    src/voronoi.cpp 
    src/geometry.cpp
    src/cytozoic.cpp
    src/neighborhood_indexer.cpp
    src/rules_dialog.cpp
    src/serialize.cpp
)

target_link_libraries(
    cytozoic PRIVATE 
    Qt6::Widgets  
    CGAL::CGAL
)

set_target_properties(cytozoic PROPERTIES
    WIN32_EXECUTABLE ON
    MACOSX_BUNDLE ON
)