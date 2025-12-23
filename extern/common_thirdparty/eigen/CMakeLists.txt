cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)
set(ROOT_PATH ..)
include(${ROOT_PATH}/../platform_core/cmake_utils/general.cmake)
project_guarded(eigen)

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.11.0")
    file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/eigen/Eigen/*)
else()
    file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/eigen/Eigen/*)
endif()

#directory_source_group(./eigen/Eigen/)

add_library(${PROJECT_NAME} INTERFACE)
add_custom_target(${PROJECT_NAME}_i SOURCES ${HEADERS})

target_include_directories(${PROJECT_NAME} INTERFACE ./eigen)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
