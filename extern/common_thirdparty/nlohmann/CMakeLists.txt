cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)
set(ROOT_PATH ..)
include(${ROOT_PATH}/../platform_core/cmake_utils/general.cmake)
project_guarded(nlohmann)

set(HEADERS
    ./json/single_include/nlohmann/json.hpp
    ./json/single_include/nlohmann/json_fwd.hpp
)

add_library(${PROJECT_NAME} INTERFACE)
add_custom_target(${PROJECT_NAME}_i SOURCES ${HEADERS})

target_include_directories(${PROJECT_NAME} INTERFACE ./json/single_include)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
