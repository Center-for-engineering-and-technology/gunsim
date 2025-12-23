cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ../..)

include(${ROOT_PATH}/extern/platform_core/cmake_utils/general.cmake)

project_guarded(internal_filter)

set(HEADERS
    ./include/internal_filter.h
)

set(SOURCES
    ./src/internal_filter.cpp
)

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES})

#generate_headers(${PROJECT_NAME} ${HEADERS})

target_include_directories(${PROJECT_NAME} PUBLIC ./include)

add_subproject(common_structs)

target_link_libraries(${PROJECT_NAME}
  PUBLIC
  common_structs
)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
