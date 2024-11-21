cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ../..)

include(${ROOT_PATH}/extern/platform_core/cmake_utils/general.cmake)

project_guarded(gund_utility)

set(HEADERS
    ./include/gund_structs.h
    ./include/gund_utility.h
)

set(SOURCES
    ./src/gund_utility.cpp
)

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES})

generate_headers(${PROJECT_NAME} ${HEADERS})

target_include_directories(${PROJECT_NAME} PUBLIC ./include)

target_link_libraries(${PROJECT_NAME}
  PRIVATE
)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)