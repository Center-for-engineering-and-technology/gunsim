cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ../..)

include(${ROOT_PATH}/extern/platform_core/cmake_utils/general.cmake)

project_guarded(gun_model)

find_package(OpenMP REQUIRED)

set(HEADERS
    ./include/gun_model.h
    ./include/sig_optimizer.h
)

set(SOURCES
    ./src/gun_model.cpp
    ./src/sig_optimizer.cpp
)

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES})

target_include_directories(${PROJECT_NAME} 
  PUBLIC 
  ./include 
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(${PROJECT_NAME}
  PUBLIC
)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
