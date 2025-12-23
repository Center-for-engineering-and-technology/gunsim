cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ../..)

include(${ROOT_PATH}/extern/platform_core/cmake_utils/general.cmake)

project_guarded(comparative_analysis)

set(HEADERS ./include/comparative_analysis.h)

set(SOURCES ./src/comparative_analysis.cpp)

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES})

#generate_headers(${PROJECT_NAME} ${HEADERS})

target_include_directories(${PROJECT_NAME} PUBLIC ./include)

add_subproject(filters)
add_subproject(gund_parser)
add_subproject(spectrum_solver)

add_subproject(extern/common_thirdparty/nlohmann)

target_link_libraries(${PROJECT_NAME} PUBLIC spectrum_solver nlohmann
                                             gund_parser filters)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
