cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ../..)

include(${ROOT_PATH}/extern/platform_core/cmake_utils/general.cmake)

project_guarded(gund_parser)

set(HEADERS
    ./include/gund_format_parser.h
    ./include/gund_json_parser.h
)

set(SOURCES
    ./src/gund_format_parser.cpp
    ./src/gund_json_parser.cpp
)

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC ./include)

add_subproject(gund_utility)
add_subproject(spectrum_solver)
add_subproject(extern/common_thirdparty/nlohmann)

target_link_libraries(${PROJECT_NAME}
  PUBLIC
  gund_utility
  spectrum_solver
  nlohmann
)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)