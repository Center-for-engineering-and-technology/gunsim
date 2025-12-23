cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ../..)

include(${ROOT_PATH}/extern/platform_core/cmake_utils/general.cmake)

project_guarded(gun_model)

find_package(OpenMP REQUIRED)

set(HEADERS
    ./include/gun_model.h ./include/sig_optimizer.h
    ./include/gun_parameters_table.h ./include/notional_source_processor.h)

set(SOURCES ./src/gun_model.cpp ./src/sig_optimizer.cpp
            ./src/gun_parameters_table.cpp ./src/notional_source_processor.cpp)

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES})

#generate_headers(${PROJECT_NAME} ${HEADERS})

target_include_directories(${PROJECT_NAME} PUBLIC ./include)

add_subproject(common_structs)
add_subproject(filters)
add_subproject(extern/common_thirdparty/nlohmann)
add_subproject(gund_parser)

# Проверяем наличие библиотеки
if(TARGET OpenMP::OpenMP_CXX)
  message(STATUS "OpenMP found.")
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
else()
  message(FATAL_ERROR "[FATAL_ERROR]: OpenMP not found!!!!")
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC common_structs nlohmann
                                             gund_parser filters)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
