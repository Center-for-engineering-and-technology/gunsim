cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ..)
include(${ROOT_PATH}/../platform_core/cmake_utils/general.cmake)
project_guarded(range-v3)

#directory_source_group(./)

#lib option handle
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(RANGE_V3_DOCS OFF)
set(RANGE_V3_TESTS OFF)
set(RANGE_V3_EXAMPLES OFF)
set(RANGE_V3_PERF OFF)
set(RANGE_V3_HEADER_CHECKS OFF)

add_subdirectory(range-v3)

target_include_directories(${PROJECT_NAME} INTERFACE ./range-v3/include)

target_link_libraries(${PROJECT_NAME}
    INTERFACE
        range-v3::range-v3
)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
