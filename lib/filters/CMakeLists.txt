cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ../..)

include(${ROOT_PATH}/extern/platform_core/cmake_utils/general.cmake)

project_guarded(filters)

find_package(Boost REQUIRED)

set(HEADERS
    ./include/common.h
    ./include/filter.h
    ./include/filterscommon.h
    ./include/q_filter.h
    ./include/wienerfilter.h
    ./include/bandpassfilter.h
    ./include/filter.h)

set(SOURCES ./src/common.cpp ./src/filterscommon.cpp ./src/q_filter.cpp
            ./src/wienerfilter.cpp ./src/bandpassfilter.cpp ./src/filter.cpp)

add_library(${PROJECT_NAME} STATIC ${SOURCES} ${HEADERS})

#generate_headers(${PROJECT_NAME} ${HEADERS})

target_include_directories(${PROJECT_NAME} PUBLIC ./include)
target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIR})

add_subproject(common_structs)
add_subproject(simple_fft_wrap)

target_link_libraries(${PROJECT_NAME} PUBLIC common_structs simple_fft_wrap)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
