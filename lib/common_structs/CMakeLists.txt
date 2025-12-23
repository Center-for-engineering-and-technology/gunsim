cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ../..)

include(${ROOT_PATH}/extern/platform_core/cmake_utils/general.cmake)

project_guarded(common_structs)

# Header-only units library
find_package(Boost REQUIRED)

set(HEADERS ./include/gund_structs.h ./include/spectrum_solver_structs.h
            ./include/gund_utility.h ./include/signature_structs.h)

set(SOURCES ./src/gund_utility.cpp)

add_subproject(extern/common_thirdparty/range-v3)
add_subproject(extern/common_thirdparty/eigen)

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES})

#generate_headers(${PROJECT_NAME} ${HEADERS})

target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                         $<INSTALL_INTERFACE:include> ${Boost_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME} PUBLIC eigen range-v3)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
