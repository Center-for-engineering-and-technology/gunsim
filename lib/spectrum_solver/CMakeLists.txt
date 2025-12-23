cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ../..)

include(${ROOT_PATH}/extern/platform_core/cmake_utils/general.cmake)

project_guarded(spectrum_solver)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED)

set(HEADERS ./include/spectrum_solver.h ./include/fourier_solver.h
            ./include/direct_diag.h ./include/energy_comp.h)

set(SOURCES ./src/spectrum_solver.cpp ./src/fourier_solver.cpp
            ./src/direct_diag.cpp ./src/energy_comp.cpp)

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES})

#generate_headers(${PROJECT_NAME} ${HEADERS})

target_include_directories(${PROJECT_NAME} PUBLIC ./include)
target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIR})

add_subproject(common_structs)
add_subproject(gund_parser)
add_subproject(gun_model)
add_subproject(simple_fft_wrap)
add_subproject(filters)

target_link_libraries(${PROJECT_NAME} PUBLIC common_structs gund_parser
                                             gun_model simple_fft_wrap filters)

# Runtime directory
if(MSVC)

  if(NOT DEFINED APP_BUILD_DIRECTORY_DEBUG)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG
        ${CMAKE_BINARY_DIR}/solver/Debug/gund_solvers)
  else()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG
        ${APP_BUILD_DIRECTORY_DEBUG}/solver/gund_solvers)
  endif()

  if(NOT DEFINED APP_BUILD_DIRECTORY_RELEASE)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE
        ${CMAKE_BINARY_DIR}/solver/Release/gund_solvers)
  else()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE
        ${APP_BUILD_DIRECTORY_RELEASE}/solver/gund_solvers)
  endif()

  if(NOT DEFINED APP_BUILD_DIRECTORY_RELWITHDEBINFO)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO
        ${CMAKE_BINARY_DIR}/solver/RelWithDebInfo/gund_solvers)
  else()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO
        ${APP_BUILD_DIRECTORY_RELWITHDEBINFO}/solver/gund_solvers)
  endif()

else(UNIX)

  if(NOT DEFINED APP_BUILD_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/solver/gund_solvers)
  else()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
        ${APP_BUILD_DIRECTORY}/solver/gund_solvers)
  endif()

endif(MSVC)

if(MSVC)

  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../../data/
       DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG}/data)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../../data/
       DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE}/data)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../../data/
       DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO}/data)

else(UNIX)

  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../../data/
       DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/data)

endif(MSVC)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
