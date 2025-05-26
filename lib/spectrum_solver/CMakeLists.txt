cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ../..)

include(${ROOT_PATH}/extern/platform_core/cmake_utils/general.cmake)

project_guarded(spectrum_solver)

set(Boost_USE_STATIC_LIBS ON)
find_package (Boost REQUIRED)

set(HEADERS
    ./include/spectrum_solver.h
    ./include/spectrum_solver_structs.h
    ./include/fourier_solver.h
    ./include/internal_filter.h
    ./include/direct_diag.h
    ./include/energy_comp.h
)

set(SOURCES
    ./src/spectrum_solver.cpp
    ./src/spectrum_solver_structs.cpp
    ./src/fourier_solver.cpp
    ./src/internal_filter.cpp
    ./src/direct_diag.cpp
    ./src/energy_comp.cpp
)

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC ./include)
target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIR})

#Simple-FFT
#add_subdirectory(${REPOSITORY_ROOT_DIR}/extern/Simple-FFT ${CMAKE_BINARY_DIR}/subproject/Simple-FFT)
# target_include_directories(${PROJECT_NAME} PUBLIC ${ROOT_PATH}/extern/Simple-FFT)

add_subproject(gund_utility)
add_subproject(gund_parser)

add_subproject(simple_fft_wrap)

target_link_libraries(${PROJECT_NAME}
  PUBLIC
  #Simple-FFT
  gund_utility
  gund_parser
  simple_fft_wrap
)

#Runtime directory
  if (MSVC)

    if (NOT DEFINED APP_BUILD_DIRECTORY_DEBUG)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/solver/Debug/gund_solvers)
    else ()
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${APP_BUILD_DIRECTORY_DEBUG}/solver/gund_solvers)
    endif()

    if (NOT DEFINED APP_BUILD_DIRECTORY_RELEASE)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/solver/Release/gund_solvers)
    else ()
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${APP_BUILD_DIRECTORY_RELEASE}/solver/gund_solvers)
    endif()

    if (NOT DEFINED APP_BUILD_DIRECTORY_RELWITHDEBINFO)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_BINARY_DIR}/solver/RelWithDebInfo/gund_solvers)
    else ()
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${APP_BUILD_DIRECTORY_RELWITHDEBINFO}/solver/gund_solvers)
    endif()

  else (UNIX)
  
    if (NOT DEFINED APP_BUILD_DIRECTORY)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/solver/gund_solvers)
    else ()
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${APP_BUILD_DIRECTORY}/solver/gund_solvers)
    endif()

  endif (MSVC)

if (MSVC)

    FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../../data/ DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG}/data)
    FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../../data/ DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE}/data)
    FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../../data/ DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO}/data)

else (UNIX)

    FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../../data/ DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/data)

endif(MSVC)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
