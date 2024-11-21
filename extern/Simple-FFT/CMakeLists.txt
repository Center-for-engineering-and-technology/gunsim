cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ..)
include(${ROOT_PATH}/platform_core/cmake_utils/general.cmake)
project_guarded(simple_fft_wrap)

set(HEADERS
    Simple-FFT/include/simple_fft/fft.h
    Simple-FFT/include/simple_fft/fft.hpp
    Simple-FFT/include/simple_fft/fft_impl.hpp
    Simple-FFT/include/simple_fft/fft_settings.h
    Simple-FFT/include/simple_fft/check_fft.hpp
    Simple-FFT/include/simple_fft/copy_array.hpp
    Simple-FFT/include/simple_fft/error_handling.hpp
)

directory_source_group(./)

add_library(${PROJECT_NAME} INTERFACE)
add_custom_target(${PROJECT_NAME}_i SOURCES ${HEADERS})

target_include_directories(${PROJECT_NAME} INTERFACE ./Simple-FFT/include)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
target_compile_definitions(${PROJECT_NAME}
    INTERFACE
        __USE_SQUARE_BRACKETS_FOR_ELEMENT_ACCESS_OPERATOR
)
