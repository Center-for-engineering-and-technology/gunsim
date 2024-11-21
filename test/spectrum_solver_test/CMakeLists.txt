cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ../..)

include(${ROOT_PATH}/extern/platform_core/cmake_utils/general.cmake)

project_guarded(spectrum_solver_test)

set(TESTS_SOURCES test.cpp)

set(Boost_USE_STATIC_LIBS ON)
find_package (Boost COMPONENTS unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_executable (${PROJECT_NAME} ${TESTS_SOURCES})

add_subproject(gund_parser)
add_subproject(spectrum_solver)

target_link_libraries(${PROJECT_NAME}
  PRIVATE
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  gund_parser
  spectrum_solver
)

macro(testcase arg)
  add_test (NAME ${arg} COMMAND ${PROJECT_NAME} -t ${arg})
endmacro()

enable_testing ()

#target_include_directories(${PROJECT_NAME} PUBLIC ../extern/gdm_thirdparty/eigen)
target_include_directories(${PROJECT_NAME} PUBLIC ${ROOT_PATH}/extern/Simple-FFT)

testcase(solver_check)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)