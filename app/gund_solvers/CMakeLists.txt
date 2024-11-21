cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ../..)

include(${ROOT_PATH}/extern/platform_core/cmake_utils/general.cmake)

project_guarded(gund_solvers)

option(WITH_TESTS "Include test project")

add_subproject(gund_utility)
add_subproject(spectrum_solver)
add_subproject(gund_parser)

if(WITH_TESTS)
  add_subproject(spectrum_solver_test)
endif(WITH_TESTS)
