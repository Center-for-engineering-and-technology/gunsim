cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)

set(ROOT_PATH ../..)

include(${ROOT_PATH}/extern/platform_core/cmake_utils/general.cmake)

project_guarded(gund_solvers)

option(WITH_TESTS "Include test project")

add_subproject(common_structs)
add_subproject(spectrum_solver)
add_subproject(gund_parser)
add_subproject(internal_filter)
add_subproject(gun_model)
add_subproject(comparative_analysis)


if(WITH_TESTS)
  #add_subproject(spectrum_solver_test)
endif(WITH_TESTS)
