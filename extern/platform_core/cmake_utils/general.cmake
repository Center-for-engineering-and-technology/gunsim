#
# Copyright © CET MIPT LLC 2019
# CET-iDESC 2.0 (Intelligent Database and Expert System Constructor) (certificate №2019664842 FEDERAL INSTITUTE OF INDUSTRIAL PROPERTY)
#

include(${CMAKE_CURRENT_LIST_DIR}/current_variables.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/compiler_flags.cmake)
include(ProcessorCount)
set(PLATFORM_CORE_ROOT_PATH ${CMAKE_CURRENT_LIST_DIR}/..)

ProcessorCount(N)

macro(project_guarded name)
  if(DEFINED ${name}_GUARD)
    if(NOT ${name}_GUARD STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
      return()
    endif()
  else()
    set(${name}_GUARD ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "${name} guard")
  endif()
  project(${name})

  # Warnings level
  cmake_policy(SET CMP0054 NEW)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")

    string(FIND "${CMAKE_CXX_FLAGS}" "/MP /EHsc /utf-8" res_find_flags)

    if (${res_find_flags} EQUAL -1)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /EHsc /utf-8")
    endif ()

    # Force to always compile with W4
    if (CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
      string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else ()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif ()
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    clang_general_flags()
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    gcc_general_flags()
  endif ()


  ### --- CDASH
  if (("${CDASH_PROJECT}" STRGREATER "") AND ("${CDASH_TOKEN}" STRGREATER ""))
    if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      clang_ctest_flags()
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      gcc_ctest_flags()
    endif()
  endif ()


  ### --- SANITIZE
  if (SANITIZE)
    if (SANITIZE STREQUAL "address")
      address_sanitizer_flags()
    elseif (SANITIZE STREQUAL "thread")
      thread_sanitizer_flags()
    elseif (SANITIZE STREQUAL "memory")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=memory -fno-omit-frame-pointer -fsanitize-memory-track-origins -fsanitize-blacklist=${ABSOLUTE_PATH_CMAKE_UTILS}/blacklist_memory_sanitizer.txt")
    elseif (SANITIZE STREQUAL "undefined")
      undefined_sanitizer_flags()
    elseif (SANITIZE STREQUAL "leak")
      leak_sanitizer_flags()
    elseif (SANITIZE STREQUAL "purify")
      set(CTEST_MEMORYCHECK_TYPE "Purify")
    elseif (SANITIZE STREQUAL "bounds")
      set(CTEST_MEMORYCHECK_TYPE "BoundsChecker")
    elseif (SANITIZE STREQUAL "drmemory")
      set(CTEST_MEMORYCHECK_TYPE "DrMemory")
    elseif (SANITIZE STREQUAL "cuda")
      set(CTEST_MEMORYCHECK_TYPE "CudaSanitizer")
    endif ()
  endif ()

  # Extensions
  set(CMAKE_CXX_EXTENSIONS OFF)

  # Optimizations
  if(MSVC)

	if (CRASH_HANDLE)
	  #remove(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/O2")
	  remove(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/Ob1")
	  remove(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG")
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /O2 /Ob2" )
      #set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /Od" )
	endif (CRASH_HANDLE)

    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 -DNDEBUG")
	
  elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG" )
  
  endif()

  # Warnings as errors
  if (MSVC)
    add_compile_options(/we4715) # 'function': not all control paths return a value
	add_compile_options(/we4172) # returning address of local variable or temporary
  endif()

  # DEBUG OUTPUT
  if (BUILD_DEBUG_OUTPUT)

    if (NOT DEFINED BUILD_SOURCE_DIRECTORIES_GUARD)

      file(WRITE ${CMAKE_BINARY_DIR}/build_source_directories.txt "")
      set(BUILD_SOURCE_DIRECTORIES_GUARD ${CMAKE_BINARY_DIR} CACHE INTERNAL "BUILD_SOURCE_DIRECTORIES_GUARD")

    endif()

	file(APPEND ${CMAKE_BINARY_DIR}/build_source_directories.txt "${CMAKE_CURRENT_SOURCE_DIR}\n")

  endif()

endmacro()

macro(add_subproject name)

  set(REPOSITORY_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
  set(REPOSITORY_ROOT_DIR_EXISTS FALSE)

  while(EXISTS ${REPOSITORY_ROOT_DIR})

    if(EXISTS "${REPOSITORY_ROOT_DIR}/cmake_project")

      set(REPOSITORY_ROOT_DIR_EXISTS TRUE)
      break()

    else()

      set(REPOSITORY_ROOT_DIR ${REPOSITORY_ROOT_DIR}/..)

    endif()

  endwhile()

  if (NOT ${REPOSITORY_ROOT_DIR_EXISTS})

    message(FATAL_ERROR cmake_project " DOES NOT EXIST FOR " ${CMAKE_CURRENT_SOURCE_DIR})
	return()

  else()
  endif()

  # DEBUG OUTPUT
  if (BUILD_DEBUG_OUTPUT)

    if (NOT DEFINED BUILD_SOURCE_MODULES_GUARD)

      file(WRITE ${CMAKE_BINARY_DIR}/build_source_modules.txt "")
	  
	  get_filename_component(ABS_REPOSITORY_ROOT_DIR ${REPOSITORY_ROOT_DIR} ABSOLUTE)
	  file(WRITE ${CMAKE_BINARY_DIR}/build_source_root.txt "${ABS_REPOSITORY_ROOT_DIR}")

      set(BUILD_SOURCE_MODULES_GUARD ${CMAKE_BINARY_DIR} CACHE INTERNAL "BUILD_SOURCE_MODULES_GUARD")

    endif()

  endif()
  
  set(SUBPROJECT_PATH "")

  # Наличие extern
  string(REPLACE "/" ";" SUBPROJECT_DIRECTORY_LIST ${name})
  list(GET SUBPROJECT_DIRECTORY_LIST 0 SUBPROJECT_NAME)
  
  if(${SUBPROJECT_NAME} STREQUAL "extern")

    list(GET SUBPROJECT_DIRECTORY_LIST 1 EXTERN_REPOSITORY_NAME)
    string(REPLACE "extern/${EXTERN_REPOSITORY_NAME}/" "" SUBPROJECT_INFO_FILE_PATH ${name})

    if(EXISTS "${REPOSITORY_ROOT_DIR}/extern")

      if(EXISTS "${REPOSITORY_ROOT_DIR}/extern/${EXTERN_REPOSITORY_NAME}/cmake_project/${SUBPROJECT_INFO_FILE_PATH}.cmake")

        set(SUBPROJECT_PATH "${REPOSITORY_ROOT_DIR}/extern/${EXTERN_REPOSITORY_NAME}/cmake_project/${SUBPROJECT_INFO_FILE_PATH}.cmake")

	    # Read subproject path
        file(STRINGS "${REPOSITORY_ROOT_DIR}/extern/${EXTERN_REPOSITORY_NAME}/cmake_project/${SUBPROJECT_INFO_FILE_PATH}.cmake" SUBPROJECT_DIRECTORY)

	    # Get subproject name
        string(REPLACE "/" ";" SUBPROJECT_DIRECTORY_LIST ${name})
        list(REVERSE SUBPROJECT_DIRECTORY_LIST)
        list(GET SUBPROJECT_DIRECTORY_LIST 0 SUBPROJECT_NAME)

        # Add subproject
        add_subdirectory(${REPOSITORY_ROOT_DIR}/extern/${EXTERN_REPOSITORY_NAME}/${SUBPROJECT_DIRECTORY} ${CMAKE_BINARY_DIR}/subproject/${SUBPROJECT_NAME}/${PROJECT_NAME})

      else()

        message(FATAL_ERROR ${REPOSITORY_ROOT_DIR}/extern/${EXTERN_REPOSITORY_NAME}/cmake_project/${SUBPROJECT_INFO_FILE_PATH}.cmake " DOES NOT EXIST")

      endif()

    else()

	  if(EXISTS "${REPOSITORY_ROOT_DIR}/../${EXTERN_REPOSITORY_NAME}/cmake_project/${SUBPROJECT_INFO_FILE_PATH}.cmake")

        set(SUBPROJECT_PATH "${REPOSITORY_ROOT_DIR}/../${EXTERN_REPOSITORY_NAME}/cmake_project/${SUBPROJECT_INFO_FILE_PATH}.cmake")

	    # Read subproject path
        file(STRINGS "${REPOSITORY_ROOT_DIR}/../${EXTERN_REPOSITORY_NAME}/cmake_project/${SUBPROJECT_INFO_FILE_PATH}.cmake" SUBPROJECT_DIRECTORY)

	    # Get subproject name
        string(REPLACE "/" ";" SUBPROJECT_DIRECTORY_LIST ${name})
        list(REVERSE SUBPROJECT_DIRECTORY_LIST)
        list(GET SUBPROJECT_DIRECTORY_LIST 0 SUBPROJECT_NAME)

        # Add subproject
        add_subdirectory(${REPOSITORY_ROOT_DIR}/../${EXTERN_REPOSITORY_NAME}/${SUBPROJECT_DIRECTORY} ${CMAKE_BINARY_DIR}/subproject/${SUBPROJECT_NAME}/${PROJECT_NAME})

	  else()

        message(FATAL_ERROR ${REPOSITORY_ROOT_DIR}/../${EXTERN_REPOSITORY_NAME}/cmake_project/${SUBPROJECT_INFO_FILE_PATH}.cmake " DOES NOT EXIST")

      endif()

    endif()

  else()

    if(EXISTS "${REPOSITORY_ROOT_DIR}/cmake_project/${name}.cmake")

      set(SUBPROJECT_PATH "${REPOSITORY_ROOT_DIR}/cmake_project/${name}.cmake")

      # Read subproject path
      file(STRINGS "${REPOSITORY_ROOT_DIR}/cmake_project/${name}.cmake" SUBPROJECT_DIRECTORY)


      # Get subproject name
      string(REPLACE "/" ";" SUBPROJECT_DIRECTORY_LIST ${name})
      list(REVERSE SUBPROJECT_DIRECTORY_LIST)
      list(GET SUBPROJECT_DIRECTORY_LIST 0 SUBPROJECT_NAME)

      # Add subproject
      add_subdirectory(${REPOSITORY_ROOT_DIR}/${SUBPROJECT_DIRECTORY} ${CMAKE_BINARY_DIR}/subproject/${SUBPROJECT_NAME}/${PROJECT_NAME})

    else()

      message(FATAL_ERROR ${REPOSITORY_ROOT_DIR}/cmake_project/${name}.cmake " DOES NOT EXIST")

    endif()

  endif()

  # DEBUG OUTPUT
  if (BUILD_DEBUG_OUTPUT)

    if (NOT DEFINED ${SUBPROJECT_NAME}_DEBUG_GUARD)

      if (EXISTS "${SUBPROJECT_PATH}")

        get_filename_component(ABS_SUBPROJECT_PATH ${SUBPROJECT_PATH} ABSOLUTE)

        file(APPEND ${CMAKE_BINARY_DIR}/build_source_modules.txt "${ABS_SUBPROJECT_PATH}\n")

      endif()

      set(${SUBPROJECT_NAME}_DEBUG_GUARD ${CMAKE_BINARY_DIR} CACHE INTERNAL "${SUBPROJECT_NAME}_DEBUG_GUARD")

    endif()

  endif()
  
endmacro()

# taken from
#https://stackoverflow.com/questions/28344564/cmake-remove-a-compile-flag-for-a-single-translation-unit
#########################
#
# Applies CMAKE_CXX_FLAGS to all targets in the current CMake directory.
# After this operation, CMAKE_CXX_FLAGS is cleared.
#
macro(apply_global_cxx_flags_to_all_targets)
  #составляем список исключений ctestexcl для CTest
  #иначе будет ошибка called with non-compilable target type
  list(APPEND left Continuous Experimental Nightly)
  list(APPEND right Build Configure Coverage MemCheck Start Submit Test Update)
  list(APPEND right MemoryCheck)
  list(APPEND ctestexcl)
  foreach(_left ${left})
    list(APPEND ctestexcl ${_left})
    foreach(_right ${right})
      list(APPEND ctestexcl ${_left}${_right})
    endforeach()
  endforeach()
  #
  separate_arguments(_global_cxx_flags_list UNIX_COMMAND ${CMAKE_CXX_FLAGS})
  get_property(_targets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
  foreach(_target ${_targets})
    list(FIND ctestexcl ${_target} rslt)
    #if((${_target} STREQUAL "Experimental") OR (${_target} STREQUAL "Nightly") OR (${_target} STREQUAL "Continuous") OR (${_target} STREQUAL "NightlyMemoryCheck"))
    if(NOT ${rslt} EQUAL -1)
      continue()
    else ()
      target_compile_options(${_target} PUBLIC ${_global_cxx_flags_list})
    endif()
  endforeach()
  unset(CMAKE_CXX_FLAGS)
  set(_flag_sync_required TRUE)
endmacro()

#
# Removes the specified compile flag from the specified target.
#   _target     - The target to remove the compile flag from
#   _flag       - The compile flag to remove
#
# Pre: apply_global_cxx_flags_to_all_targets() must be invoked.
#
macro(remove_flag_from_target _target _flag)
    get_target_property(_target_cxx_flags ${_target} COMPILE_OPTIONS)
    if(_target_cxx_flags)
        list(REMOVE_ITEM _target_cxx_flags ${_flag})
        set_target_properties(${_target} PROPERTIES COMPILE_OPTIONS "${_target_cxx_flags}")
    endif()
endmacro()
#########################
