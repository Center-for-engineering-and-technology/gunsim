#
# Copyright © CET MIPT LLC 2019
# CET-iDESC 2.0 (Intelligent Database and Expert System Constructor) (certificate №2019664842 FEDERAL INSTITUTE OF INDUSTRIAL PROPERTY)
#

### --- GCC
set(GCC_GENERAL_FLAGS "-Wall -Wno-long-long -pedantic")
set(GCC_CTEST_FLAGS "-Werror=return-type")
### --- Clang
set(CLANG_GENERAL_FLAGS "-Wall")
set(CLANG_CTEST_FLAGS "-Wextra -pedantic")
set(CLANG_ANALYZER_FLAGS "-Wno-unused-command-line-argument --analyze")
### --- Sanitizer
set(ADDRESS_SANITIZER_FLAGS "-fsanitize=address")
set(THREAD_SANITIZER_FLAGS "-fsanitize=thread")
set(MEMORY_SANITIZER_FLAGS "-fsanitize=memory")
set(UNDEFINED_SANITIZER_FLAGS "-fsanitize=undefined")
set(LEAK_SANITIZER_FLAGS "-fsanitize=leak")


### --- GCC
macro(gcc_general_flags)
  string(FIND "${CMAKE_CXX_FLAGS}" "${GCC_GENERAL_FLAGS}" res_find_flags)
  if (${res_find_flags} EQUAL -1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_GENERAL_FLAGS}")
  endif ()
endmacro()

macro(gcc_ctest_flags)
  string(FIND "${CMAKE_CXX_FLAGS}" "${GCC_CTEST_FLAGS}" res_find_flags)
  if (${res_find_flags} EQUAL -1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_CTEST_FLAGS}")
  endif ()
endmacro()
### --------------------------------------------------------------------


### --- Clang
macro(clang_general_flags)
  string(FIND "${CMAKE_CXX_FLAGS}" "${CLANG_GENERAL_FLAGS}" res_find_flags)
  if (${res_find_flags} EQUAL -1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CLANG_GENERAL_FLAGS}")
  endif ()
endmacro()

macro(clang_ctest_flags)
  string(FIND "${CMAKE_CXX_FLAGS}" "${CLANG_CTEST_FLAGS}" res_find_flags)
  if (${res_find_flags} EQUAL -1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CLANG_CTEST_FLAGS}")
  endif ()

  string(FIND "${CMAKE_CXX_FLAGS}" "${CLANG_ANALYZER_FLAGS}" res_find_flags)
  if (CLANG_ANALYZER AND ${res_find_flags} EQUAL -1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CLANG_ANALYZER_FLAGS}")
  endif ()
endmacro()
### --------------------------------------------------------------------


### --- Sanitizer
macro(address_sanitizer_flags)
  string(FIND "${CMAKE_CXX_FLAGS}" "${ADDRESS_SANITIZER_FLAGS}" res_find_flags)
  if (${res_find_flags} EQUAL -1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDRESS_SANITIZER_FLAGS}")
  endif ()
endmacro()

macro(thread_sanitizer_flags)
  string(FIND "${CMAKE_CXX_FLAGS}" "${THREAD_SANITIZER_FLAGS}" res_find_flags)
  if (${res_find_flags} EQUAL -1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${THREAD_SANITIZER_FLAGS}")
  endif ()
endmacro()

macro(memory_sanitizer_flags)
  string(FIND "${CMAKE_CXX_FLAGS}" "${MEMORY_SANITIZER_FLAGS}" res_find_flags)
  if (${res_find_flags} EQUAL -1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MEMORY_SANITIZER_FLAGS}")
  endif ()
endmacro()

macro(undefined_sanitizer_flags)
  string(FIND "${CMAKE_CXX_FLAGS}" "${UNDEFINED_SANITIZER_FLAGS}" res_find_flags)
  if (${res_find_flags} EQUAL -1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UNDEFINED_SANITIZER_FLAGS}")
  endif ()
endmacro()

macro(leak_sanitizer_flags)
  string(FIND "${CMAKE_CXX_FLAGS}" "${LEAK_SANITIZER_FLAGS}" res_find_flags)
  if (${res_find_flags} EQUAL -1)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LEAK_SANITIZER_FLAGS}")
  endif ()
endmacro()
### --------------------------------------------------------------------
