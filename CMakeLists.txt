cmake_minimum_required(VERSION 3.13)
include(CheckIncludeFiles)

function(group FILES)
 foreach(FILE ${FILES})
  string(LENGTH ${CMAKE_CURRENT_SOURCE_DIR} LEN)
  string(LENGTH ${FILE} SUBLEN)
  if(NOT LEN GREATER SUBLEN)
   string(SUBSTRING ${FILE} 0 ${LEN} PREFIX)
  endif()
  if(PREFIX STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
   source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${FILE})
  else()
   source_group("" FILES ${FILE})
  endif()
 endforeach()
endfunction()

function(add_to_sources FILES)
 foreach(FILE ${FILES})
  get_filename_component(ABSOLUTE_PATH ${FILE} ABSOLUTE)
  set(sources "${sources}" ${ABSOLUTE_PATH} CACHE INTERNAL "")
 endforeach()
endfunction()

function(add_to_avx_sources FILES)
 foreach(FILE ${FILES})
  get_filename_component(ABSOLUTE_PATH ${FILE} ABSOLUTE)
  set(avx_sources "${avx_sources}" ${ABSOLUTE_PATH} CACHE INTERNAL "")
 endforeach()
endfunction()

function(set_compile_options FILES OPTIONS)
 foreach(FILE ${FILES})
  get_filename_component(ABSOLUTE_PATH ${FILE} ABSOLUTE)
   set_source_files_properties(${ABSOLUTE_PATH} PROPERTIES COMPILE_OPTIONS "${OPTIONS}")
 endforeach()
endfunction()

project(vp CXX C ASM_NASM)

find_package(Git QUIET)

if(GIT_EXECUTABLE AND EXISTS ${CMAKE_SOURCE_DIR}/.git)
 execute_process(
  COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE LATEST_TAG
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
 )
 execute_process(
  COMMAND ${GIT_EXECUTABLE} rev-list ${LATEST_TAG}.. --count --first-parent
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE TAG_DISTANCE
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
 )
 execute_process(
  COMMAND ${GIT_EXECUTABLE} update-index -q --refresh
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_QUIET
  ERROR_QUIET
 )
 execute_process(
  COMMAND ${GIT_EXECUTABLE} diff-index --name-only HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_DIFF_INDEX
  ERROR_QUIET
 )
endif()

if(DEFINED LATEST_TAG AND NOT LATEST_TAG STREQUAL "")
 set(PROJECT_VERSION "${LATEST_TAG}")
 if(TAG_DISTANCE)
  string(APPEND PROJECT_VERSION ".${TAG_DISTANCE}")
 endif()
 if(NOT GIT_DIFF_INDEX STREQUAL "")
  string(APPEND PROJECT_VERSION "-dirty")
 endif()
else()
 set(PROJECT_VERSION "unknown")
endif()

list(APPEND CMAKE_LIBRARY_PATH ${LIBRARY_PATH} $ENV{LIBRARY_PATH})
list(APPEND CMAKE_INCLUDE_PATH ${INCLUDE_PATH} $ENV{INCLUDE_PATH})

find_path(VAPOURSYNTH_INCLUDE_DIR NAMES vapoursynth.h VapourSynth.h
 PATH_SUFFIXES vapoursynth)

find_path(XXHASH_INCLUDE_DIR NAMES xxhash.h)
find_library(XXHASH_LIB NAMES xxhash libxxhash)

add_subdirectory(sources)

message(STATUS "vp version ${PROJECT_VERSION}")
