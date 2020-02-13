unset(SOURCES CACHE)

set(CURRENT_SRCS vp.cxx filter.h uf.h cpu.cxx cpu.h config.h)
add_to_sources("${CURRENT_SRCS}")

add_subdirectory(cvtcolor)

if(MSVC)
 configure_file(vp.h.in ${CMAKE_BINARY_DIR}/generated/vp.h)
 add_to_sources(${CMAKE_BINARY_DIR}/generated/vp.h)
 group("${SOURCES}")
endif()

add_library(vp SHARED ${SOURCES})

set_target_properties(vp PROPERTIES CXX_STANDARD 17 VISIBILITY_INLINES_HIDDEN 1 CXX_VISIBILITY_PRESET hidden C_VISIBILITY_PRESET hidden)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/generated)

if(VAPOURSYNTH_INCLUDE_DIR)
 include_directories(${VAPOURSYNTH_INCLUDE_DIR})
endif()

if(XXHASH_INCLUDE_DIR)
 include_directories(${XXHASH_INCLUDE_DIR})
endif()
