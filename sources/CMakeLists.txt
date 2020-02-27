unset(avx_sources CACHE)
unset(sources CACHE)

add_to_sources("vp.cxx;filter.h;uf.h;cpu.cxx;cpu.h;config.h")

add_subdirectory(addgrain)
add_subdirectory(cvtcolor)

if(MSVC)
 configure_file(vp.h.in ${CMAKE_BINARY_DIR}/generated/vp.h)
 add_to_sources(${CMAKE_BINARY_DIR}/generated/vp.h)
 group("${sources};${avx_sources}")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
 set_compile_options("${avx_sources}" "-mavx;-maes")
else()
 set_compile_options("${avx_sources}" "/arch:AVX")
endif()

add_library(vp SHARED ${sources} ${avx_sources})

set_target_properties(vp PROPERTIES CXX_STANDARD 17 VISIBILITY_INLINES_HIDDEN 1 CXX_VISIBILITY_PRESET hidden C_VISIBILITY_PRESET hidden)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/generated)

if(VAPOURSYNTH_INCLUDE_DIR)
 include_directories(${VAPOURSYNTH_INCLUDE_DIR})
endif()

if(XXHASH_INCLUDE_DIR)
 include_directories(${XXHASH_INCLUDE_DIR})
endif()
