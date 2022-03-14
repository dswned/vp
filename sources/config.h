#pragma once
#if defined _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS
#pragma warning(disable: 26812) // prefer 'enum class' over 'enum'
#if defined __cplusplus && _MSVC_LANG < 201703L
#pragma message("min stdc++ is 17 there")
#endif
#endif

#if defined __cplusplus
#include <cstdint> // c++ definitions
#else
#include <stdint.h>
#endif

#if defined _MSC_VER
#define FORCE_INLINE __forceinline
#define NO_INLINE __declspec(noinline)
#else
#define FORCE_INLINE inline __attribute__((always_inline))
#define NO_INLINE __attribute__((noinline))
#endif

#define VP_ENABLE_INTRIN

#if defined VP_ENABLE_INTRIN
#if defined __GNUC__ || defined __clang__
#include <immintrin.h>
#if !defined __AVX__ || !defined __F16C__ || !defined __AES__
#include <pmmintrin.h>
#include <smmintrin.h>
#include <avxintrin.h>
#include <wmmintrin.h>
#include <f16cintrin.h>
#endif
#define A_F16C __attribute__((target("f16c")))
#define A_AES __attribute__((target("aes")))
#define A_SSE4 __attribute__((target("sse4")))
#define A_AVX __attribute__((target("avx")))
#else
#include <intrin.h>
#define A_F16C
#define A_AES
#define A_SSE4
#define A_AVX
#endif
#endif

#if defined __cplusplus
#define ALIGNAS(N) alignas(N)
#define VP_EXTERN_C extern "C"
#elif defined _MSC_VER
#define ALIGNAS(N) __declspec(align(N))
#else
#define ALIGNAS(N) __attribute__((aligned(N)))
#endif

#ifdef __cpp_lib_is_constant_evaluated
#define CONST_EVAL_CONSTEXPR constexpr
#else
#define CONST_EVAL_CONSTEXPR
#endif
#if defined __GNUC__ || defined __clang__
#define EXPECT(x, y) __builtin_expect(x, y)
#endif
#if defined __clang__
#define ASSUME(x) __builtin_assume(x)
#elif defined __GNUC__
#define ASSUME(x) do { if (!(x)) __builtin_unreachable(); } while(0)
#else
#define ASSUME(x) __assume(x)
#endif

#define XXH_STATIC_LINKING_ONLY
#define FFTW_DLL
#define BOOST_AUTO_LINK_SYSTEM
#if defined _HAS_CXX20 && _HAS_CXX20
#define BOOST_COMPUTE_USE_CPP11
#endif

//#define BUILD_EXAMPLE
#define BUILD_CVTCOLOR
#define BUILD_RAWS
#define BUILD_ADDGRAIN
#define BUILD_DFTTEST
#define BUILD_NNEDI
