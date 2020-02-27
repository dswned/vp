#pragma once
#define XXH_STATIC_LINKING_ONLY
#define FFTW_DLL
#define BOOST_AUTO_LINK_SYSTEM

#define BUILD_ADDGRAIN
#define BUILD_CVTCOLOR

#if defined(_WIN32)
#define strtok_r strtok_s
#endif

#if defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#define NO_INLINE __declspec(noinline)
#else
#define FORCE_INLINE inline __attribute__((always_inline))
#define NO_INLINE __attribute__((noinline))
#endif

#if defined(__cplusplus)
#define ALIGNAS(N) alignas(N)
#elif defined(_MSC_VER)
#define ALIGNAS(N) __declspec(align(N))
#else
#define ALIGNAS(N) __attribute__((aligned(N)))
#endif