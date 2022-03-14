#pragma once
#include "config.h"

#include <type_traits>
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <numeric>
#include <utility>
#include <algorithm>
#include <iterator>
#include <mutex>
#include <condition_variable>

#include <cstring>
#include <cstdarg>
#include <cstdio>

#include <xxhash.h>

namespace uf {

constexpr int bsr_const(uint64_t n)
{
	if (!n)
		return -1;
	uint64_t a = n, i = 0;
	int b = 0, j = 64;
	do
	{
		j >>= 1;
		i = 1ull << j;
		if (a >= i)
		{
			a >>= j;
			b += j;
		}
	} while (j);
	return b;
}

CONST_EVAL_CONSTEXPR inline int bsr(uint64_t n)
{
#ifdef __cpp_lib_is_constant_evaluated
	if (std::is_constant_evaluated())
		return bsr_const(n);
#endif
#if defined(_MSC_VER) && !defined(__clang__)
	unsigned long x;
	if (!_BitScanReverse64(&x, n))
		return -1;
#else
	uint64_t x = -1;
	__asm("bsr %0, %1" : "=r"(x) : "r"(n) : );
#endif
	return x;
}

inline int bsf(uint64_t n)
{
#if 0
	return _tzcnt_u64(n);
#else
#if defined(_MSC_VER) && !defined(__clang__)
	unsigned long x;
	if (!_BitScanForward64(&x, n))
		return 64;
#else
	uint64_t x = 64;
	__asm("bsf %0, %1" : "=r"(x) : "r"(n) : );
#endif
	return x;
#endif
}

template< typename T > constexpr bool bt(T x, unsigned n) { return x & T(1) << n; }
template< typename T > constexpr T ceil2_n(T x, unsigned n) { return x + (n - 1) & ~T(n - 1); }
template< typename T > constexpr T ceil_n(T x, unsigned n) { return (x + (n - 1)) / n * n; }
template< typename T > constexpr T floor2_n(T x, unsigned n) { return x & ~T(n - 1); }
template< typename T > constexpr T floor_n(T x, unsigned n) { return x - x % n; }
template< typename T > constexpr int sign(T x, std::false_type is_signed) { return T(0) < x; }
template< typename T > constexpr int sign(T x, std::true_type is_signed) { return (T(0) < x) - (x < T(0)); }
template< typename T > constexpr int sign(T x) { return sign(x, std::is_signed<T>()); }

template< unsigned B, typename T>
constexpr T signext(T x)
{
	const struct { std::make_signed_t<T> i : B; } s = { x };
	return s.i;
}

template< typename T, typename... Ts>
inline constexpr bool any_of_v = std::disjunction_v<std::is_same<T, Ts>...>;

union sf32
{
	uint32_t u;
	int32_t i;
	float f;
	struct
	{
		uint32_t m : 23;
		uint32_t e : 8;
		uint32_t s : 1;
	};
	sf32()
		: u(0)
	{
	}
	sf32(float f)
		: f(f)
	{
	}
	operator float() const
	{
		return f;
	}
	sf32& bit_and(uint32_t x)
	{
		u &= x;
		return *this;
	}
	sf32& bit_or(uint32_t x)
	{
		u |= x;
		return *this;
	}
};

union sf16
{
	uint16_t u;
	int16_t i;
	struct
	{
		uint16_t m : 10;
		uint16_t e : 5;
		uint16_t s : 1;
	};
	sf16()
		: u(0)
	{
	}
	sf16(float f)
	{
		sf32 x(f);
		m = 0;
		e = 0;
		s = x.s;
		x.s = 0;
		if (x.u >= 0x47800000)
			u |= x.u > 0x7f800000 ? 0x7e00 : 0x7c00;
		else
		{
			if (x.u < 0x38800000)
			{
				x.f += .5f;
				u |= x.u - 0x3f000000;
			}
			else
			{
				unsigned t = x.u >> 13 & 1;
				u |= x.u + 0xc8000fff + t >> 13;
			}
		}
	}
	operator float() const
	{
		sf32 x;
		unsigned exp = (u & 0x7c00);
		unsigned sign = (u & 0x8000) << 16;
		unsigned t = (u & 0x7fff) << 13;
		if (exp == 0x7c00)
			x.u = t + 0x70000000;
		else if (exp != 0)
			x.u = t + 0x38000000;
		else
		{
			x.u = t + 0x38800000;
			x.f -= 6.103515625e-05f;
		}
		x.u |= sign;
		return x.f;
	}
};

template< typename T>
inline T saturate_cast(float f)
{
	if constexpr (std::is_signed_v<T>)
		f += sf32(f).bit_and(0x80000000).bit_or(0x3f000000);
	else
		f += .5f;
	using int_t = std::conditional_t<std::is_same_v<T, int32_t> || std::is_same_v<T, uint32_t>, int64_t, int32_t>;
	constexpr int_t min = std::numeric_limits<T>::min();
	constexpr int_t max = std::numeric_limits<T>::max();
	int_t i = static_cast<int_t>(f);
	return static_cast<T>(std::clamp(i, min, max));
}

#if defined VP_ENABLE_INTRIN && !defined VP_DISABLE_SATURATE_CAST_INTRIN

template<>
A_SSE4 inline uint16_t saturate_cast(float f)
{
	__m128i i = _mm_cvtps_epi32(_mm_load_ss(&f));
	return static_cast<uint16_t>(_mm_cvtsi128_si32(_mm_packus_epi32(i, i)));
}

template<>
inline uint8_t saturate_cast(float f)
{
	__m128i i = _mm_cvtps_epi32(_mm_load_ss(&f));
	i = _mm_packs_epi32(i, i);
	return static_cast<uint8_t>(_mm_cvtsi128_si32(_mm_packus_epi16(i, i)));
}
#endif

template< typename T = uint8_t>
inline T* aligned_alloc(size_t count, unsigned alignment)
{
	size_t size = ceil2_n(count * sizeof(T), alignment);
	void* ptr;
#if defined(_WIN32)
	ptr = _aligned_malloc(size, alignment);
#else
	if (posix_memalign(&ptr, alignment, size))
		ptr = 0;
#endif
	return static_cast<T*>(ptr);
}
inline void aligned_free(void* ptr)
{
#if defined(_WIN32)
	_aligned_free(ptr);
#else
	free(ptr);
#endif
}
struct aligned_deleter
{
	void operator()(void* ptr)
	{
		aligned_free(ptr);
	}
};
template< typename T = uint8_t>
using aligned_unique = std::unique_ptr<T, aligned_deleter>;
template< typename T = uint8_t>
inline auto make_aligned_unique(size_t count, unsigned alignment)
{
	return aligned_unique<T>(aligned_alloc<T>(count, alignment), aligned_deleter());
}

struct file_h
{
	int h = -1;
	intmax_t size = 0;
	bool open(const char* path);
	int read(void* buf, size_t nbyte);
	int read(void* buf, intmax_t off, size_t nbyte);
	~file_h();
};

template< typename T, typename... R>
void free(void* ptr, R...)
{
	delete static_cast<T*>(ptr);
}

inline uint64_t xxh_value(const void* ptr, size_t size)
{
	XXH64_state_t state;
	XXH64_reset(&state, 0);
	XXH64_update(&state, ptr, size);
	return XXH64_digest(&state);
}

template< typename ST, typename DT>
struct conv
{
	DT operator ()(ST x);
	static void transform_n(const ST*, size_t, DT*);
};

template<>
void conv<sf16, float>::transform_n(const sf16* s, size_t n, float* d);

template<>
inline float conv<uint8_t, float>::operator()(uint8_t x)
{
	return x * (1 / 255.f);
}

template<>
inline uint8_t conv<float, uint8_t>::operator()(float x)
{
	return saturate_cast<uint8_t>(x * 255.f);
}

template<>
inline float conv<uint16_t, float>::operator()(uint16_t x)
{
	return x * (1 / 65535.f);
}

template<>
inline uint16_t conv<float, uint16_t>::operator()(float x)
{
	return saturate_cast<uint16_t>(x * 65535.f);
}

template< typename ST, typename DT = ST, template< typename, typename > typename F = conv>
inline void copy(size_t dst_h, size_t dst_w, ptrdiff_t dst_stride1, void* __restrict dst_p,
	size_t src_h, size_t src_w, ptrdiff_t src_stride1, const void* src_p)
{
	const ST* src = static_cast<const ST*>(src_p);
	DT* dst = static_cast<DT*>(dst_p);
	ptrdiff_t src_stride = src_stride1 ? src_stride1 / sizeof(ST) : src_w;
	ptrdiff_t dst_stride = dst_stride1 ? dst_stride1 / sizeof(DT) : dst_w;
	for (size_t i = 0; i < src_h; i++)
	{
		auto ptr = src + i * src_stride;
		if constexpr (std::is_same_v<ST, DT>)
			std::copy(ptr, ptr + src_w, dst + i * dst_stride);
		else
			std::transform(ptr, ptr + src_w, dst + i * dst_stride, F<ST, DT>());
	}
}

template< typename T>
inline void pad_replicate(size_t dst_h, size_t dst_w, ptrdiff_t dst_stride1, void* __restrict dst_p,
	size_t src_h, size_t src_w, ptrdiff_t src_stride1, const void* src_p, size_t oy, size_t ox)
{
	T* dst = static_cast<T*>(dst_p);
	ptrdiff_t dst_stride = dst_stride1 ? dst_stride1 / sizeof(T) : dst_w;
	for (size_t i = oy; i < src_h + oy; i++)
	{
		auto ptr = dst + i * dst_stride;
		ptr = std::fill_n(ptr, ox, ptr[ox]) + src_w;
		std::fill_n(ptr, dst_w - src_w - ox, ptr[-1]);
	}
	for (size_t i = 0, y = oy; i < oy; i++)
		std::copy_n(dst + y * dst_stride, dst_w, dst + i * dst_stride);
	for (size_t i = src_h + oy, y = i - 1; i < dst_h; i++)
		std::copy_n(dst + y * dst_stride, dst_w, dst + i * dst_stride);
}

template< typename T>
inline void pad_reflect(size_t dst_h, size_t dst_w, ptrdiff_t dst_stride1, void* __restrict dst_p,
	size_t src_h, size_t src_w, ptrdiff_t src_stride1, const void* src_p, size_t oy, size_t ox)
{
	T* dst = static_cast<T*>(dst_p);
	ptrdiff_t dst_stride = dst_stride1 ? dst_stride1 / sizeof(T) : dst_w;
	for (size_t i = oy; i < src_h + oy; i++)
	{
		auto ptr = dst + i * dst_stride;
		ptr = std::reverse_copy(ptr + ox + 1, ptr + ox + ox + 1, ptr) + src_w;
		std::reverse_copy(ptr - dst_w + src_w + ox - 1, ptr - 1, ptr);
	}
	for (size_t i = 0, y = oy + oy; i < oy; i++, y--)
		std::copy_n(dst + y * dst_stride, dst_w, dst + i * dst_stride);
	for (size_t i = src_h + oy, y = i - 2; i < dst_h; i++, y--)
		std::copy_n(dst + y * dst_stride, dst_w, dst + i * dst_stride);
}

template< typename ST, typename DT = ST, template< typename, typename > typename F = conv>
inline void copy_pad_replicate(size_t dst_h, size_t dst_w, ptrdiff_t dst_stride1, void* __restrict dst_p,
	size_t src_h, size_t src_w, ptrdiff_t src_stride1, const void* src_p, size_t oy, size_t ox)
{
	DT* dst = static_cast<DT*>(dst_p);
	ptrdiff_t dst_stride = dst_stride1 ? dst_stride1 / sizeof(DT) : dst_w;
	copy<ST, DT, F>(dst_h, dst_w, dst_stride1, dst + oy * dst_stride + ox, src_h, src_w, src_stride1, src_p);
	pad_replicate<DT>(dst_h, dst_w, dst_stride1, dst_p, src_h, src_w, src_stride1, src_p, oy, ox);
}

template< typename ST, typename DT = ST, template< typename, typename > typename F = conv>
inline void copy_pad_reflect(size_t dst_h, size_t dst_w, ptrdiff_t dst_stride1, void* __restrict dst_p,
	size_t src_h, size_t src_w, ptrdiff_t src_stride1, const void* src_p, size_t oy, size_t ox)
{
	DT* dst = static_cast<DT*>(dst_p);
	ptrdiff_t dst_stride = dst_stride1 ? dst_stride1 / sizeof(DT) : dst_w;
	copy<ST, DT, F>(dst_h, dst_w, dst_stride1, dst + oy * dst_stride + ox, src_h, src_w, src_stride1, src_p);
	pad_reflect<DT>(dst_h, dst_w, dst_stride1, dst_p, src_h, src_w, src_stride1, src_p, oy, ox);
}

typedef void (*copy_t)(size_t dst_h, size_t dst_w, ptrdiff_t dst_stride1, void* __restrict dst_p,
	size_t src_h, size_t src_w, ptrdiff_t src_stride1, const void* src_p);
typedef void (*copy_pad_t)(size_t dst_h, size_t dst_w, ptrdiff_t dst_stride1, void* __restrict dst_p,
	size_t src_h, size_t src_w, ptrdiff_t src_stride1, const void* src_p, size_t oy, size_t ox);

template< class V, class F, class K>
struct vfk_s
{
	static_assert(std::is_integral_v<F> | std::is_enum_v<F>, "");
	V v;
	F f;
	K k;
	template< class... As>
	constexpr vfk_s(V v, F f, As... kargs)
		: v(v)
		, f(f)
		, k(std::forward<As>(kargs)...)
	{
	}
	template< class IT>
	static constexpr V select(IT first, IT last, int64_t f, K k)
	{
		for (; first != last; first++)
		{
			if ((f & first->f) == first->f && first->k == k)
				return first->v;
		}
		return nullptr;
	}
};

static std::string format(const char* f, ...)
{
	va_list arg;
	va_start(arg, f);
	int n = std::vsnprintf(nullptr, 0, f, arg);
	if (n < 0)
	{
		va_end(arg);
		throw std::runtime_error("vsnprintf() failed");
	}
	std::string str(++n, 0);
	std::vsnprintf(str.data(), str.size(), f, arg);
	va_end(arg);
	return str;
}

}
