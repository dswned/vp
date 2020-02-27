#pragma once
#include "config.h"

#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#else
#include <immintrin.h>
#endif
#include <cstdint>
#include <cstring>
#include <utility>

namespace uf {
inline int bsr(uint64_t n)
{
#if defined(_MSC_VER) && !defined(__clang__)
	unsigned long x = -1;
	if (n)
		_BitScanReverse64(&x, n);
#else
	uint64_t x = -1;
	if (n)
		__asm("bsr %0, %1" : "=r"(x) : "r"(n) : );
#endif
	return x;
}
constexpr int bsr_const(uint64_t n)
{
	if (!n)
		return -1;
	uint64_t a = n, b = 0, j = 64, k = 0;
	do
	{
		j >>= 1;
		k = 1ULL << j;
		if (a >= k)
		{
			a >>= j;
			b += j;
		}
	} while (j > 0);
	return b;
}
template< typename T > constexpr int sign(T x, std::false_type is_signed) { return T(0) < x; }
template< typename T > constexpr int sign(T x, std::true_type is_signed) { return (T(0) < x) - (x < T(0)); }
template< typename T > constexpr int sign(T x) { return sign(x, std::is_signed<T>()); }
template< typename T > constexpr bool bt(T x, int n) { return x & T(1) << n; }
template< unsigned B, typename T >
constexpr T signext(T x)
{
	const struct { std::make_signed_t<T> x : B; } s = { x };
	return s.x;
}
template< typename T > constexpr T ceil_n(T x, T n) { return x + n - 1 & ~(n - 1); }
template< typename T >
void shr(T* data, size_t size, size_t n)
{
	constexpr size_t wsize = sizeof(T) * 8;
	size_t wshift = n / wsize;
	if (wshift)
	{
		for (size_t i = 0; i < size; ++i)
			data[i] = wshift < size - i ? data[i + wshift] : 0;
	}
	if (n %= wsize)
	{
		for (size_t i = 0; i < size - 1; ++i)
			data[i] = (data[i] >> n) | (data[i + 1] << (wsize - n));
		data[size - 1] >>= n;
	}
}
template< typename ResultType, typename... Args>
struct function_traits_s
{
	static constexpr size_t arity = sizeof...(Args);
	using result_type = ResultType;
	template< size_t i>
	struct arg
	{
		using type = typename std::tuple_element<i, std::tuple<Args...>>::type;
	};
};
template< typename T >
struct function_traits_impl {};
template< typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(Args...)>
	: function_traits_s<ReturnType, Args...> {};
template< typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(*)(Args...)>
	: function_traits_s<ReturnType, Args...> {};
template< typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...)>
	: function_traits_s<ReturnType, Args...> {};
template< typename ClassType, typename ReturnType, typename... Args>
struct function_traits_impl<ReturnType(ClassType::*)(Args...) const>
	: function_traits_s<ReturnType, Args...> {};
template< typename T, typename V = void>
struct function_traits
	: function_traits_impl<T> {};
template< typename T >
struct function_traits<T, decltype((void)&T::operator())>
	: function_traits_impl<decltype(&T::operator())> {};
}