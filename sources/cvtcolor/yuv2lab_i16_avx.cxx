#include "cvtcolor.h"
#if defined __GNUC__ || defined __clang__
#include <immintrin.h>
#else
#include <intrin.h>
#endif
#include <cmath>

namespace {

void spline_build(size_t n, const float* f, float* tab)
{
	float cn = 0;
	tab[0] = tab[1] = 0;
	for (size_t i = 1; i < n; i++)
	{
		float t = (f[i + 1] - f[i] * 2 + f[i - 1]) * 3;
		float l = 1.f / (4 - tab[(i - 1) * 4]);
		tab[i * 4] = l;
		tab[i * 4 + 1] = (t - tab[(i - 1) * 4 + 1]) * l;
	}
	for (size_t i = n - 1; i != -1; i--)
	{
		float c = tab[i * 4 + 1] - tab[i * 4] * cn;
		float b = f[i + 1] - f[i] - (cn + c * 2) / 3;
		float d = (cn - c) / 3;
		tab[i * 4] = f[i];
		tab[i * 4 + 1] = b;
		tab[i * 4 + 2] = c;
		tab[i * 4 + 3] = d;
		cn = c;
	}
}

inline void transpose4x4(__m128& b0, __m128& b1, __m128& b2, __m128& b3,
	const __m128& a0, const __m128& a1, const __m128& a2, const __m128& a3)
{
	__m128i t0 = _mm_castps_si128(_mm_unpacklo_ps(a0, a1));
	__m128i t1 = _mm_castps_si128(_mm_unpacklo_ps(a2, a3));
	__m128i t2 = _mm_castps_si128(_mm_unpackhi_ps(a0, a1));
	__m128i t3 = _mm_castps_si128(_mm_unpackhi_ps(a2, a3));
	b0 = _mm_castsi128_ps(_mm_unpacklo_epi64(t0, t1));
	b1 = _mm_castsi128_ps(_mm_unpackhi_epi64(t0, t1));
	b2 = _mm_castsi128_ps(_mm_unpacklo_epi64(t2, t3));
	b3 = _mm_castsi128_ps(_mm_unpackhi_epi64(t2, t3));
}

inline __m128 spline_interpolate(const __m128& x, const float* tab)
{
	__m128i ix = _mm_and_si128(_mm_cvttps_epi32(x), _mm_set1_epi32(LAB_CBRT_TAB_SIZE - 1));
	__m128 xx = _mm_sub_ps(x, _mm_cvtepi32_ps(ix));
	ix = _mm_slli_epi32(ix, 2);
	__m128 t[4], tt[4];
	t[0] = _mm_loadu_ps(tab + _mm_cvtsi128_si32(ix));
	t[1] = _mm_loadu_ps(tab + _mm_extract_epi32(ix, 1));
	t[2] = _mm_loadu_ps(tab + _mm_extract_epi32(ix, 2));
	t[3] = _mm_loadu_ps(tab + _mm_extract_epi32(ix, 3));
	transpose4x4(tt[0], tt[1], tt[2], tt[3], t[0], t[1], t[2], t[3]);
	return _mm_add_ps(_mm_mul_ps(_mm_add_ps(_mm_mul_ps(_mm_add_ps(_mm_mul_ps(tt[3], xx), tt[2]), xx), tt[1]), xx), tt[0]);
}

}

VP_EXTERN_C void yuv2lab_i16_avx(size_t n, void* lab_l, void* lab_a, void* lab_b, void* yuv_y, void* yuv_u, void* yuv_v)
{
	const float* coeffs = lab_yuv_i16_data, * tab = lab_cbrt_tab;
	const __m256 k_sa = _mm256_broadcast_ss(&::k_sa),
		k_sb = _mm256_broadcast_ss(&::k_sb),
		k_116_100 = _mm256_broadcast_ss(&::k_116_100),
		k_16_100 = _mm256_broadcast_ss(&::k_16_100),
		k_512_65535 = _mm256_broadcast_ss(&::k_512_65535),
		k_1_65536 = _mm256_broadcast_ss(&::k_1_65536),
		k_65535 = _mm256_broadcast_ss(&::k_65535),
		k_65536 = _mm256_broadcast_ss(&::k_65536),
		k_32768 = _mm256_broadcast_ss(&::k_32768),
		k_half = _mm256_broadcast_ss(&::k_half);
	const __m128i k_zero = _mm_setzero_si128();
	for (size_t i = 0; i < n; i += 16)
	{
		__m128i y_8i16 = _mm_loadu_si128((__m128i const*)((uint8_t*)yuv_y + i));
		__m128i u_8i16 = _mm_loadu_si128((__m128i const*)((uint8_t*)yuv_u + i));
		__m128i v_8i16 = _mm_loadu_si128((__m128i const*)((uint8_t*)yuv_v + i));
		__m256i y_i, u_i, v_i;
		y_i = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_unpacklo_epi16(y_8i16, k_zero)), _mm_unpackhi_epi16(y_8i16, k_zero), 1);
		u_i = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_unpacklo_epi16(u_8i16, k_zero)), _mm_unpackhi_epi16(u_8i16, k_zero), 1);
		v_i = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_unpacklo_epi16(v_8i16, k_zero)), _mm_unpackhi_epi16(v_8i16, k_zero), 1);
		__m256 y, u, v;
		y = _mm256_mul_ps(_mm256_cvtepi32_ps(y_i), k_512_65535);
		u = _mm256_sub_ps(_mm256_mul_ps(_mm256_cvtepi32_ps(u_i), k_1_65536), k_half);
		v = _mm256_sub_ps(_mm256_mul_ps(_mm256_cvtepi32_ps(v_i), k_1_65536), k_half);
		__m256 X, Z;
		X = _mm256_add_ps(_mm256_add_ps(y, _mm256_mul_ps(u, _mm256_broadcast_ss(coeffs + 1))), _mm256_mul_ps(v, _mm256_broadcast_ss(coeffs + 2)));
		Z = _mm256_add_ps(_mm256_add_ps(y, _mm256_mul_ps(u, _mm256_broadcast_ss(coeffs + 7))), _mm256_mul_ps(v, _mm256_broadcast_ss(coeffs + 8)));
		__m128 fx0, fx1, fy0, fy1, fz0, fz1;
		fx0 = spline_interpolate(_mm256_castps256_ps128(X), tab);
		fx1 = spline_interpolate(_mm256_extractf128_ps(X, 1), tab);
		fy0 = spline_interpolate(_mm256_castps256_ps128(y), tab);
		fy1 = spline_interpolate(_mm256_extractf128_ps(y, 1), tab);
		fz0 = spline_interpolate(_mm256_castps256_ps128(Z), tab);
		fz1 = spline_interpolate(_mm256_extractf128_ps(Z, 1), tab);
		__m256 fx, fy, fz;
		fx = _mm256_insertf128_ps(_mm256_castps128_ps256(fx0), fx1, 1);
		fy = _mm256_insertf128_ps(_mm256_castps128_ps256(fy0), fy1, 1);
		fz = _mm256_insertf128_ps(_mm256_castps128_ps256(fz0), fz1, 1);
		__m256 L, a, b;
		L = _mm256_sub_ps(_mm256_mul_ps(fy, k_116_100), k_16_100);
		a = _mm256_mul_ps(_mm256_sub_ps(fx, fy), k_sa);
		b = _mm256_mul_ps(_mm256_sub_ps(fy, fz), k_sb);
		__m256i L_i, a_i, b_i;
		L_i = _mm256_cvtps_epi32(_mm256_mul_ps(L, k_65535));
		a_i = _mm256_cvtps_epi32(_mm256_add_ps(_mm256_mul_ps(a, k_65536), k_32768));
		b_i = _mm256_cvtps_epi32(_mm256_add_ps(_mm256_mul_ps(b, k_65536), k_32768));
		_mm_storeu_si128((__m128i*)((uint8_t*)lab_l + i), _mm_packus_epi32(_mm256_castsi256_si128(L_i), _mm256_extractf128_si256(L_i, 1)));
		_mm_storeu_si128((__m128i*)((uint8_t*)lab_a + i), _mm_packus_epi32(_mm256_castsi256_si128(a_i), _mm256_extractf128_si256(a_i, 1)));
		_mm_storeu_si128((__m128i*)((uint8_t*)lab_b + i), _mm_packus_epi32(_mm256_castsi256_si128(b_i), _mm256_extractf128_si256(b_i, 1)));
	}
}

VP_EXTERN_C void init_lab_cbrt_tab()
{
	static bool initialized = 0;
	if (!initialized)
	{
		float f[LAB_CBRT_TAB_SIZE + 1];
		const float tscale = 2. / LAB_CBRT_TAB_SIZE;
		const float lthresh = 216 / 24389.f;
		const float lscale = 841 / 108.f;
		const float lbias = 16 / 116.f;
		for (int i = 0; i <= LAB_CBRT_TAB_SIZE; i++)
		{
			float x = tscale * i;
			f[i] = x > lthresh ? std::cbrt(x) : x * lscale + lbias;
		}
		spline_build(LAB_CBRT_TAB_SIZE, f, lab_cbrt_tab);
		initialized = true;
	}
}
