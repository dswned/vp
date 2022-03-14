#include "cvtcolor.h"
#if defined __GNUC__ || defined __clang__
#include <immintrin.h>
#else
#include <intrin.h>
#endif

VP_EXTERN_C void lab2yuv_i16_avx(size_t n, void* yuv_y, void* yuv_u, void* yuv_v, void* lab_l, void* lab_a, void* lab_b)
{
	const float* coeffs = lab_yuv_i16_data + 12;
	const __m256 k_1_65535 = _mm256_broadcast_ss(&::k_1_65535),
		k_1_65536 = _mm256_broadcast_ss(&::k_1_65536),
		k_65535 = _mm256_broadcast_ss(&::k_65535),
		k_65536 = _mm256_broadcast_ss(&::k_65536),
		k_32768 = _mm256_broadcast_ss(&::k_32768),
		k_half = _mm256_broadcast_ss(&::k_half),
		k_100_116 = _mm256_broadcast_ss(&::k_100_116),
		k_16_116 = _mm256_broadcast_ss(&::k_16_116),
		k_sai = _mm256_broadcast_ss(&::k_sai),
		k_sbi = _mm256_broadcast_ss(&::k_sbi),
		k_2700_24389 = _mm256_broadcast_ss(&::k_2700_24389),
		k_008 = _mm256_broadcast_ss(&::k_008),
		k_6_29 = _mm256_broadcast_ss(&::k_6_29),
		k_sai_3132_24389 = _mm256_broadcast_ss(&::k_sai_3132_24389),
		k_sbi_3132_24389 = _mm256_broadcast_ss(&::k_sbi_3132_24389);
	const __m128i k_zero = _mm_setzero_si128();
	for (size_t i = 0; i < n; i += 16)
	{
		__m128i L_8i16 = _mm_loadu_si128((__m128i const*)((uint8_t*)lab_l + i));
		__m128i a_8i16 = _mm_loadu_si128((__m128i const*)((uint8_t*)lab_a + i));
		__m128i b_8i16 = _mm_loadu_si128((__m128i const*)((uint8_t*)lab_b + i));
		__m256i L_i, a_i, b_i;
		L_i = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_unpacklo_epi16(L_8i16, k_zero)), _mm_unpackhi_epi16(L_8i16, k_zero), 1);
		a_i = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_unpacklo_epi16(a_8i16, k_zero)), _mm_unpackhi_epi16(a_8i16, k_zero), 1);
		b_i = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_unpacklo_epi16(b_8i16, k_zero)), _mm_unpackhi_epi16(b_8i16, k_zero), 1);
		__m256 L, a, b;
		L = _mm256_mul_ps(_mm256_cvtepi32_ps(L_i), k_1_65535);
		a = _mm256_sub_ps(_mm256_mul_ps(_mm256_cvtepi32_ps(a_i), k_1_65536), k_half);
		b = _mm256_sub_ps(_mm256_mul_ps(_mm256_cvtepi32_ps(b_i), k_1_65536), k_half);
		__m256 xf, yf, zf;
		yf = _mm256_add_ps(_mm256_mul_ps(L, k_100_116), k_16_116);
		xf = _mm256_add_ps(yf, _mm256_mul_ps(a, k_sai));
		zf = _mm256_add_ps(yf, _mm256_mul_ps(b, k_sbi));
		__m256 X, Y, Z, fy;
		fy = _mm256_mul_ps(L, k_2700_24389);
		Y = _mm256_blendv_ps(fy, _mm256_mul_ps(_mm256_mul_ps(yf, yf), yf), _mm256_cmp_ps(L, k_008, _CMP_GT_OQ));
		X = _mm256_blendv_ps(_mm256_add_ps(fy, _mm256_mul_ps(a, k_sai_3132_24389)), _mm256_mul_ps(_mm256_mul_ps(xf, xf), xf), _mm256_cmp_ps(xf, k_6_29, _CMP_GT_OQ));
		Z = _mm256_blendv_ps(_mm256_add_ps(fy, _mm256_mul_ps(b, k_sbi_3132_24389)), _mm256_mul_ps(_mm256_mul_ps(zf, zf), zf), _mm256_cmp_ps(zf, k_6_29, _CMP_GT_OQ));
		__m256 u, v;
		u = _mm256_add_ps(_mm256_add_ps(
			_mm256_mul_ps(_mm256_broadcast_ss(coeffs + 3), X),
			_mm256_mul_ps(_mm256_broadcast_ss(coeffs + 4), Y)),
			_mm256_mul_ps(_mm256_broadcast_ss(coeffs + 5), Z));
		v = _mm256_add_ps(_mm256_add_ps(
			_mm256_mul_ps(_mm256_broadcast_ss(coeffs + 6), X),
			_mm256_mul_ps(_mm256_broadcast_ss(coeffs + 7), Y)),
			_mm256_mul_ps(_mm256_broadcast_ss(coeffs + 8), Z));
		__m256i y_i, u_i, v_i;
		y_i = _mm256_cvtps_epi32(_mm256_mul_ps(Y, k_65535));
		u_i = _mm256_cvtps_epi32(_mm256_add_ps(_mm256_mul_ps(u, k_65536), k_32768));
		v_i = _mm256_cvtps_epi32(_mm256_add_ps(_mm256_mul_ps(v, k_65536), k_32768));
		_mm_storeu_si128((__m128i*)((uint8_t*)yuv_y + i), _mm_packus_epi32(_mm256_castsi256_si128(y_i), _mm256_extractf128_si256(y_i, 1)));
		_mm_storeu_si128((__m128i*)((uint8_t*)yuv_u + i), _mm_packus_epi32(_mm256_castsi256_si128(u_i), _mm256_extractf128_si256(u_i, 1)));
		_mm_storeu_si128((__m128i*)((uint8_t*)yuv_v + i), _mm_packus_epi32(_mm256_castsi256_si128(v_i), _mm256_extractf128_si256(v_i, 1)));
	}
}
