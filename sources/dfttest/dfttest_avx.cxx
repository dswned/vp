#include "config.h"
#if defined(BUILD_DFTTEST)

#include "dfttest.h"

#define INSTRSET 7
#include <vectorclass2/vectorclass.h>
#include <vectorclass2/vectormath_exp.h>

#define ORG

namespace dfttest {

namespace {

template< typename T>
inline void load(float* __restrict d, const float* w, size_t sbsize, const T* src, ptrdiff_t src_stride);

template<>
inline void load(float* __restrict d, const float* w, size_t sbsize, const uint8_t* src, ptrdiff_t src_stride)
{
	__m128 s0, s1;
	for (int i = 0; i < sbsize; i++)
	{
		for (int j = 0; j < sbsize; j += 4)
		{
			s0 = _mm_cvtepi32_ps(_mm_cvtepu8_epi32(_mm_cvtsi32_si128(*(const int*)(src + j))));
			s1 = _mm_loadu_ps(w + j);
			_mm_storeu_ps(d + j, _mm_mul_ps(s0, s1));
		}
		d += sbsize;
		w += sbsize;
		src += src_stride;
	}
}

template<>
inline void load(float* __restrict d, const float* w, size_t sbsize, const uint16_t* src, ptrdiff_t src_stride)
{
	__m128 s0, s1, scale = _mm_broadcast_ss(&scale_0<uint16_t>);
	for (int i = 0; i < sbsize; i++)
	{
		for (int j = 0; j < sbsize; j += 4)
		{
			s0 = _mm_cvtepi32_ps(_mm_cvtepu16_epi32(_mm_loadl_epi64((const __m128i*)(src + j))));
			s1 = _mm_loadu_ps(w + j);
#if 1
			_mm_storeu_ps(d + j, _mm_mul_ps(_mm_mul_ps(s1, scale), s0));
#else
			_mm_storeu_ps(d + j, _mm_mul_ps(s1, s0));
#endif
		}
		d += sbsize;
		w += sbsize;
		src += src_stride;
	}
}

template<>
inline void load(float* __restrict d, const float* w, size_t sbsize, const float* src, ptrdiff_t src_stride)
{
	__m128 s0, s1, scale = _mm_broadcast_ss(&scale_0<float>);
	for (int i = 0; i < sbsize; i++)
	{
		for (int j = 0; j < sbsize; j += 4)
		{
			s0 = _mm_loadu_ps(src + j);
			s1 = _mm_loadu_ps(w + j);
			_mm_storeu_ps(d + j, _mm_mul_ps(_mm_mul_ps(s1, s0), scale));
		}
		d += sbsize;
		w += sbsize;
		src += src_stride;
	}
}

inline void store(float* __restrict dst, ptrdiff_t dst_stride, const float* s, const float* w, size_t sbsize)
{
	Vec4f s0, s1;
	if (int n = sbsize & 3; !n)
	{
		for (int i = 0; i < sbsize; i++)
		{
			for (int j = 0; j < sbsize; j += 4)
			{
				s0.load(s + j);
				s1.load(w + j);
				mul_add(s0, s1, Vec4f().load(dst + j)).store(dst + j);
			}
			dst += dst_stride;
			s += sbsize;
			w += sbsize;
		}
	}
	else
	{
		for (int i = 0, j; i < sbsize; i++)
		{
			for (j = 0; j < sbsize - 4; j += 4)
			{
				s0.load(s + j);
				s1.load(w + j);
				mul_add(s0, s1, Vec4f().load(dst + j)).store(dst + j);
			}
			s0.load(s + j);
			s1.load(w + j);
			mul_add(s0, s1, Vec4f().load(dst + j)).store_partial(n, dst + j);
			dst += dst_stride;
			s += sbsize;
			w += sbsize;
		}
	}
}

template< typename T>
void cast(T* __restrict dst, size_t dst_h, size_t dst_w, ptrdiff_t dst_stride, const float* src, ptrdiff_t src_stride);

template<>
void cast(uint8_t* __restrict dst, size_t dst_h, size_t dst_w, ptrdiff_t dst_stride, const float* src, ptrdiff_t src_stride)
{
	for (int i = 0; i < dst_h; i++)
	{
		for (int j = 0; j < dst_w; j += 16)
		{
			__m128i i4_0 = _mm_cvtps_epi32(_mm_loadu_ps(src + j));
			__m128i i4_1 = _mm_cvtps_epi32(_mm_loadu_ps(src + j + 4));
			__m128i i4_2 = _mm_cvtps_epi32(_mm_loadu_ps(src + j + 8));
			__m128i i4_3 = _mm_cvtps_epi32(_mm_loadu_ps(src + j + 12));
			__m128i s8_0 = _mm_packs_epi32(i4_0, i4_1);
			__m128i s8_1 = _mm_packs_epi32(i4_2, i4_3);
			_mm_storeu_ps((float*)(dst + j), _mm_castsi128_ps(_mm_packus_epi16(s8_0, s8_1)));
		}
		dst += dst_stride;
		src += src_stride;
	}
}

template<>
void cast(uint16_t* __restrict dst, size_t dst_h, size_t dst_w, ptrdiff_t dst_stride, const float* src, ptrdiff_t src_stride)
{
	__m128 scale = _mm_broadcast_ss(&scale_1<uint16_t>);
	for (int i = 0; i < dst_h; i++)
	{
		for (int j = 0; j < dst_w; j += 8)
		{
#ifndef ORG
#if 1
			__m128i i4_0 = _mm_cvtps_epi32(_mm_mul_ps(_mm_loadu_ps(src + j), scale));
			__m128i i4_1 = _mm_cvtps_epi32(_mm_mul_ps(_mm_loadu_ps(src + j + 4), scale));
#else
			__m128i i4_0 = _mm_cvtps_epi32(_mm_loadu_ps(src + j));
			__m128i i4_1 = _mm_cvtps_epi32(_mm_loadu_ps(src + j + 4));
#endif
#else
			__m128i i4_0 = _mm_cvttps_epi32(_mm_add_ps(_mm_mul_ps(_mm_loadu_ps(src + j), scale), _mm_set1_ps(.5f)));
			__m128i i4_1 = _mm_cvttps_epi32(_mm_add_ps(_mm_mul_ps(_mm_loadu_ps(src + j + 4), scale), _mm_set1_ps(.5f)));
#endif
			_mm_storeu_ps((float*)(dst + j), _mm_castsi128_ps(_mm_packus_epi32(i4_0, i4_1)));
		}
		dst += dst_stride;
		src += src_stride;
	}
}

template<>
void cast(float* __restrict dst, size_t dst_h, size_t dst_w, ptrdiff_t dst_stride, const float* src, ptrdiff_t src_stride)
{
	__m256 scale = _mm256_broadcast_ss(&scale_1<float>);
	for (int i = 0; i < dst_h; i++)
	{
		for (int j = 0; j < dst_w; j += 8)
		{
			__m256 f8 = _mm256_mul_ps(_mm256_loadu_ps(src + j), scale);
			_mm256_storeu_ps(dst + j, f8);
		}
		dst += dst_stride;
		src += src_stride;
	}
}

void remove_mean(size_t ccnt, fftwf_complex* __restrict dftc, const fftwf_complex* dftgc, fftwf_complex* __restrict dftc2)
{
	__m256 c0 = _mm256_loadu_ps((float*)dftc);
	__m256 c1 = _mm256_loadu_ps((float*)dftgc);
	__m128 g4 = _mm_div_ss(_mm256_castps256_ps128(c0), _mm256_castps256_ps128(c1));
	g4 = _mm_shuffle_ps(g4, g4, 0);
	__m256 g = _mm256_insertf128_ps(_mm256_castps128_ps256(g4), g4, 1);
	for (int i = 0;;)
	{
		c1 = _mm256_mul_ps(c1, g);
		_mm256_storeu_ps((float*)(dftc2 + i), c1);
		_mm256_storeu_ps((float*)(dftc + i), _mm256_sub_ps(c0, c1));
		if (i = i + 4; i >= ccnt)
			break;
		c0 = _mm256_loadu_ps((float*)(dftc + i));
		c1 = _mm256_loadu_ps((float*)(dftgc + i));
	}
}

void add_mean(size_t ccnt, fftwf_complex* __restrict dftc, const fftwf_complex* dftc2)
{
	for (int i = 0; i < ccnt; i += 4)
	{
		__m256 c0 = _mm256_loadu_ps((float*)(dftc + i));
		__m256 c1 = _mm256_loadu_ps((float*)(dftc2 + i));
		_mm256_storeu_ps((float*)(dftc + i), _mm256_add_ps(c0, c1));
	}
}

}

void filter_0_avx(size_t ccnt, fftwf_complex* __restrict dftc, const float* sigmas, const void*)
{
	for (int i = 0; i < ccnt; i += 4)
	{
		__m128 lo = _mm_loadu_ps((float*)(dftc + i));
		__m128 hi = _mm_loadu_ps((float*)(dftc + i + 2));
		__m128 real = _mm_shuffle_ps(lo, hi, _MM_SHUFFLE(2, 0, 2, 0));
		__m128 imag = _mm_shuffle_ps(lo, hi, _MM_SHUFFLE(3, 1, 3, 1));
		__m128 psd = _mm_add_ps(_mm_mul_ps(real, real), _mm_mul_ps(imag, imag));
#ifndef ORG
		__m128 diff = _mm_sub_ps(psd, _mm_loadu_ps(sigmas + i));
		__m128 mult = _mm_mul_ps(_mm_max_ps(diff, _mm_setzero_ps()), _mm_rcp_ps(psd));
#else
		__m128 diff = _mm_sub_ps(psd, _mm_loadu_ps(sigmas + i));
		__m128 mult = _mm_div_ps(_mm_max_ps(diff, _mm_setzero_ps()), _mm_add_ps(psd, _mm_set1_ps(1e-15f)));
#endif
		real = _mm_mul_ps(real, mult);
		imag = _mm_mul_ps(imag, mult);
		lo = _mm_unpacklo_ps(real, imag);
		hi = _mm_unpackhi_ps(real, imag);
		_mm_storeu_ps((float*)(dftc + i), lo);
		_mm_storeu_ps((float*)(dftc + i + 2), hi);
	}
}

void filter_0_sqrt_avx(size_t ccnt, fftwf_complex* __restrict dftc, const float* sigmas, const void*)
{
	for (int i = 0; i < ccnt; i += 4)
	{
		__m128 lo = _mm_loadu_ps((float*)(dftc + i));
		__m128 hi = _mm_loadu_ps((float*)(dftc + i + 2));
		__m128 real = _mm_shuffle_ps(lo, hi, _MM_SHUFFLE(2, 0, 2, 0));
		__m128 imag = _mm_shuffle_ps(lo, hi, _MM_SHUFFLE(3, 1, 3, 1));
		__m128 psd = _mm_add_ps(_mm_mul_ps(real, real), _mm_mul_ps(imag, imag));
		__m128 diff = _mm_sub_ps(psd, _mm_loadu_ps(sigmas + i));
		__m128 mult = _mm_mul_ps(_mm_max_ps(diff, _mm_setzero_ps()), _mm_rcp_ps(psd));
		mult = _mm_sqrt_ps(mult);
		real = _mm_mul_ps(real, mult);
		imag = _mm_mul_ps(imag, mult);
		lo = _mm_unpacklo_ps(real, imag);
		hi = _mm_unpackhi_ps(real, imag);
		_mm_storeu_ps((float*)(dftc + i), lo);
		_mm_storeu_ps((float*)(dftc + i + 2), hi);
	}
}

void filter_0_pow_avx(size_t ccnt, fftwf_complex* __restrict dftc, const float* sigmas, const void* beta)
{
	for (int i = 0; i < ccnt; i += 4)
	{
		__m128 lo = _mm_loadu_ps((float*)(dftc + i));
		__m128 hi = _mm_loadu_ps((float*)(dftc + i + 2));
		__m128 real = _mm_shuffle_ps(lo, hi, _MM_SHUFFLE(2, 0, 2, 0));
		__m128 imag = _mm_shuffle_ps(lo, hi, _MM_SHUFFLE(3, 1, 3, 1));
		__m128 psd = _mm_add_ps(_mm_mul_ps(real, real), _mm_mul_ps(imag, imag));
		__m128 diff = _mm_sub_ps(psd, _mm_loadu_ps(sigmas + i));
		__m128 mult = _mm_mul_ps(_mm_max_ps(diff, _mm_setzero_ps()), _mm_rcp_ps(psd));
		mult = pow(Vec4f(mult), *reinterpret_cast<const float*>(beta));
		real = _mm_mul_ps(real, mult);
		imag = _mm_mul_ps(imag, mult);
		lo = _mm_unpacklo_ps(real, imag);
		hi = _mm_unpackhi_ps(real, imag);
		_mm_storeu_ps((float*)(dftc + i), lo);
		_mm_storeu_ps((float*)(dftc + i + 2), hi);
	}
}

void filter_1_avx(size_t ccnt, fftwf_complex* __restrict dftc, const float* sigmas, const void*)
{
	for (int i = 0; i < ccnt; i += 4)
	{
		__m128 lo = _mm_loadu_ps((float*)(dftc + i));
		__m128 hi = _mm_loadu_ps((float*)(dftc + i + 2));
		__m128 real = _mm_shuffle_ps(lo, hi, _MM_SHUFFLE(2, 0, 2, 0));
		__m128 imag = _mm_shuffle_ps(lo, hi, _MM_SHUFFLE(3, 1, 3, 1));
		__m128 psd = _mm_add_ps(_mm_mul_ps(real, real), _mm_mul_ps(imag, imag));
#if 1
		__m128 mask = _mm_cmp_ps(psd, _mm_loadu_ps(sigmas + i), _CMP_GE_OS);
		real = _mm_and_ps(mask, real);
		imag = _mm_and_ps(mask, imag);
#else
		__m128i mask = _mm_cmplt_epi32(_mm_castps_si128(psd), _mm_castps_si128(_mm_loadu_ps(sigmas + i)));
		real = _mm_castsi128_ps(_mm_andnot_si128(mask, _mm_castps_si128(real)));
		imag = _mm_castsi128_ps(_mm_andnot_si128(mask, _mm_castps_si128(imag)));
#endif
		lo = _mm_unpacklo_ps(real, imag);
		hi = _mm_unpackhi_ps(real, imag);
		_mm_storeu_ps((float*)(dftc + i), lo);
		_mm_storeu_ps((float*)(dftc + i + 2), hi);
	}
}

template< typename T>
void proc_avx(const filter* d, context* ctx)
{
	fftwf_complex* dftc = ctx->dftc[0].get(), * dftc2 = ctx->dftc[1].get();
	uint8_t* tmp_p = ctx->tmp_p.front();
	float* dftr = ctx->dftr.get(), * buf = ctx->buf.front();
	for (int p = 0; p < d->np; p++)
	{
		auto& srcp = ctx->fb.front()->plane[p], & dstp = ctx->fb.back()->plane[p];
		if (!d->process[p])
			srcp.copy_to(dstp);
		else
		{
			size_t tmp_h = d->padded_h[p], tmp_w = d->padded_w[p];
			size_t dst_h = dstp.h, dst_w = dstp.w, oy = (tmp_h - dst_h) / 2, ox = (tmp_w - dst_w) / 2;
			uf::copy_pad_reflect<T>(tmp_h, tmp_w, tmp_w * sizeof(T), tmp_p, dst_h, dst_w, srcp.stride, srcp.p, oy, ox);
			std::fill_n(buf, tmp_h * tmp_w, 0.f);
			for (int y = 0; y <= tmp_h - d->sbsize; y += d->step)
			{
				for (int x = 0; x <= tmp_w - d->sbsize; x += d->step)
				{
					load(dftr, d->hw, d->sbsize, reinterpret_cast<T*>(tmp_p) + y * tmp_w + x, tmp_w);
					fftwf_execute_dft_r2c(d->ft, dftr, dftc);
					if (d->zmean)
						remove_mean(d->ccnt, dftc, d->dftgc, dftc2);
					d->filter_coeffs(d->ccnt, dftc, d->sigmas[p], &d->f0beta);
					if (d->zmean)
						add_mean(d->ccnt, dftc, dftc2);
					fftwf_execute_dft_c2r(d->fti, dftc, dftr);
					if (d->smode)
						store(buf + y * tmp_w + x, tmp_w, dftr, d->hw, d->sbsize);
					else
						buf[y * tmp_w + x] = dftr[oy * d->sbsize + ox] * d->hw[oy * d->sbsize + ox];
				}
			}
			cast(reinterpret_cast<T*>(dstp.p), dst_h, dst_w, dstp.stride / sizeof(T), buf + (oy * tmp_w + ox) * !!d->smode, tmp_w);
		}
	}
}

template void proc_avx<uint8_t>(const filter*, context*);
template void proc_avx<uint16_t>(const filter*, context*);
template void proc_avx<float>(const filter*, context*);

template< typename T>
void proc_temporal_avx(const filter* d, context* ctx, int n)
{
	int k = !d->tmode ? d->tbsize / 2 : n % d->tbsize, pos = n / d->tbsize;
	fftwf_complex* dftc = ctx->dftc[0].get(), * dftc2 = ctx->dftc[1].get();
	float* dftr = ctx->dftr.get();
	for (int p = 0; p < d->np; p++)
	{
		auto& dstp = ctx->fb.back()->plane[p];
		if (!d->process[p])
			ctx->fb[k]->plane[p].copy_to(dstp);
		else
		{
			size_t tmp_h = d->padded_h[p], tmp_w = d->padded_w[p];
			size_t dst_h = dstp.h, dst_w = dstp.w, oy = (tmp_h - dst_h) / 2, ox = (tmp_w - dst_w) / 2;
			size_t stride_z = tmp_h * tmp_w, sbsize2 = d->sbsize * d->sbsize;
			float* buf = ctx->buf[p];
			if (!d->tmode || (ctx->pos != pos))
			{
				for (int i = 0; i < d->tbsize; i++)
					uf::copy_pad_reflect<T>(tmp_h, tmp_w, tmp_w * sizeof(T), ctx->tmp_p[i],
						dst_h, dst_w, ctx->fb[i]->stride(p), ctx->fb[i]->ptr(p), oy, ox);
				std::fill_n(buf, stride_z * (!d->tmode ? 1 : d->tbsize), 0.f);
				for (int y = 0; y <= tmp_h - d->sbsize; y += d->step)
				{
					for (int x = 0; x <= tmp_w - d->sbsize; x += d->step)
					{
						for (int z = 0; z < d->tbsize; z++)
							load(dftr + z * sbsize2, d->hw + z * sbsize2, d->sbsize, reinterpret_cast<T*>(ctx->tmp_p[z]) + y * tmp_w + x, tmp_w);
						fftwf_execute_dft_r2c(d->ft, dftr, dftc);
						if (d->zmean)
							remove_mean(d->ccnt, dftc, d->dftgc, dftc2);
						d->filter_coeffs(d->ccnt, dftc, d->sigmas[p], &d->f0beta);
						if (d->zmean)
							add_mean(d->ccnt, dftc, dftc2);
						fftwf_execute_dft_c2r(d->fti, dftc, dftr);
						if (!d->tmode)
						{
							if (d->smode)
								store(buf + y * tmp_w + x, tmp_w, dftr + k * sbsize2, d->hw + k * sbsize2, d->sbsize);
							else
								buf[y * tmp_w + x] = dftr[k * sbsize2 + oy * d->sbsize + ox] * d->hw[k * sbsize2 + oy * d->sbsize + ox];
						}
						else
						{
							if (d->smode)
							{
								for (int z = 0; z < d->tbsize; z++)
									store(buf + z * stride_z + y * tmp_w + x, tmp_w, dftr + z * sbsize2, d->hw + z * sbsize2, d->sbsize);
							}
							else
							{
								for (int z = 0; z < d->tbsize; z++)
								{
									buf[z * stride_z + y * tmp_w + x] =
										dftr[z * sbsize2 + oy * d->sbsize + ox] * d->hw[z * sbsize2 + oy * d->sbsize + ox];
								}
							}
						}
					}
				}
			}
			cast(reinterpret_cast<T*>(dstp.p), dst_h, dst_w, dstp.stride / sizeof(T),
				buf + k * stride_z * !!d->tmode + (oy * tmp_w + ox) * !!d->smode, tmp_w);
		}
	}
	ctx->pos = pos;
}

template void proc_temporal_avx<uint8_t>(const filter*, context*, int);
template void proc_temporal_avx<uint16_t>(const filter*, context*, int);
template void proc_temporal_avx<float>(const filter*, context*, int);

}
#endif