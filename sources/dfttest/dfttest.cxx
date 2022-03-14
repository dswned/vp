/*
**                    dfttest for Avisynth+
**
**   2D/3D frequency domain denoiser.
**
**   Copyright (C) 2007-2010 Kevin Stone, 2017 (C) DJATOM
**             (C) 2020 pinterf
**
**   This program is free software; you can redistribute it and/or modify
**   it under the terms of the GNU General Public License as published by
**   the Free Software Foundation; either version 2 of the License, or
**   (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**   but WITHOUT ANY WARRANTY; without even the implied warranty of
**   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**   GNU General Public License for more details.
**
**   You should have received a copy of the GNU General Public License
**   along with this program; if not, write to the Free Software
**   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
// v1.9.6 adapted to vp

#include "config.h"
#if defined(BUILD_DFTTEST)

#if defined _HAS_CXX20 && _HAS_CXX20
#include <numbers>
#endif
#include <cmath>

#define ORG
/* break bitexact: overlaps fix, round 'away from zero' to 'nearest even', // to *, op reorder, etc. */

#ifdef ORG
#define VP_DISABLE_SATURATE_CAST_INTRIN
#endif

#include "cpu.h"
#include "dfttest.h"

namespace dfttest {

template< typename T>
void load(float* __restrict d, const float* w, size_t sbsize, const T* src, ptrdiff_t src_stride)
{
	for (size_t i = 0; i < sbsize; i++)
	{
		for (size_t j = 0; j < sbsize; j++)
			d[j] = src[j] * w[j] * scale_0<T>;
		d += sbsize;
		w += sbsize;
		src += src_stride;
	}
}

void store(float* __restrict dst, ptrdiff_t dst_stride, const float* s, const float* w, size_t sbsize)
{
	for (size_t i = 0; i < sbsize; i++)
	{
		for (size_t j = 0; j < sbsize; j++)
			dst[j] += s[j] * w[j];
		dst += dst_stride;
		s += sbsize;
		w += sbsize;
	}
}

template< typename T>
std::enable_if_t<std::is_integral_v<T>> cast(
	T* __restrict dst, size_t dst_h, size_t dst_w, ptrdiff_t dst_stride, const float* src, ptrdiff_t src_stride)
{
	for (size_t i = 0; i < dst_h; i++)
	{
		for (size_t j = 0; j < dst_w; j++)
			dst[j] = uf::saturate_cast<T>(src[j] * scale_1<T>);
		dst += dst_stride;
		src += src_stride;
	}
}

template< typename T>
std::enable_if_t<std::is_floating_point_v<T>> cast(
	T* __restrict dst, size_t dst_h, size_t dst_w, ptrdiff_t dst_stride, const float* src, ptrdiff_t src_stride)
{
	for (size_t i = 0; i < dst_h; i++)
	{
		for (size_t j = 0; j < dst_w; j++)
			dst[j] = src[j] * scale_1<T>;
		dst += dst_stride;
		src += src_stride;
	}
}

void remove_mean(size_t ccnt, fftwf_complex* __restrict dftc, const fftwf_complex* dftgc, fftwf_complex* __restrict dftc2)
{
	float gf = dftc[0][0] / dftgc[0][0];
	for (size_t i = 0; i < ccnt; i++)
	{
		dftc2[i][0] = gf * dftgc[i][0];
		dftc2[i][1] = gf * dftgc[i][1];
		dftc[i][0] -= dftc2[i][0];
		dftc[i][1] -= dftc2[i][1];
	}
}

void add_mean(size_t ccnt, fftwf_complex* __restrict dftc, const fftwf_complex* dftc2)
{
	for (size_t i = 0; i < ccnt; i++)
	{
		dftc[i][0] += dftc2[i][0];
		dftc[i][1] += dftc2[i][1];
	}
}

void filter_0(size_t ccnt, fftwf_complex* __restrict dftc, const float* sigmas, const void* /*beta = 1.0*/)
{
	for (size_t i = 0; i < ccnt; i++)
	{
		float psd = dftc[i][0] * dftc[i][0] + dftc[i][1] * dftc[i][1];
		float mult = std::max(psd - sigmas[i], 0.f) / (psd + 1e-15f);
		dftc[i][0] *= mult;
		dftc[i][1] *= mult;
	}
}

void filter_0_sqrt(size_t ccnt, fftwf_complex* __restrict dftc, const float* sigmas, const void* /*beta = 0.5*/)
{
	for (size_t i = 0; i < ccnt; i++)
	{
		float psd = dftc[i][0] * dftc[i][0] + dftc[i][1] * dftc[i][1];
		float mult = std::sqrt(std::max(psd - sigmas[i], 0.f) / (psd + 1e-15f));
		dftc[i][0] *= mult;
		dftc[i][1] *= mult;
	}
}

void filter_0_pow(size_t ccnt, fftwf_complex* __restrict dftc, const float* sigmas, const void* beta)
{
	for (size_t i = 0; i < ccnt; i++)
	{
		float psd = dftc[i][0] * dftc[i][0] + dftc[i][1] * dftc[i][1];
		float mult = std::pow(std::max(psd - sigmas[i], 0.f) / (psd + 1e-15f), *static_cast<const float*>(beta));
		dftc[i][0] *= mult;
		dftc[i][1] *= mult;
	}
}

void filter_1(size_t ccnt, fftwf_complex* __restrict dftc, const float* sigmas, const void*)
{
	for (size_t i = 0; i < ccnt; i++)
	{
		float psd = dftc[i][0] * dftc[i][0] + dftc[i][1] * dftc[i][1];
		if (psd < sigmas[i])
			dftc[i][0] = dftc[i][1] = 0.f;
	}
}

void filter_2(size_t ccnt, fftwf_complex* __restrict dftc, const float* sigmas, const void*)
{
	for (size_t i = 0; i < ccnt; i++)
	{
		dftc[i][0] *= sigmas[i];
		dftc[i][1] *= sigmas[i];
	}
}

template< typename T>
void proc(const filter* d, context* ctx)
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

template< typename T>
void proc_temporal(const filter* d, context* ctx, int n)
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
			size_t stride_z = tmp_h * tmp_w, sbsize2 = (size_t)d->sbsize * d->sbsize;
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

void select_functions(filter* d, int ftype, int opt)
{
	if (ftype == 0)
	{
		if (d->f0beta == 1.f)
			d->filter_coeffs = filter_0;
		else if (d->f0beta == .5f)
			d->filter_coeffs = filter_0_sqrt;
		else
			d->filter_coeffs = filter_0_pow;
	}
	else if (ftype == 1)
		d->filter_coeffs = filter_1;
	else if (ftype == 2)
		d->filter_coeffs = filter_2;
	switch (d->srcc->css)
	{
	case 0:
		d->proc = proc<uint8_t>;
		d->proc_temporal = proc_temporal<uint8_t>;
		break;
	case 1:
		d->proc = proc<uint16_t>;
		d->proc_temporal = proc_temporal<uint16_t>;
		break;
	case 2:
		d->proc = proc<float>;
		d->proc_temporal = proc_temporal<float>;
		break;
	}
#if 1
	uf::cpu cpu = opt < 0 ? uf::query_x86_capabilities() : uf::cpu(opt);
	if (cpu.flags & uf::cpuflags::avx)
	{
		if (ftype == 0)
		{
			if (d->f0beta == 1.f)
				d->filter_coeffs = filter_0_avx;
			else if (d->f0beta == .5f)
				d->filter_coeffs = filter_0_sqrt_avx;
			else
				d->filter_coeffs = filter_0_pow_avx;
		}
		else if (ftype == 1)
			d->filter_coeffs = filter_1_avx;
		switch (d->srcc->css)
		{
		case 0:
			d->proc = proc_avx<uint8_t>;
			d->proc_temporal = proc_temporal_avx<uint8_t>;
			break;
		case 1:
			d->proc = proc_avx<uint16_t>;
			d->proc_temporal = proc_temporal_avx<uint16_t>;
			break;
		case 2:
			d->proc = proc_avx<float>;
			d->proc_temporal = proc_temporal_avx<float>;
			break;
		}
	}
#endif
}

double besselI0(double p)
{
	p /= 2.;
	double n = 1., t = 1., d = 1.;
	int k = 1;
	double v;
	do
	{
		n *= p;
		d *= k;
		v = n / d;
		t += v * v;
	} while (++k < 15 && v > 1e-8);
	return t;
}

double get_win_value(double n, double size, int win, double beta)
{
	constexpr double pi = 3.141592653589793;
	switch (win)
	{
	case 0: // hanning
		return 0.5 - 0.5 * std::cos(2. * pi * n / size);
	case 1: // hamming
		return 0.53836 - 0.46164 * std::cos(2. * pi * n / size);
	case 2: // blackman
		return 0.42 - 0.5 * std::cos(2. * pi * n / size) + 0.08 * std::cos(4. * pi * n / size);
	case 3: // 4 term blackman-harris
		return 0.35875 - 0.48829 * std::cos(2. * pi * n / size) + 0.14128 * std::cos(4. * pi * n / size) - 0.01168 * std::cos(6. * pi * n / size);
	case 4: // kaiser-bessel
	{
		const double v = 2. * n / size - 1.;
		return besselI0(pi * beta * std::sqrt(1. - v * v)) / besselI0(pi * beta);
	}
	case 5: // 7 term blackman-harris
		return 0.27105140069342415 -
			0.433297939234486060 * std::cos(2. * pi * n / size) +
			0.218122999543110620 * std::cos(4. * pi * n / size) -
			0.065925446388030898 * std::cos(6. * pi * n / size) +
			0.010811742098372268 * std::cos(8. * pi * n / size) -
			7.7658482522509342E-4 * std::cos(10. * pi * n / size) +
			1.3887217350903198E-5 * std::cos(12. * pi * n / size);
	case 6: // flat top
		return 0.2810639 - 0.5208972 * std::cos(2. * pi * n / size) + 0.1980399 * std::cos(4. * pi * n / size);
	case 7: // rectangular
		return 1.;
	case 8: // Bartlett
		return 2. / size * (size / 2. - std::abs(n - size / 2.));
	case 9: // Bartlett-Hann
		return 0.62 - 0.48 * (n / size - 0.5) - 0.38 * std::cos(2. * pi * n / size);
	case 10: // Nuttall
		return 0.355768 - 0.487396 * std::cos(2. * pi * n / size) + 0.144232 * std::cos(4. * pi * n / size) - 0.012604 * std::cos(6. * pi * n / size);
	case 11: // Blackman-Nuttall
		return 0.3635819 - 0.4891775 * std::cos(2. * pi * n / size) + 0.1365995 * std::cos(4. * pi * n / size) - 0.0106411 * std::cos(6. * pi * n / size);
	default:
		return 0.;
	}
}

void normalize_for_overlap_add(double* hw, int bsize, int osize)
{
	std::vector<double> nw(bsize);
	for (int i = 0; i < bsize; i++)
	{
		for (int step = bsize - osize, j = i % step; j < bsize; j += step)
			nw[i] += hw[j] * hw[j];
	}
	for (int i = 0; i < bsize; i++)
		hw[i] /= std::sqrt(nw[i]);
}

void make_window(filter* f, float* hw)
{
	std::vector<double> tw(f->tbsize), sw(f->sbsize);
	for (int i = 0; i < f->tbsize; i++)
		tw[i] = get_win_value(i + .5, f->tbsize, f->twin, f->tbeta);
	for (int i = 0; i < f->sbsize; i++)
		sw[i] = get_win_value(i + .5, f->sbsize, f->swin, f->sbeta);
	if (f->smode)
		normalize_for_overlap_add(sw.data(), f->sbsize, f->sosize);
	double nscale = 1 / std::sqrt(f->bvolume);
	for (int z = 0; z < f->tbsize; z++)
		for (int y = 0; y < f->sbsize; y++)
			for (int x = 0; x < f->sbsize; x++)
				hw[(z * f->sbsize + y) * f->sbsize + x] = static_cast<float>(tw[z] * sw[y] * sw[x] * nscale);
}

context::context(const filter* f, vf_t&& dstf)
{
	fb.resize((size_t)f->tbsize + 1);
	fb.back() = std::move(dstf);
	uint8_t* first = 0, * second = 0;
	for (int i = 0; i < f->tbsize; i++)
		tmp_p.push_back(second), second += (size_t)f->padded_h[0] * f->padded_w[0] << f->srcc->css;
	if (!f->tmode)
	{
		buf.resize(f->np, reinterpret_cast<float*>(second));
		second += sizeof(float) * f->padded_h[0] * f->padded_w[0];
	}
	else
	{
		for (int p = 0; p < f->np; p++)
		{
			buf.push_back(reinterpret_cast<float*>(second));
			second += sizeof(float) * f->padded_h[p] * f->padded_w[p] * f->tbsize;
		}
	}
	ptrdiff_t n_alloc = std::distance(first, second);
	ptr = std::make_unique<uint8_t[]>(n_alloc);
	if (!ptr)
		throw uf::format("dfttest::context(): failed to allocate %zu bytes", n_alloc);
	for (auto& p : tmp_p)
		p = ptr.get() + std::distance(first, p);
	for (auto& p : buf)
		p = reinterpret_cast<float*>(ptr.get() + std::distance(first, reinterpret_cast<uint8_t*>(p)));
	dftr = uf::make_aligned_unique<float>((size_t)f->bvolume + 64, 64);
	dftc[0] = uf::make_aligned_unique<fftwf_complex>((size_t)f->ccnt + 32, 64);
	dftc[1] = uf::make_aligned_unique<fftwf_complex>((size_t)f->ccnt + 32, 64);
}

filter::filter(vmap&& map) try
{
	VSNode* node = map.get<VSNode*>("clip");
	vi = *vsapi->getVideoInfo(node);
	srcc = std::make_unique<vvfcache>(node);
	np = srcc->np;
	int ftype = map.get("ftype", 0);
	smode = map.get("smode", 1);
	sbsize = map.get("sbsize", 16);
	sosize = smode ? map.get("sosize", 12) : sbsize / 2;
	tmode = map.get("tmode", 0);
	tbsize = map.get("tbsize", 3);
	swin = smode && !sosize ? 7 : map.get("swin", 0);
	twin = tmode ? 7 : map.get("twin", 7);
	sbeta = map.get("sbeta", 2.5f);
	tbeta = map.get("tbeta", 2.5f);
	zmean = map.get("zmean", true);
	f0beta = map.get("f0beta", 1.f);
	std::vector<float> sigma;
	if (!map.get(sigma, "sigma"))
		sigma.resize(np, 8.f);
	else
		sigma.resize(np, sigma.back());
	int opt = map.get("opt", -1);
	bool format_check = 1;
	switch (srcc->st)
	{
	case sample_type_e::BYTE:
	case sample_type_e::WORD:
		if (srcc->bps != 8 && srcc->bps != 16)
			format_check = 0;
		break;
	case sample_type_e::FLOAT:
		if (srcc->bps != 32)
			format_check = 0;
		zmean = true;
		break;
	default:
		format_check = 0;
	}
	if (!format_check)
		throw "only 8 or 16 bit integer and 32 bit float input supported"s;
	if (ftype < 0 || ftype > 2)
		throw "ftype must be set to 0, 1, or 2!"s;
	if (twin < 0 || twin > 11)
		throw "twin must be between 0 and 11 (inclusive)!"s;
	if (swin < 0 || swin > 11)
		throw "swin must be between 0 and 11 (inclusive)!"s;
	if (tbsize < 1 || tbsize > srcc->nf)
		throw "tbsize must be between 1 and number of frames in the clip (inclusive)"s;
	if (sbsize < 1)
		throw "sbsize must be greater than or equal to 1"s;
	if (sosize < 0 || sosize >= sbsize)
		throw "sosize must be between 0 and sbsize-1 (inclusive)"s;
	if (!tmode && ~tbsize & 1)
		throw "tbsize must be odd when using tmode=0!"s;
	if (!smode && ~sbsize & 1)
		throw "sbsize must be odd when using smode=0!"s;
	if (sosize > (sbsize - sosize) && sosize % (sbsize - sosize))
		throw "spatial overlap greater than 50% requires that (sbsize-sosize) is a divisor of sbsize!"s;
	std::vector<float> ssx, ssy, ssz;
	auto ss_resize = [](std::vector<float>& x, size_t size)
	{
		size_t n = size >> 1, a = n + (size & 1), b = n + 1;
		if (x.size() != b)
			return -1;
		x.resize(size);
		std::reverse_copy(x.data() + 1, x.data() + a, x.data() + b);
		return 0;
	};
	if (map.get(ssx, "ssx") && ssx.size() != (size_t)sbsize / 2 + 1)
		throw "ssx.size() != sbsize/2+1"s;
	if (map.get(ssy, "ssy") && ss_resize(ssy, sbsize))
		throw "ssy.size() != sbsize/2+1"s;
	if (map.get(ssz, "ssz") && ss_resize(ssz, tbsize))
		throw "ssz.size() != tbsize/2+1"s;
	select_functions(this, ftype, opt);
	bvolume = sbsize * sbsize * tbsize;
	ccnt = (sbsize / 2 + 1) * sbsize * tbsize;
	step = smode ? sbsize - sosize : 1;
	hw = uf::aligned_alloc<float>((size_t)bvolume + 64, 64);
	make_window(this, hw);
	dftgc = uf::aligned_alloc<fftwf_complex>((size_t)ccnt + 32, 64);
	float* dftgr = uf::aligned_alloc<float>(bvolume, 64);
	if (tbsize > 1)
	{
		ft = fftwf_plan_dft_r2c_3d(tbsize, sbsize, sbsize, dftgr, dftgc, FFTW_MEASURE | FFTW_DESTROY_INPUT);
		fti = fftwf_plan_dft_c2r_3d(tbsize, sbsize, sbsize, dftgc, dftgr, FFTW_MEASURE | FFTW_DESTROY_INPUT);
	}
	else
	{
		ft = fftwf_plan_dft_r2c_2d(sbsize, sbsize, dftgr, dftgc, FFTW_MEASURE | FFTW_DESTROY_INPUT);
		fti = fftwf_plan_dft_c2r_2d(sbsize, sbsize, dftgc, dftgr, FFTW_MEASURE | FFTW_DESTROY_INPUT);
	}
	float wscale = 0.f;
	for (int i = 0; i < bvolume; i++)
	{
		dftgr[i] = 255.f * hw[i];
		wscale += hw[i] * hw[i];
	}
	fftwf_execute_dft_r2c(ft, dftgr, dftgc);
	uf::aligned_free(dftgr);
	for (int p = 0; p < np; p++)
	{
		process[p] = sigma[p] > 0.f;
#ifndef ORG
		padded_h[p] = uf::ceil_n(srcc->h(p) + sosize, step) + sosize;
		padded_w[p] = uf::ceil_n(srcc->w(p) + sosize, step) + sosize;
		sigma[p] *= wscale;
#else
#define EXTRA(a, b) (((a) % (b)) ? ((b) - ((a) % (b))) : 0)
		const int width = srcc->w(p);
		const int height = srcc->h(p);
		if (smode == 0)
		{
			const int ae = (sbsize >> 1) << 1;
			padded_w[p] = width + ae;
			padded_h[p] = height + ae;
		}
		else
		{
			const int ae = std::max(sbsize - sosize, sosize) * 2;
			padded_w[p] = width + EXTRA(width, sbsize) + ae;
			padded_h[p] = height + EXTRA(height, sbsize) + ae;
		}
		sigma[p] /= 1 / wscale;
#undef EXTRA
#endif
		sigmas[p] = uf::aligned_alloc<float>((size_t)ccnt + 64, 64);
	}
	if (!ssx.empty() || !ssy.empty() || !ssz.empty())
	{
		int sz = tbsize, sy = sbsize, sx = sbsize / 2 + 1;
		if (ssx.empty())
			ssx.resize(sx, 1);
		if (ssy.empty())
			ssy.resize(sy, 1);
		if (ssz.empty())
			ssz.resize(sz, 1);
		for (int p = 0; p < np; p++)
			for (int z = 0; z < sz; z++)
				for (int y = 0; y < sy; y++)
					for (int x = 0; x < sx; x++)
						sigmas[p][(z * sy + y) * sx + x] = ssx[x] * ssy[y] * ssz[z] * sigma[p];
	}
	else
	{
		for (int p = 0; p < np; p++)
			std::fill_n(sigmas[p], ccnt, sigma[p]);
	}
	deps[0] = { node, rpGeneral };
}
catch (...)
{
	// possible 'that' be destroyed?
	for (int p = 0; p < np && sigmas[p]; p++)
		uf::aligned_free(sigmas[p]);
	if (hw)
		uf::aligned_free(hw);
	if (dftgc)
		uf::aligned_free(dftgc);
	if (ft)
		fftwf_destroy_plan(ft);
	if (fti)
		fftwf_destroy_plan(fti);
}

filter::~filter()
{
	for (int p = 0; p < np; p++)
		uf::aligned_free(sigmas[p]);
	uf::aligned_free(hw);
	uf::aligned_free(dftgc);
	fftwf_destroy_plan(ft);
	fftwf_destroy_plan(fti);
}

std::unique_ptr<context> filter::ctx()
{
	return std::make_unique<context>(this, std::make_shared<vvf>());
}

const VSFrame* get(int n, int activationReason,
	void* instanceData, void**, VSFrameContext* frameCtx, VSCore* core, const VSAPI* vsapi)
{
	filter* d = static_cast<filter*>(instanceData);
	if (activationReason == arInitial)
		d->srcc->request(n, frameCtx);
	if (activationReason != arAllFramesReady)
		return 0;
	std::unique_ptr<context> ctx;
	std::unique_lock lock(d->mutex);
	if (!d->cc.empty())
		ctx = std::move(d->cc.back()), d->cc.pop_back();
	lock.unlock();
	const VSFrame* dst = 0;
	try
	{
		if (!ctx)
			ctx = d->ctx();
		vvf* dstf = static_cast<vvf*>(ctx->fb.back().get());
		d->srcc->get(ctx->fb.front(), n, frameCtx);
		dstf->create(*ctx->fb.front(), core);
		d->proc(d, ctx.get());
		dst = dstf->release();
	}
	catch (const std::string& e)
	{
		vsapi->setFilterError(e.c_str(), frameCtx);
	}
	if (ctx)
	{
		lock.lock();
		d->cc.emplace_back(std::move(ctx));
		lock.unlock();
	}
	return dst;
}

const VSFrame* get_temporal(int n, int activationReason,
	void* instanceData, void**, VSFrameContext* frameCtx, VSCore* core, const VSAPI* vsapi)
{
	filter* d = static_cast<filter*>(instanceData);
	if (activationReason == arInitial)
	{
		int k = !d->tmode ? d->tbsize / 2 : n % d->tbsize, first = n - k;
		for (int i = std::max(first, 0), last = std::min(first + d->tbsize, d->srcc->nf); i < last; i++)
			d->srcc->request(i, frameCtx);
	}
	if (activationReason != arAllFramesReady)
		return 0;
	std::unique_ptr<context> ctx;
	std::unique_lock lock(d->mutex);
	if (!d->cc.empty())
		ctx = std::move(d->cc.back()), d->cc.pop_back();
	lock.unlock();
	const VSFrame* dst = 0;
	try
	{
		if (!ctx)
			ctx = d->ctx();
		vvf* dstf = static_cast<vvf*>(ctx->fb.back().get());
		int k = !d->tmode ? d->tbsize / 2 : n % d->tbsize, first = n - k;
		for (int i = first, last = first + d->tbsize; i < last; i++)
			d->srcc->get(ctx->fb[i - first], std::clamp(i, 0, d->srcc->nf - 1), frameCtx);
		dstf->create(*ctx->fb[k], core);
		d->proc_temporal(d, ctx.get(), n);
		dst = dstf->release();
	}
	catch (const std::string& e)
	{
		vsapi->setFilterError(e.c_str(), frameCtx);
	}
	if (ctx)
	{
		lock.lock();
		d->cc.emplace_back(std::move(ctx));
		lock.unlock();
	}
	return dst;
}

void create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi)
{
	if (!::vsapi)
		::vsapi = vsapi;
	try
	{
		filter* d = new filter(in);
		return vsapi->createVideoFilter(out, filter::name, &d->vi,
			d->tbsize == 1 ? get : get_temporal, uf::free<filter>, !d->tmode ? fmParallel : fmParallelRequests,
			d->deps, std::size(d->deps), d, core);
	}
	catch (const std::string& e)
	{
		vsapi->mapSetError(out, e.c_str());
	}
}

const char filter::name[] = "dfttest";
const char filter::args[] =
"clip:vnode;"
"sigma:float[]:opt;"
"tmode:int:opt;"
"tbsize:int:opt;"
"smode:int:opt;"
"sbsize:int:opt;"
"sosize:int:opt;"
"zmean:int:opt;"
"ftype:int:opt;"
"f0beta:float:opt;"
"twin:int:opt;"
"swin:int:opt;"
"tbeta:float:opt;"
"sbeta:float:opt;"
"ssx:float[]:opt;"
"ssy:float[]:opt;"
"ssz:float[]:opt;"
"opt:int:opt;";

void reg_f(VSPlugin* plugin, const VSPLUGINAPI* vspapi)
{
	vspapi->registerFunction(filter::name, filter::args, "clip:vnode;", create, 0, plugin);
}

PUSH_REG_F(reg_f);
}
#endif