#include "config.h"
#if defined(BUILD_DFTTEST)
#ifndef VP_CXX

#include <fftw3.h>
#include "filter.h"

namespace dfttest {

struct context;

struct filter
{
	static const char name[], args[];
	VSVideoInfo vi;
	VSFilterDependency deps[1];
	std::mutex mutex;
	std::vector<std::unique_ptr<context>> cc;
	vvfcache_t srcc;
	filter(vmap&&);
	~filter();
	std::unique_ptr<context> ctx();
	int np, padded_h[3], padded_w[3], process[3];
	int smode, tmode, sbsize, sosize, tbsize, swin, twin;
	float sbeta, tbeta, f0beta;
	bool zmean;
	float* hw = 0, * sigmas[3] = {};
	fftwf_complex* dftgc = 0;
	fftwf_plan ft = 0, fti = 0;
	int bvolume, ccnt, step;
	void(*filter_coeffs)(size_t, fftwf_complex* __restrict, const float*, const void*);
	void(*proc)(const filter*, context*);
	void(*proc_temporal)(const filter*, context*, int);
};

struct context
{
	std::vector<vf_t> fb;
	std::vector<uint8_t*> tmp_p;
	std::vector<float*> buf;
	std::unique_ptr<uint8_t[]> ptr;
	uf::aligned_unique<float> dftr;
	uf::aligned_unique<fftwf_complex> dftc[2];
	int pos = -1;
	context(const filter*, vf_t&&);
};

template< typename T>
constexpr float scale_0 = std::is_same_v<T, float> ? 255.f : std::is_same_v<T, uint8_t> ? 1 : 1 / 256.f;
template< typename T>
constexpr float scale_1 = std::is_same_v<T, float> ? 1 / 255.f : std::is_same_v<T, uint8_t> ? 1 : 256.f;

void filter_0_avx(size_t, fftwf_complex* __restrict, const float*, const void*);
void filter_0_sqrt_avx(size_t, fftwf_complex* __restrict, const float*, const void*);
void filter_0_pow_avx(size_t, fftwf_complex* __restrict, const float*, const void*);
void filter_1_avx(size_t, fftwf_complex* __restrict, const float*, const void*);
template<typename T> void proc_avx(const filter*, context*);
template<typename T> void proc_temporal_avx(const filter*, context*, int);

}
#else

#ifndef USE_FFTW
#define USE_FFTW
#endif

#endif
#endif