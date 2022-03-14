#include "config.h"
#if defined(BUILD_ADDGRAIN)

#include <random>

#include "filter.h"
#include "cpu.h"
#include "addgrain.h"

namespace addgrain {

// [0,1]
float uniform01(uint32_t u)
{
	constexpr float x = 1.f / 0x100000000;
	return u * x;
}

// [0,1)
float uniform01_24(uint32_t u)
{
	uf::sf32 u1, u2;
	u1.i = 0b00111111100000000000000000000000 - (u & 1);
	u2.i = (u >> 9) | 0b00111111100000000000000000000000;
	return u2.f - u1.f;
}

// [0,1)
float uniform01_23(uint32_t u)
{
	uf::sf32 u2;
	u2.i = 0b00111111100000000000000000000000 | (u >> 9);
	return u2.f - 1.f;
}

// [-1,1)
float uniform11_23(uint32_t u)
{
	uf::sf32 u2;
	u2.i = 0b01000000000000000000000000000000 | (u >> 9);
	return u2.f - 3.f;
}

__m128 uniform11_23(__m128i u)
{
	u = _mm_srli_epi32(u, 9);
	__m128 f = _mm_or_ps(_mm_castsi128_ps(u), _mm_set1_ps(2.f));
	f = _mm_sub_ps(f, _mm_set1_ps(3.f));
	return f;
}

constexpr double pi = 3.141592653589793238462643383279503;
constexpr double exp1 = 2.718281828459045235360287471352662;

struct uniform : dist
{
	float operator()(drbg& g)
	{
		return uniform11_23(g());
	}
	A_AES __m128 get(drbg& g)
	{
		return uniform11_23(g.get<__m128i>());
	}
	static float mean() { return .5f; }
};

struct zignor : dist
{
	zignor()
	{
		double d = r, e = std::exp(-.5 * d * d), _d = v / e;
		wn[0] = (float)(_d / scale);
		for (int i = 1; i < 256; i++)
		{
			kn[i - 1] = (unsigned)(d / _d * scale);
			fn[i - 1] = (float)e;
			wn[i] = (float)(d / scale);
			_d = d, d = std::sqrt(-2. * std::log(v / d + e));
			e = std::exp(-.5 * d * d);
		}
		fn[255] = 1;
		kn[255] = 0;
	}
	float operator()(drbg& g)
	{
		int hz = g();
		unsigned iz = hz & 255;
		if ((unsigned)std::abs(hz) < kn[iz])
			return hz * wn[iz];
		return nfix(g, hz, iz);
	}
	static float mean() { return (float)std::sqrt(2. / pi); }
private:
	NO_INLINE float nfix(drbg& g, int hz, unsigned iz)
	{
#define UNI(g) (uniform01(g()))
		float x, y;
		for (;;)
		{
			if (!iz)
			{
				do
				{
					x = std::log(UNI(g)) * (1 / float(r));
					y = -std::log(UNI(g));
				} while (y + y < x * x);
				return (hz > 0) ? float(r) - x : x - float(r);
			}
			x = hz * wn[iz];
			if (fn[iz] + UNI(g) * (fn[iz - 1] - fn[iz]) < std::exp(-.5f * x * x))
				return x;
			hz = g();
			iz = hz & 255;
			if ((unsigned)std::abs(hz) < kn[iz])
				return hz * wn[iz];
		}
#undef UNI
	}
	unsigned kn[256];
	float wn[256], fn[256];
	static constexpr double r = 3.654152885361008771645429720399516;
	static constexpr double v = 4.928673233974655347361775402336028e-03;
	static constexpr double scale = 2147483648.;
};

struct beta : dist
{
	using result_type = float;
#define UNI(g) (uniform01(g()))
	beta(float a = 1.f, float b = 1.f)
		: a(a), b(b)
	{
		use_gamma_dist = a >= 1 || b >= 1;
		if (use_gamma_dist)
		{
			p[0] = (float)(exp1 / (a + exp1));
			sqrt[0] = std::sqrt(2 * a - 1);
			p[1] = (float)(exp1 / (b + exp1));
			sqrt[1] = std::sqrt(2 * b - 1);
		}
		else
			p[0] = p[1] = sqrt[0] = sqrt[1] = 0;
	}
	float operator()(drbg& g)
	{
		float x, y;
		if (use_gamma_dist)
		{
			x = gamma_dist(g, a, p[0], sqrt[0]);
			y = gamma_dist(g, b, p[1], sqrt[1]);
		}
		else
		{
			do
			{
				x = std::pow(UNI(g), 1 / a);
				y = std::pow(UNI(g), 1 / b);
			} while (x + y > 1);
		}
		return uf::sf32(x).bit_or(g() & 0x80000000) / (x + y);
	}
	static float mean(float a, float b) { return a / (a + b); }
private:
	NO_INLINE float gamma_dist(drbg& g, float a, float p, float sqrt)
	{
		constexpr float b = 1;
		float x, y, u, v, q;
		if (a < 1)
		{
			for (;;)
			{
				u = UNI(g);
				do
				{
					v = UNI(g);
				} while (v == 0);
				if (u < p)
				{
					x = std::pow(v, 1 / a);
					q = std::exp(-x);
				}
				else
				{
					x = 1 - std::log(v);
					q = std::pow(x, a - 1);
				}
				if (UNI(g) < q)
					return b * x;
			}
		}
		if (a == 1)
			return b * -std::log(1 - UNI(g));
		if (int count = (int)a; count == a && count < 20)
		{
			y = UNI(g);
			while (--count)
			{
				do
				{
					u = UNI(g);
				} while (u == 0);
				y *= u;
			}
			return b * -std::log(y);
		}
		for (;;)
		{
			y = std::tan(float(pi) * UNI(g));
			x = sqrt * y + a - 1;
			if (x > 0 && UNI(g) <= (1 + y * y) * std::exp((a - 1) * std::log(x / (a - 1)) - sqrt * y))
				return b * x;
		}
	}
	float p[2], sqrt[2], a, b;
	bool use_gamma_dist;
#undef UNI
};

template< typename T>
void add(size_t h, size_t w, ptrdiff_t stride1, void* dst_p, drbg& g, dist& d, float sigma)
{
	T* dst = static_cast<T*>(dst_p);
	ptrdiff_t stride = stride1 / sizeof(T);
	if constexpr (std::is_floating_point_v<T>)
	{
		for (int i = 0; i < h; i++)
		{
			for (int j = 0; j < w; j++)
				dst[i * stride + j] += d(g) * sigma;
		}
	}
	else
	{
		for (int i = 0; i < h; i++)
		{
			for (int j = 0; j < w; j++)
				dst[i * stride + j] = uf::saturate_cast<T>(dst[i * stride + j] + d(g) * sigma);
		}
	}
}

template< typename T, std::enable_if_t<std::is_same_v<T, float>, int> = 0>
A_AVX void add_avx(size_t h, size_t w, ptrdiff_t stride1, void* dst_p, drbg& g, dist& d, float sigma)
{
	T* dst = static_cast<T*>(dst_p);
	ptrdiff_t stride = stride1 / sizeof(T);
	const __m128 k_sigma = _mm_broadcast_ss(&sigma);
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j += 4)
			_mm_storeu_ps(dst + j, _mm_add_ps(_mm_mul_ps(d.get(g), k_sigma), _mm_loadu_ps(dst + j)));
		dst += stride;
	}
}

template< typename T, std::enable_if_t<std::is_same_v<T, uint16_t>, int> = 0>
A_AVX void add_avx(size_t h, size_t w, ptrdiff_t stride1, void* dst_p, drbg& g, dist& d, float sigma)
{
	T* dst = static_cast<T*>(dst_p);
	ptrdiff_t stride = stride1 / sizeof(T);
	const __m128 k_sigma = _mm_broadcast_ss(&sigma);
	const __m128i k_zero = _mm_setzero_si128();
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j += 4)
		{
			__m128i x = _mm_loadl_epi64((const __m128i*)(dst + j));
			__m128 lo = _mm_cvtepi32_ps(_mm_unpacklo_epi16(x, k_zero));
			x = _mm_cvtps_epi32(_mm_add_ps(_mm_mul_ps(d.get(g), k_sigma), lo));
			_mm_storel_epi64((__m128i*)(dst + j), _mm_packus_epi32(x, x));
		}
		dst += stride;
	}
}

enum { UNIFORM, NORMAL, BETA };

struct context
{
	drbg g;
	int n;
};

impl::impl(int css, int dt, float a, float b)
{
	float scale = 1;
	switch (css)
	{
	case 0:
		add = addgrain::add<uint8_t>;
		break;
	case 1:
		add = addgrain::add<uint16_t>;
		scale = 256.f;
		break;
	case 2:
		add = addgrain::add<float>;
		scale = 1 / 256.f;
		break;
	default:
		ASSUME(0);
	}
	switch (dt)
	{
	case UNIFORM:
		d = std::make_unique<uniform>();
		sigma = scale / uniform::mean();
		break;
	case NORMAL:
		d = std::make_unique<zignor>();
		sigma = scale / zignor::mean();
		break;
	case BETA:
		d = std::make_unique<beta>(a, b);
		sigma = scale / beta::mean(a, b);
		break;
	default:
		ASSUME(0);
	}
	uf::cpu cpu = uf::query_x86_capabilities();
	if (cpu.flags & uf::cpuflags::avx)
	{
		if (dt == UNIFORM && css == 1)
			add = addgrain::add_avx<uint16_t>;
		if (dt == UNIFORM && css == 2)
			add = addgrain::add_avx<float>;
	}
}

struct filter
{
	static const char name[], args[];
	VSVideoInfo vi;
	VSFilterDependency deps[1];
	std::unique_ptr<impl> addgrain;
	std::vector<float> sigma;
	std::unique_ptr<vvfcache> srcc;
	filter(vmap&& map)
	{
		VSNode* node = map.get<VSNode*>("clip");
		vi = *vsapi->getVideoInfo(node);
		srcc = std::make_unique<vvfcache>(node);
		int dt;
		float a = 0, b = 0;
		if (std::string d; map.get(d, "d"))
		{
			char* y = 0, * x = strtok_r(d.data(), ":", &y);
			if (!strcmp(x, "normal"))
				dt = NORMAL;
			else if (!strcmp(x, "beta"))
			{
				dt = BETA;
				if (!y || sscanf(y, "%f,%f", &a, &b) != 2)
					a = 1, b = 1;
			}
			else if (!strcmp(x, "uniform"))
				dt = UNIFORM;
			else
				throw ""s;
		}
		if (!map.get(sigma, "sigma"))
			sigma.push_back(1);
		sigma.resize(srcc->np, sigma.back());
		addgrain = std::make_unique<impl>(srcc->css, dt, a, b);
		deps[0] = { node, rpStrictSpatial };
	}
	void proc(context* ctx, vf* dstf);
};

void filter::proc(context* ctx, vf* dstf)
{
	for (int p = 0; p < srcc->np; p++)
	{
		if (sigma[p] > 0)
		{
			ctx->g.seed((uint64_t)p << 32 | ctx->n);
			addgrain->add(dstf->h(p), dstf->w(p), dstf->stride(p), dstf->ptr(p), ctx->g, *addgrain->d, sigma[p] * addgrain->sigma);
		}
	}
}

const VSFrame* get(int n, int activationReason,
	void* instanceData, void**, VSFrameContext* frameCtx, VSCore* core, const VSAPI*)
{
	filter* d = static_cast<filter*>(instanceData);
	if (activationReason == arInitial)
		d->srcc->request(n, frameCtx);
	if (activationReason != arAllFramesReady)
		return 0;
	thread_local context ctx;
	vvf dstf, srcf;
	d->srcc->get(srcf, n, frameCtx);
	dstf.create_copy(srcf, core);
	ctx.n = n;
	d->proc(&ctx, &dstf);
	return dstf.release();
}

void create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi)
{
	if (!::vsapi)
		::vsapi = vsapi;
	try
	{
		filter* d = new filter(in);
		vsapi->createVideoFilter(out, filter::name, d->srcc->vi, get, uf::free<filter>, fmParallel, d->deps, std::size(d->deps), d, core);
	}
	catch (const std::string& e)
	{
		vsapi->mapSetError(out, e.c_str());
	}
}

const char filter::name[] = "addgrain";
const char filter::args[] =
"clip:vnode;"
"sigma:float[]:opt;"
"d:data:opt;"
"seed:int:opt;"
"constant:int:opt;";

void reg_f(VSPlugin* plugin, const VSPLUGINAPI* vspapi)
{
	vspapi->registerFunction(filter::name, filter::args, "clip:vnode;", create, 0, plugin);
}

PUSH_REG_F(reg_f);
}
#endif