#include "config.h"
#if defined(BUILD_ADDGRAIN)

#include <random>
#include "filter.h"

namespace addgrain {

struct drbg_ctr_aes
{
	using result_type = uint32_t;
	drbg_ctr_aes(uint64_t ctr = 0, uint64_t x = 0, const void* y = "0123456789abcdef")
	{
		seed(ctr, x);
		memcpy(k, y, 16);
	}
	void seed(uint64_t ctr, uint64_t x)
	{
		v[0] = ctr, v[1] = x, bi = -1;
	}
	uint32_t operator()()
	{
		if (bi >= 4)
		{
			const __m128i _k = _mm_loadu_si128((__m128i*)k), one = { 1 };
			__m128i _v = _mm_loadu_si128((__m128i*)v);
			_mm_storeu_si128((__m128i*)buf, _mm_aesenc_si128(_mm_aesenc_si128(_mm_aesenc_si128(_v, _k), _k), _k));
			_mm_storeu_si128((__m128i*)v, _mm_add_epi64(_v, one));
			bi = 0;
		}
		return buf[bi++];
	}
	NO_INLINE void fill_n(uint32_t* dst, size_t n)
	{
		const __m128i _k = _mm_loadu_si128((__m128i*)k), one = { 1 };
		__m128i _v = _mm_loadu_si128((__m128i*)v);
		for (uint32_t* p = dst, *end = p + n; p != end; p += 4)
		{
			_mm_storeu_si128((__m128i*)p, _mm_aesenc_si128(_mm_aesenc_si128(_mm_aesenc_si128(_v, _k), _k), _k));
			_v = _mm_add_epi64(_v, one);
		}
		_mm_storeu_si128((__m128i*)v, _v);
	}
private:
	uint32_t buf[4];
	uint64_t k[2], v[2];
	unsigned bi;
};

union float32_u
{
	unsigned i;
	float f;
	struct
	{
		unsigned m : 23;
		unsigned e : 8;
		unsigned s : 1;
	} s;
};

struct rsign
{
	float operator()(drbg_ctr_aes& g, float x)
	{
		float32_u u;
		if (i >= 32)
			t = g(), i = 0;
		u.f = x, u.i |= t << 31, t >>= 1, ++i;
		return u.f;
	}
private:
	unsigned i = -1, t = 0;
};

// [0,1] q=1.0,0.5,0.5 b=32
float uniform01_a(uint32_t r)
{
	constexpr float x = 1. / 4294967296ULL;
	return r * x;
}

// [0,1) q=1.0
float uniform01_24(uint32_t r)
{
	float32_u u1, u2;
	u1.i = 0b00111111100000000000000000000000 - (r & 1);
	u2.i = (r >> 9) | 0b00111111100000000000000000000000;
	return u2.f - u1.f;
}

// [-1,1) q=1.0 b=22
float uniform11_23(uint32_t r)
{
	float32_u u2;
	u2.i = 0b01000000000000000000000000000000 | (r >> 9);
	return u2.f - 3.f;
}

constexpr long double pi = 3.141592653589793238462643383279503L;
constexpr long double exp1 = 2.718281828459045235360287471352662L;

struct dist
{
	virtual ~dist() = default;
	virtual float operator()(drbg_ctr_aes& g) = 0;
};

struct uniform : dist
{
	float operator()(drbg_ctr_aes& g)
	{
		return uniform11_23(g());
	}
	static float mean() { return .5f; }
};

struct zignor : dist
{
	using result_type = float;
#define UNI(g) (uniform01_24(g()))
#define C 256
	zignor()
	{
		double d = r, e = std::exp(-.5 * d * d), _d = v / e;
		wn[0] = _d / scale;
		for (int i = 1; i < C; i++)
		{
			kn[i - 1] = d / _d * scale;
			fn[i - 1] = e;
			wn[i] = d / scale;
			_d = d, d = std::sqrt(-2. * std::log(v / d + e));
			e = std::exp(-.5 * d * d);
		}
		fn[C - 1] = 1;
		kn[C - 1] = 0;
	}
	float operator()(drbg_ctr_aes& g)
	{
		int hz = g();
		size_t iz = hz & (C - 1);
		if (std::abs(hz) < kn[iz])
			return hz * wn[iz];
		return nfix(g, hz, iz);
	}
	static float mean() { return std::sqrt(2 / pi); }
private:
	NO_INLINE float nfix(drbg_ctr_aes& g, int hz, size_t iz)
	{
		float x, y;
		for (;;)
		{
			if (iz == 0)
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
			iz = hz & (C - 1);
			if (std::abs(hz) < kn[iz])
				return hz * wn[iz];
		}
	}
	unsigned kn[C];
	float wn[C], fn[C];
#if C == 256
	static constexpr long double r = 3.654152885361008771645429720399516L;
	static constexpr long double v = 4.928673233974655347361775402336028e-03L;
#endif
#undef C
	static constexpr long double scale = 2147483648ULL;
#undef UNI
};

struct beta_dist : dist
{
	using result_type = float;
#define UNI(g) (uniform01_a(g()))
	beta_dist(float a = 1.f, float b = 1.f)
		: a(a), b(b)
	{
		use_gamma_dist = a >= 1 || b >= 1;
		if (use_gamma_dist)
		{
			p[0] = exp1 / (a + exp1);
			sqrt[0] = std::sqrt(2 * a - 1);
			p[1] = exp1 / (b + exp1);
			sqrt[1] = std::sqrt(2 * b - 1);
		}
		else
			p[0] = p[1] = sqrt[0] = sqrt[1] = 0;
	}
#ifdef FILTER_DEBUG
	~beta_dist()
	{
		fputc('\n', stderr);
	}
#endif
	float operator()(drbg_ctr_aes& g)
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
		return sg(g, x) / (x + y);
	}
	static float mean(float a, float b) { return a / (a + b); }
private:
	NO_INLINE float gamma_dist(drbg_ctr_aes& g, float a, float p, float sqrt)
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
		if (int count = a; count == a && count < 20)
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
	rsign sg;
	float p[2], sqrt[2], a, b;
	bool use_gamma_dist;
#undef UNI
};

template < typename T >
void add(void* dst, size_t h, size_t w, size_t pitch, drbg_ctr_aes& g, dist& d, float sigma)
{
	uint8_t* pdst = (uint8_t*)dst;
	if constexpr (std::is_floating_point_v<T>)
	{
		for (size_t i = 0; i < h; i++)
		{
			for (size_t j = 0; j < w; j++)
				reinterpret_cast<T*>(pdst + i * pitch)[j] += d(g) * sigma;
		}
	}
	else
	{
		for (size_t i = 0; i < h; i++)
		{
			for (size_t j = 0; j < w; j++)
				reinterpret_cast<T*>(pdst + i * pitch)[j] =
				std::clamp<int>(reinterpret_cast<T*>(pdst + i * pitch)[j] + (int)std::round(d(g) * sigma), 0, std::numeric_limits<T>::max());
		}
	}
}

enum { UNIFORM, NORMAL, BETA };

struct context;

struct filter
{
	VSVideoInfo vi;
	static const char name[], args[];
	filter(vmap&& map)
	{
		VSNodeRef* node = map.get<VSNodeRef*>("clip");
		vi = *vsapi->getVideoInfo(node);
		srcc = std::make_unique<vvfcache>(node);
		std::string d;
		map.get<std::string>(d, "d");
		char* y, * x = 0;
		float ss = 1;
		if (!d.empty())
		{
			x = strtok_r(d.data(), ":", &y);
			if (!strcmp(x, "normal"))
			{
				dt = NORMAL;
				ss /= zignor::mean();
			}
			else if (!strcmp(x, "beta"))
			{
				dt = BETA;
				if (sscanf(y, "%f,%f", &a, &b) != 2)
					a = 1, b = 1;
				ss /= beta_dist::mean(a, b);
			}
			else
				x = 0;
		}
		if (!x)
		{
			dt = UNIFORM;
			ss /= uniform::mean();
		}
		switch (srcc->css)
		{
		case 0:
			add = addgrain::add<uint8_t>;
			break;
		case 1:
			add = addgrain::add<uint16_t>;
			ss *= 257;
			break;
		case 2:
			add = addgrain::add<float>;
			ss /= 256;
			break;
		}
		if (map.get(sigma, "sigma"))
			sigma.push_back(1);
		sigma.resize(srcc->np, sigma.back());
		for (auto& sigma : sigma)
			sigma *= ss;
	}
	void proc(context* ctx, vf* dstf);
	std::unique_ptr<vvfcache> srcc;
	void(*add)(void*, size_t, size_t, size_t, drbg_ctr_aes&, dist&, float);
	int dt;
	float a, b;//move
	std::vector<float> sigma;
};

struct context
{
	drbg_ctr_aes g;
	std::unique_ptr<dist> d;
	uint64_t _seed = 0;
	context(filter* f)
	{
		switch (f->dt)
		{
		case UNIFORM:
			d = std::make_unique<uniform>();
			break;
		case NORMAL:
			d = std::make_unique<zignor>();
			break;
		case BETA:
			d = std::make_unique<beta_dist>(f->a, f->b);
			break;
		}
	}
};

void filter::proc(context* ctx, vf* dstf)
{
	for (int p = 0; p < srcc->np; p++)
	{
		if (sigma[p] > 0)
		{
			ctx->g.seed(0, ctx->_seed ^ p);
			add(dstf->p(p), dstf->h(p), dstf->w(p), dstf->pitch(p), ctx->g, *ctx->d, sigma[p]);
		}
	}
}

const VSFrameRef* get(int n, int activationReason, void** instanceData, void** frameData, VSFrameContext* frameCtx, VSCore* core, const VSAPI* vsapi)
{
	filter* d = static_cast<filter*>(*instanceData);
	if (activationReason == arInitial)
		d->srcc->request(n, frameCtx);
	if (activationReason != arAllFramesReady)
		return 0;
	thread_local context ctx(d);
	vf_t dstf, srcf;
	d->srcc->get(srcf, n, frameCtx);
	VSFrameRef* dst = uf::copy(dstf, srcf, core);
	ctx._seed = n;
	d->proc(&ctx, dstf.get());
	return dst;
}

void init(VSMap*, VSMap*, void** instanceData, VSNode* node, VSCore*, const VSAPI* vsapi)
{
	filter* d = static_cast<filter*>(*instanceData);
	vsapi->setVideoInfo(&d->vi, 1, node);
}

void create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi)
{
	if (!::vsapi)
		::vsapi = vsapi;
	try
	{
		filter* d = new filter(vmap(in));
		vsapi->createFilter(in, out, filter::name, init, get, uf::free<filter>, fmParallel, 0, d, core);
	}
	catch (const char* error)
	{
		vsapi->setError(out, error);
	}
}

const char filter::name[] = "addgrain";
const char filter::args[] = "clip:clip;sigma:float[]:opt;g:data:opt;d:data:opt;seed:int:opt;constant:int:opt;";

struct reg
{
	reg() {
		vregf.emplace_back([](VSRegisterFunction registerFunc, VSPlugin* plugin) {
			registerFunc(filter::name, filter::args, create, plugin, plugin); });
	}
} _;
}
#endif