#include "config.h"
#if defined(BUILD_ADDGRAIN)

#include "filter.h"

namespace addgrain {

struct drbg
{
	drbg(const std::string& seed_str = "0123456789abcdef")
	{
		std::copy_n(seed_str.data(), 16, m_k);
	}
	void seed(uint64_t x)
	{
		m_v[0] = 0, m_v[1] = x, m_i = 0;
	}
	uint32_t operator()()
	{
		int i = m_i++ & 3;
		if (!i)
			fill_buf();
		return buf[i];
	}
	template< typename T>
	A_AES std::enable_if_t<std::is_same_v<T, __m128i>, __m128i> get()
	{
		__m128i v = _mm_loadu_si128((__m128i*)m_v), k = _mm_loadu_si128((__m128i*)m_k);
		__m128i u = _mm_aesenc_si128(_mm_aesenc_si128(_mm_aesenc_si128(v, k), k), k);
		_mm_storeu_si128((__m128i*)m_v, _mm_add_epi64(v, _mm_set_epi64x(0, 1)));
		return u;
	}
protected:
	NO_INLINE A_AES void fill_buf()
	{
		_mm_storeu_si128((__m128i*)buf, get<__m128i>());
	}
private:
	uint32_t buf[4];
	uint64_t m_v[2];
	uint8_t m_k[16];
	int m_i = 0;
};

struct dist
{
	virtual ~dist() = default;
	virtual float operator()(drbg& g) = 0;
	virtual __m128 get(drbg& g) { ASSUME(0);/*?*/ }
};

struct impl
{
	std::unique_ptr<dist> d;
	impl(int, int, float, float);
	void(*add)(size_t, size_t, ptrdiff_t, void*, drbg&, dist&, float);
	float sigma;
};

}
#endif
