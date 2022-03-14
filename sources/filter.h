#pragma once
#include "config.h"
#include "uf.h"

#include <vapoursynth4.h>

#if defined _WIN32
#define strtok_r strtok_s
#define strcasecmp stricmp
#endif

using namespace std::literals::string_literals;

extern const VSPlugin* plugin;
extern const VSAPI* vsapi;
extern std::vector<void(*)(VSPlugin*, const VSPLUGINAPI*)> v_reg_f;

inline std::string get_plugins_path()
{
	std::string plugin_path = vsapi->getPluginPath(plugin);
	return plugin_path.substr(0, plugin_path.rfind('/') + 1);
}

template< typename, typename = void>
struct _vmap;

template<>
struct _vmap<VSNode*>
{
	using element_type = VSNode*;
	static element_type get(const VSMap* map, const char* key, int index, int* error) noexcept
	{
		return vsapi->mapGetNode(map, key, index, error);
	}
};

template<>
struct _vmap<VSFunction*>
{
	using element_type = VSFunction*;
	static element_type get(const VSMap* map, const char* key, int index, int* error) noexcept
	{
		return vsapi->mapGetFunction(map, key, index, error);
	}
};

template< typename T>
struct _vmap<T, std::enable_if_t<std::is_same_v<T, const char*>>>
{
	using element_type = const char*;
	static element_type get(const VSMap* map, const char* key, int index, int* error) noexcept
	{
		return vsapi->mapGetData(map, key, index, error);
	}
};

template< typename T>
struct _vmap<T, std::enable_if_t<std::is_floating_point_v<T>>>
{
	using element_type = double;
	static element_type get(const VSMap* map, const char* key, int index, int* error) noexcept
	{
		return vsapi->mapGetFloat(map, key, index, error);
	}
	static const element_type* get(const VSMap* map, const char* key, int* error) noexcept
	{
		return vsapi->mapGetFloatArray(map, key, error);
	}
};

template< typename T>
struct _vmap<T, std::enable_if_t<std::is_integral_v<T>>>
{
	using element_type = int64_t;
	static element_type get(const VSMap* map, const char* key, int index, int* error) noexcept
	{
		return vsapi->mapGetInt(map, key, index, error);
	}
	static const element_type* get(const VSMap* map, const char* key, int* error) noexcept
	{
		return vsapi->mapGetIntArray(map, key, error);
	}
	static_assert(sizeof(element_type) >= sizeof(T));
};

struct vmap
{
	vmap(const VSMap* map)
		: m_map(map)
	{
	}
	bool contains(const char* key) const noexcept
	{
		return vsapi->mapNumElements(m_map, key) >= 0;
	}
	int num_elements(const char* key) const noexcept
	{
		return vsapi->mapNumElements(m_map, key);
	}
	/* nothrow? */
	template< typename T>
	std::enable_if_t<std::is_same_v<T, VSNode*>, VSNode*> get(const char* key, int index) const
	{
		int error;
		auto x = _vmap<VSNode*>::get(m_map, key, index, &error);
		return x;
	}
	template< typename T>
	std::enable_if_t<std::is_same_v<T, std::string>, std::string_view> get(const char* key, int index = 0) const
	{
		int error;
		size_t n = get_data_size(key, index, &error);
		if (error)
			return std::string_view();
		auto p = _vmap<const char*>::get(m_map, key, index, 0);
		return std::string_view(p, n);
	}
	template< typename T = int64_t>
	std::enable_if_t<!std::is_convertible_v<T, std::string_view>, T> get(const char* key) const
	{
		int error;
		auto x = _vmap<T>::get(m_map, key, 0, &error);
		return static_cast<T>(x);
	}
	/**/
	template< typename T>
	[[nodiscard]] std::enable_if_t<std::is_convertible_v<T, std::string_view>, std::string_view> get(const char* key, T&& def) const
	{
		int error;
		size_t n = get_data_size(key, 0, &error);
		if (error)
			return def;
		auto p = _vmap<const char*>::get(m_map, key, 0, 0);
		return std::string_view(p, n);
	}
	template< typename T>
	[[nodiscard]] std::enable_if_t<!std::is_convertible_v<T, std::string_view>, T> get(const char* key, T&& def) const
	{
		int error;
		auto x = _vmap<T>::get(m_map, key, 0, &error);
		return error ? def : static_cast<T>(x);
	}
	/**/
	template< typename T>
	std::enable_if_t<std::is_same_v<T, std::string>, bool> get(T& dst, const char* key) const
	{
		int error;
		size_t n = get_data_size(key, 0, &error);
		if (error)
			return false;
		auto ptr = _vmap<const char*>::get(m_map, key, 0, 0);
		dst.assign(ptr, n);
		return true;
	}
	template< typename T>
	std::enable_if_t<!std::is_convertible_v<T, std::string_view>, bool> get(T& dst, const char* key) const
	{
		int error;
		auto x = _vmap<T>::get(m_map, key, 0, &error);
		if (error)
			return false;
		dst = static_cast<T>(x);
		return true;
	}
	/**/
	template< typename T>
	bool get(std::vector<T>& dst, const char* key) const
	{
		int error, n = num_elements(key);
		auto p = _vmap<T>::get(m_map, key, &error);
		if (error)
			return false;
		if constexpr (sizeof(typename _vmap<T>::element_type) != sizeof(T))
		{
			if (!dst.empty())
				dst.clear();
			dst.reserve(n);
			std::transform(p, p + n, std::back_inserter(dst), [](auto x) { return static_cast<T>(x); });
		}
		else
			dst.assign(p, p + n);
		return true;
	}
protected:
	int get_data_size(const char* key, int index, int* error) const noexcept
	{
		return vsapi->mapGetDataSize(m_map, key, index, error);
	}
private:
	const VSMap* m_map;
};

struct vfunc
{
	vfunc(VSFunction* f)
		: m_f(f)
		, m_map(vsapi->createMap())
	{
	}
	~vfunc()
	{
		vsapi->freeMap(m_map);
		vsapi->freeFunction(m_f);
	}
	vmap call()
	{
		vsapi->callFunction(m_f, m_map, m_map);
		return m_map;
	}
private:
	VSFunction* m_f;
	VSMap* m_map;
};

struct vfplane
{
	uint8_t* p;
	int h, w, css, rowsize;
	ptrdiff_t stride;
	vfplane()
		: p(0), h(0), w(0), css(0), rowsize(0), stride(0)
	{
	}
	vfplane(int h, int w, int css)
		: p(0), h(h), w(w), css(css), rowsize(w << css), stride(rowsize)
	{
	}
	vfplane(int h, int w, int css, int align)
		: p(0), h(h), w(w), css(css), rowsize(w << css), stride(uf::ceil_n(rowsize, align))
	{
	}
	vfplane(uint8_t* p, int h, int w, int css, ptrdiff_t stride)
		: p(p), h(h), w(w), css(css), rowsize(w << css), stride(stride)
	{
	}
	template< typename T>
	T& at(int i, int j) const { return reinterpret_cast<T*>(p + i * stride)[j]; }
	NO_INLINE void copy_to(vfplane& dstp) const
	{
		if (w != dstp.w || h != dstp.h || css != dstp.css)
			return;
		if (stride != dstp.stride)
		{
			for (int i = 0; i < h; i++)
				std::copy_n(p + i * stride, rowsize, dstp.p + i * dstp.stride);
		}
		else
			std::copy_n(p, h * stride, dstp.p);
	}
	NO_INLINE void copy_to(vfplane& dstp, bool is_odd) const
	{
		if (w != dstp.w || h != dstp.h || css != dstp.css)
			return;
		for (int i = is_odd; i < h; i += 2)
			std::copy_n(p + i * stride, rowsize, dstp.p + i * dstp.stride);
	}
};

typedef vfplane& vfplane_ref;

struct vf
{
	std::vector<vfplane> plane;
	vf() {}
	vf(const std::vector<vfplane>& plane)
		: plane(plane) {}
	operator const VSFrame* () const;
	operator const VSVideoFormat* () const;
	int np() const { return static_cast<int>(plane.size()); }
	uint8_t* ptr(int p) const { return plane[p].p; }
	template< typename T = uint8_t>
	ptrdiff_t stride(int p) const { return plane[p].stride / sizeof(T); }
	int h(int p = 0) const { return plane[p].h; }
	int w(int p = 0) const { return plane[p].w; }
	void copy_to(vf* dstf) const
	{
		for (int p = 0; p < plane.size(); p++)
			plane[p].copy_to(dstf->plane[p]);
	}
};

typedef std::shared_ptr<vf> vf_t;

struct vvf
	: vf
{
	const VSFrame* f;
	const VSVideoFormat* ff;
	vvf()
		: f(0)
		, ff(0)
	{
	}
	vvf(const VSFrame* f)
		: f(f)
	{
		ff = vsapi->getVideoFrameFormat(f);
		for (int p = 0; p < ff->numPlanes; p++)
		{
			plane.emplace_back(
				const_cast<uint8_t*>(vsapi->getReadPtr(f, p)),
				vsapi->getFrameHeight(f, p), vsapi->getFrameWidth(f, p),
				uf::bsr(ff->bytesPerSample), vsapi->getStride(f, p));
		}
	}
	vvf(VSFrame* f)
		: f(f)
	{
		ff = vsapi->getVideoFrameFormat(f);
		for (int p = 0; p < ff->numPlanes; p++)
		{
			plane.emplace_back(
				vsapi->getWritePtr(f, p),
				vsapi->getFrameHeight(f, p), vsapi->getFrameWidth(f, p),
				uf::bsr(ff->bytesPerSample), vsapi->getStride(f, p));
		}
	}
	~vvf()
	{
		if (f)
			vsapi->freeFrame(f);
	}
	const VSFrame* release()
	{
		return std::exchange(f, nullptr);
	}
	void reset()
	{
		if (f)
			vsapi->freeFrame(f);
		f = 0;
	}
	void reset(const VSFrame* _f)
	{
		if (f)
			vsapi->freeFrame(f);
		f = _f;
		if (plane.empty())
		{
			ff = vsapi->getVideoFrameFormat(_f);
			for (int p = 0; p < ff->numPlanes; p++)
			{
				plane.emplace_back(
					const_cast<uint8_t*>(vsapi->getReadPtr(_f, p)),
					vsapi->getFrameHeight(_f, p), vsapi->getFrameWidth(_f, p),
					uf::bsr(ff->bytesPerSample), vsapi->getStride(_f, p));
			}
		}
		else
		{
			for (int p = 0; p < plane.size(); p++)
				plane[p].p = const_cast<uint8_t*>(vsapi->getReadPtr(_f, p));
		}
	}
	void reset(VSFrame* _f)
	{
		f = _f;
		if (!plane.empty())
		{
			for (int p = 0; p < plane.size(); p++)
				plane[p].p = vsapi->getWritePtr(_f, p);
		}
		else
		{
			ff = vsapi->getVideoFrameFormat(_f);
			for (int p = 0; p < ff->numPlanes; p++)
			{
				plane.emplace_back(
					vsapi->getWritePtr(_f, p),
					vsapi->getFrameHeight(_f, p), vsapi->getFrameWidth(_f, p),
					uf::bsr(ff->bytesPerSample), vsapi->getStride(_f, p));
			}
		}
	}
	void create(const vf& srcf, VSCore* core)
	{
		VSFrame* f = vsapi->newVideoFrame(srcf, srcf.w(), srcf.h(), srcf, core);
		if (!f)
			throw "vvf create failed"s;
		reset(f);
	}
	void create(const VSVideoInfo* vi, const VSFrame* srcf, VSCore* core)
	{
		VSFrame* f = vsapi->newVideoFrame(&vi->format, vi->width, vi->height, srcf, core);
		if (!f)
			throw "vvf create failed"s;
		reset(f);
	}
	void create_copy(const vf& srcf, VSCore* core)
	{
		VSFrame* f = vsapi->copyFrame(srcf, core);
		if (!f)
			throw "vvf create failed"s;
		reset(f);
	}
};

inline vf::operator const VSFrame* () const
{
	return static_cast<const vvf*>(this)->f;
}

inline vf::operator const VSVideoFormat* () const
{
	return static_cast<const vvf*>(this)->ff;
}

typedef std::shared_ptr<vvf> vvf_t;

enum class sample_type_e : uint8_t { BYTE, WORD, HALF, FLOAT, UNKNOWN };

struct vvfcache
{
	const VSVideoInfo* vi;
	int nf, np, ssh, ssw, bps, css;
	sample_type_e st;
	vvfcache(VSNode* node)
		: m_node(node)
	{
		vi = vsapi->getVideoInfo(node);
		nf = vi->numFrames;
		np = vi->format.numPlanes;
		m_h = vi->height;
		m_w = vi->width;
		ssh = vi->format.subSamplingH;
		ssw = vi->format.subSamplingW;
		bps = vi->format.bitsPerSample;
		css = uf::bsr(vi->format.bytesPerSample);
		if (int t = vi->format.sampleType; t == stInteger && css == 0)
			st = sample_type_e::BYTE;
		else if (t == stInteger && css == 1)
			st = sample_type_e::WORD;
		else if (t == stFloat && css == 1)
			st = sample_type_e::HALF;
		else if (t == stFloat && css == 2)
			st = sample_type_e::FLOAT;
		else
			st = sample_type_e::UNKNOWN;
	}
	~vvfcache()
	{
		vsapi->freeNode(m_node);
	}
	void request(int n, VSFrameContext* frameCtx)
	{
		vsapi->requestFrameFilter(n, m_node, frameCtx);
	}
	void get(vf_t& dst, int n, VSFrameContext* frameCtx)
	{
		auto f = vsapi->getFrameFilter(n, m_node, frameCtx);
		if (!f)
			throw ""s;
		dst = std::make_shared<vvf>(f);
	}
	void get(vvf& dst, int n, VSFrameContext* frameCtx)
	{
		auto f = vsapi->getFrameFilter(n, m_node, frameCtx);
		if (!f)
			throw ""s;
		dst.reset(f);
	}
	operator VSNode* () const { return m_node; }
	int h(int p = 0) const { return m_h >> (p ? ssh : 0); }
	int w(int p = 0) const { return m_w >> (p ? ssw : 0); }
private:
	VSNode* m_node;
	int m_h, m_w;
};

typedef std::unique_ptr<vvfcache> vvfcache_t;

#define PUSH_REG_F(f) struct reg_s { reg_s() { v_reg_f.push_back(f); } } _;
