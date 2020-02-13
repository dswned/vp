#pragma once
#include "config.h"

#include <mutex>
#include <algorithm>
#include <memory>
#include <condition_variable>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <iterator>
#include <xxhash.h>
#include <vapoursynth.h>

#include "uf.h"

extern const VSAPI* vsapi;
extern std::vector<void(*)(VSRegisterFunction, VSPlugin*)> vregf;

struct vprop
{
	template< typename T>
	static typename std::enable_if_t<std::is_same_v<T, std::string> ||
		std::is_same_v<T, const char*>, const char*>
		get(const VSMap* map, const char* key, int index, int* error)
	{
		return vsapi->propGetData(map, key, index, error);
	}
	template< typename T>
	static typename std::enable_if_t<std::is_integral_v<T>, int64_t>
		get(const VSMap* map, const char* key, int index, int* error)
	{
		return vsapi->propGetInt(map, key, index, error);
	}
	template< typename T>
	static typename std::enable_if_t<std::is_integral_v<T>, const int64_t*>
		get(const VSMap* map, const char* key, int* error)
	{
		return vsapi->propGetIntArray(map, key, error);
	}
	template< typename T>
	static typename std::enable_if_t<std::is_floating_point_v<T>, double>
		get(const VSMap* map, const char* key, int index, int* error)
	{
		return vsapi->propGetFloat(map, key, index, error);
	}
	template< typename T>
	static typename std::enable_if_t<std::is_floating_point_v<T>, const double*>
		get(const VSMap* map, const char* key, int* error)
	{
		return vsapi->propGetFloatArray(map, key, error);
	}
	template< typename T>
	static typename std::enable_if_t<std::is_same_v<T, VSNodeRef*>, T>
		get(const VSMap* map, const char* key, int index, int* error)
	{
		return vsapi->propGetNode(map, key, index, error);
	}
	template< typename T>
	static typename std::enable_if_t<std::is_same_v<T, VSFuncRef*>, T>
		get(const VSMap* map, const char* key, int index, int* error)
	{
		return vsapi->propGetFunc(map, key, index, error);
	}
};

struct vmap
{
	vmap(const VSMap* map)
		: map(map) {}
	bool contains(const char* key)
	{
		return vsapi->propNumElements(map, key) >= 0;
	}
	int num_elements(const char* key)
	{
		return vsapi->propNumElements(map, key);
	}
	int data_size(const char* key, int index, int* error)
	{
		return vsapi->propGetDataSize(map, key, index, error);
	}
	template< typename T = int64_t>
	T get(const char* key)
	{
		int error;
		auto _ = vprop::get<T>(map, key, 0, &error);
		return T(_);
	}
	template< typename T>
	T get(const char* key, T&& def)
	{
		int error;
		auto _ = vprop::get<T>(map, key, 0, &error);
		return error ? def : T(_);
	}
	template< typename T>
	T get(const char* key, T& def)
	{
		int error;
		auto _ = vprop::get<T>(map, key, 0, &error);
		return error ? def : T(_);
	}
	template< typename T>
	int get(T& dst, const char* key)
	{
		int error;
		auto _ = vprop::get<T>(map, key, 0, &error);
		if (!error)
			dst = T(_);
		return error;
	}
	template< typename T>
	int get(std::vector<T>& dst, const char* key)
	{
		int error, n = num_elements(key);
		auto p = vprop::get<T>(map, key, &error);
		if (error)
			return error;
		if constexpr (sizeof(*p) != sizeof(T))
			dst.clear(), std::copy_n(p, n, std::back_inserter(dst));
		else
			dst.assign(p, p + n);
		return error;
	}
	int get(std::vector<std::string>& dst, const char* key)
	{
		int error = -1, k = num_elements(key);
		if (k)
			dst.clear();
		for (int i = 0; i < k; i++)
		{
			int n = data_size(key, i, &error);
			if (error)
				break;
			auto p = vprop::get<const char*>(map, key, i, &error);
			if (error)
				break;
			dst.emplace_back(std::string(p, n));
		}
		return error;
	}
private:
	const VSMap* map;
};

struct vfunc
{
	vfunc(VSFuncRef* f)
		: f(f)
		, map(vsapi->createMap())
	{
	}
	~vfunc()
	{
		vsapi->freeMap(map);
		vsapi->freeFunc(f);
	}
	vmap call()
	{
		vsapi->callFunc(f, 0, map, 0, 0);
		return vmap(map);
	}
private:
	VSFuncRef* f;
	VSMap* map;
};

struct vfplane
{
	uint8_t* p;
	int h, w, css, rowsize, pitch;
	vfplane()
		: p(0), h(0), w(0), css(0), rowsize(0), pitch(0) {}
	vfplane(int h, int w, int css)
		: p(0), h(h), w(w), css(css), rowsize(w << css), pitch(rowsize) {}
	vfplane(uint8_t* p, int h, int w, int css, int pitch)
		: p(p), h(h), w(w), css(css), rowsize(w << css), pitch(pitch) {}
	template< typename T>
	T& at(int i, int j) const { return reinterpret_cast<T*>(p + (size_t)i * pitch)[j]; }
#if 1
	template< typename T, typename F>
	void foreach(F f)
	{
		uint8_t* pdst = p;
		for (int i = 0; i < h; i++)
		{
			for (int j = 0; j < w; j++)
			{
				T& x = reinterpret_cast<T*>(pdst)[j];
				x = f(x);
			}
			pdst += pitch;
		}
	}
#else
	template< typename T, typename F>
	void foreach(F f)
	{
		for (int i = 0; i < h; i++)
		{
			for (int j = 0; j < w; j++)
				at<T>(i, j) = f(this, i, j);
		}
	}
#endif
	void copy_to(vfplane& dstp) const
	{
		if (w != dstp.w || h != dstp.h || css != dstp.css)
			return;
		if (pitch != dstp.pitch)
		{
			for (int i = 0; i < h; i++)
				std::memcpy(dstp.p + (size_t)i * dstp.pitch, p + (size_t)i * pitch, rowsize);
		}
		else
			std::memcpy(dstp.p, p, (size_t)h * pitch);
	}
	void copy_to(vfplane& dstp, bool is_even) const
	{
		if (w != dstp.w || h != dstp.h || css != dstp.css)
			return;
		for (int i = is_even; i < h; i += 2)
			std::memcpy(dstp.p + (size_t)i * dstp.pitch, p + (size_t)i * pitch, rowsize);
	}
};

struct vf
{
	std::vector<vfplane> plane;
	vf() {}
	vf(const std::vector<vfplane>& plane)
		: plane(plane) {}
	template< typename T>
	operator T();
	uint8_t* p(int p) const { return plane[p].p; }
	int pitch(int p) const { return plane[p].pitch; }
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
	const VSFrameRef* f;
	vvf()
		: f(0) {}
	vvf(const VSFrameRef* f)
		: f(f)
	{
		auto ff = vsapi->getFrameFormat(f);
		for (int p = 0; p < ff->numPlanes; p++)
		{
			plane.emplace_back(
				const_cast<uint8_t*>(vsapi->getReadPtr(f, p)),
				vsapi->getFrameHeight(f, p), vsapi->getFrameWidth(f, p),
				uf::bsr(ff->bytesPerSample), vsapi->getStride(f, p));
		}
	}
	vvf(VSFrameRef* f)
		: f(0)
	{
		auto ff = vsapi->getFrameFormat(f);
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
	void reset(const VSFrameRef* f)
	{
		for (int p = 0; p < plane.size(); p++)
			plane[p].p = (uint8_t*)vsapi->getReadPtr(f, p);
		std::swap(this->f, f);
		if (f)
			vsapi->freeFrame(f);
	}
};

template<>
inline vf::operator const VSFrameRef* ()
{
	return static_cast<vvf*>(this)->f;
}

typedef std::shared_ptr<vvf> vvf_t;

struct vvfcache
{
	const VSVideoInfo* vi;
	int nf, np, h, w, ssh, ssw, css;
	vvfcache(VSNodeRef* node)
		: node(node)
	{
		vi = vsapi->getVideoInfo(node);
		nf = vi->numFrames;
		np = vi->format->numPlanes;
		h = vi->height;
		w = vi->width;
		ssh = vi->format->subSamplingH;
		ssw = vi->format->subSamplingW;
		css = uf::bsr(vi->format->bytesPerSample);
	}
	~vvfcache()
	{
		vsapi->freeNode(node);
	}
	void request(int n, VSFrameContext* frameCtx)
	{
		vsapi->requestFrameFilter(n, node, frameCtx);
	}
	void get(vf_t& dst, int n, VSFrameContext* frameCtx)
	{
		auto f = vsapi->getFrameFilter(n, node, frameCtx);
		if (!f)
			throw "";
		dst = std::make_shared<vvf>(f);
	}
	operator VSNodeRef* () const { return node; }
private:
	VSNodeRef* node;
};

typedef std::unique_ptr<vvfcache> vvfcache_t;

namespace uf {
template< typename T, typename... U>
void free(void* p, U...)
{
	delete static_cast<T*>(p);
}
inline VSFrameRef* copy(vf_t& dst, vf_t& src, VSCore* core)
{
	auto f = vsapi->copyFrame(*src, core);
	if (!f)
		throw "";
	dst = std::make_shared<vvf>(f);
	return f;
}
inline VSFrameRef* create(vf_t& dst, int h, int w, const VSFormat* ff, const VSFrameRef* src, VSCore* core)
{
	auto f = vsapi->newVideoFrame(ff, w, h, src, core);
	if (!f)
		throw "";
	dst = std::make_shared<vvf>(f);
	return f;
}
inline uint64_t xxh_value(const char* p, const char* end)
{
	XXH64_state_t state;
	XXH64_reset(&state, 0);
	XXH64_update(&state, p, end - p);
	return XXH64_digest(&state);
}
inline uint64_t xxh_value(vfplane& srcp)
{
	XXH64_state_t state;
	XXH64_reset(&state, 0);
	if (srcp.rowsize != srcp.pitch)
	{
		for (int i = 0; i < srcp.h; i++)
			XXH64_update(&state, srcp.p + (size_t)i * srcp.pitch, srcp.rowsize);
	}
	else
		XXH64_update(&state, srcp.p, (size_t)srcp.h * srcp.pitch);
	return XXH64_digest(&state);
}
}