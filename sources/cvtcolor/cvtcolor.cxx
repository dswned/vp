#include "config.h"
#if defined(BUILD_CVTCOLOR)

#include "cpu.h"
#include "filter.h"

typedef void(*cvt_f)(size_t, void*, void*, void*, void*, void*, void*);
extern "C" std::remove_pointer_t<cvt_f> yuv2lab_709_32f_ivb, lab2yuv_709_32f_ivb;

namespace cvtcolor {

enum { NONE, YUV, RGB, LAB };

cvt_f get_cvt_f(int x, int y, uf::cpu& cpu)
{
	const std::map<std::tuple<int, int>, cvt_f> map =
	{
		{ std::make_tuple(YUV, LAB), yuv2lab_709_32f_ivb },
		{ std::make_tuple(LAB, YUV), lab2yuv_709_32f_ivb }
	};
	auto k = std::make_tuple(x, y);
	if (map.count(k))
		return map.at(k);
	return nullptr;
}

struct filter
{
	VSVideoInfo vi;
	static const char name[], args[];
	filter(vmap&& map)
	{
		auto cs = [](const char* p)
		{
			if (!strcmp("yuv", p))
				return YUV;
			if (!strcmp("lab", p))
				return LAB;
			return NONE;
		};
		auto params = map.get<std::string>("params");
		char* y, * x = strtok_r(params.data(), ">", &y);
		if (!x)
			throw "";
		auto cpu = uf::query_x86_capabilities();
		cvt = get_cvt_f(cs(x), cs(y), cpu);
		if (!cvt)
			throw "";
		VSNodeRef* node = map.get<VSNodeRef*>("clip");
		vi = *vsapi->getVideoInfo(node);
		if (vi.format->numPlanes != 3 || vi.format->sampleType != stFloat || vi.format->bitsPerSample != 32)
		{
			vsapi->freeNode(node);
			throw "";
		}
		srcc = std::make_unique<vvfcache>(node);
	}
	void proc(vf* dstf, vf* srcf)
	{
		size_t h = srcf->h(), pitch = srcf->pitch(0);
		uint8_t* sy = srcf->p(0), * su = srcf->p(1), * sv = srcf->p(2);
		uint8_t* dy = dstf->p(0), * du = dstf->p(1), * dv = dstf->p(2);
		cvt(h * pitch, dy, du, dv, sy, su, sv);
	}
	std::unique_ptr<vvfcache> srcc;
	cvt_f cvt;
};

const VSFrameRef* get(int n, int activationReason, void** instanceData, void** frameData, VSFrameContext* frameCtx, VSCore* core, const VSAPI* vsapi)
{
	filter* d = static_cast<filter*>(*instanceData);
	if (activationReason == arInitial)
		d->srcc->request(n, frameCtx);
	if (activationReason != arAllFramesReady)
		return 0;
	vf_t dstf, srcf;
	d->srcc->get(srcf, n, frameCtx);
	VSFrameRef* dst = uf::create(dstf, d->vi.height, d->vi.width, d->vi.format, *srcf, core);
	d->proc(dstf.get(), srcf.get());
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

const char filter::name[] = "cvtcolor";
const char filter::args[] = "clip:clip;params:data;";

struct reg
{
	reg() {
		vregf.emplace_back([](VSRegisterFunction registerFunc, VSPlugin* plugin) {
			registerFunc(filter::name, filter::args, create, plugin, plugin); });
	}
} _;
}
#endif