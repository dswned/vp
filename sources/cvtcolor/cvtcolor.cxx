#include "config.h"
#if defined(BUILD_CVTCOLOR)

#include "filter.h"
#include "cpu.h"
#include "cvtcolor.h"

namespace cvtcolor {

enum { CS_NONE, CS_YUV, CS_LAB, CS_COMP };

cvt_t get_cvt_f(int cs_in, int cs_out, sample_type_e st, uf::cpu cpu)
{
	typedef uf::vfk_s<cvt_t, uf::cpuflags,/*+initfunc,etc*/ std::tuple<int, int, sample_type_e>> vfk;
	constexpr vfk map[] =
	{
		vfk(lab2yuv_f32_avx512, uf::cpuflags::avx512f, CS_LAB, CS_YUV, sample_type_e::FLOAT),
		vfk(yuv2lab_f32_ivb, uf::cpuflags::avx, CS_YUV, CS_LAB, sample_type_e::FLOAT),
		vfk(lab2yuv_f32_ivb, uf::cpuflags::avx, CS_LAB, CS_YUV, sample_type_e::FLOAT),
		vfk(yuv2lab_i16_avx, uf::cpuflags::avx, CS_YUV, CS_LAB, sample_type_e::WORD),
		vfk(lab2yuv_i16_avx, uf::cpuflags::avx, CS_LAB, CS_YUV, sample_type_e::WORD),
	};
	return vfk::select(std::begin(map), std::end(map), cpu.flags, std::make_tuple(cs_in, cs_out, st));
}

struct filter
{
	static const char name[], args[];
	VSVideoInfo vi;
	VSFilterDependency deps[1];
	std::unique_ptr<vvfcache> srcc;
	cvt_t cvt;
	filter(vmap&& map)
	{
		auto translate_cs = [](const char* p)
		{
			if (!strcmp("yuv", p))
				return CS_YUV;
			if (!strcmp("lab", p))
				return CS_LAB;
			if (!strcmp("comp", p))
				return CS_COMP;
			return CS_NONE;
		};
		int opt = map.get("opt", -1);
		std::string params = std::string(map.get<std::string>("params"));
		char* y = 0, * x = strtok_r(params.data(), ">", &y);
		if (!x || y && !*y)
			throw "params parse failure { format '%s>%s' }"s;
		VSNode* node = map.get<VSNode*>("clip");
		vi = *vsapi->getVideoInfo(node);
		srcc = std::make_unique<vvfcache>(node);
		int cs_in = translate_cs(x), cs_out = translate_cs(y);
		switch (cs_in)
		{
		case CS_YUV:
		case CS_LAB:
			if (srcc->np != 3 || srcc->ssh || srcc->ssw)
				throw "input must be 3-plane w/o chroma subsampling"s;
			break;
		case CS_COMP:
			break;
		}
		uf::cpu cpu = opt < 0 ? uf::query_x86_capabilities() : uf::cpu(opt);
		cvt = get_cvt_f(cs_in, cs_out, srcc->st, cpu);
		if (!cvt)
			throw "cvt function not found"s;
		if (cvt == yuv2lab_i16_avx)
			init_lab_cbrt_tab();
		deps[0] = { node, rpStrictSpatial };
	}
	void proc(vf* dstf, vf* srcf)
	{
		size_t n = srcf->h() * srcf->stride(0);
		uint8_t* sy = srcf->ptr(0), * su = srcf->ptr(1), * sv = srcf->ptr(2);
		uint8_t* dy = dstf->ptr(0), * du = dstf->ptr(1), * dv = dstf->ptr(2);
		cvt(n, dy, du, dv, sy, su, sv);
	}
};

const VSFrame* get(int n, int activationReason,
	void* instanceData, void**, VSFrameContext* frameCtx, VSCore* core, const VSAPI*)
{
	filter* d = static_cast<filter*>(instanceData);
	if (activationReason == arInitial)
		d->srcc->request(n, frameCtx);
	if (activationReason != arAllFramesReady)
		return 0;
	vvf dstf, srcf;
	d->srcc->get(srcf, n, frameCtx);
	dstf.create(srcf, core);
	d->proc(&dstf, &srcf);
	return dstf.release();
}

void create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi)
{
	if (!::vsapi)
		::vsapi = vsapi;
	try
	{
		filter* d = new filter(in);
		vsapi->createVideoFilter(out, filter::name, &d->vi, get, uf::free<filter>, fmParallel, d->deps, std::size(d->deps), d, core);
	}
	catch (const std::string& e)
	{
		vsapi->mapSetError(out, e.c_str());
	}
}

const char filter::name[] = "cvtcolor";
const char filter::args[] = "clip:vnode;params:data;opt:int:opt;";

void reg_f(VSPlugin* plugin, const VSPLUGINAPI* vspapi)
{
	vspapi->registerFunction(filter::name, filter::args, "clip:vnode;", create, 0, plugin);
}

PUSH_REG_F(reg_f);
}
#endif