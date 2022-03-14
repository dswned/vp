#include "config.h"
#if defined(BUILD_RAWS)
#include "filter.h"

#include <string.h>

namespace raws {

struct filter
{
	static const char name[], args[];
	VSVideoInfo vi;
	filter(vmap&&, VSCore*);
	void proc(int n, vf* dstf);
	uf::file_h fh;
	size_t frame_size;
	int order[4];
	uint8_t** pptr, * ptr;
	std::vector<uint64_t> index;
	std::unique_ptr<uint8_t[]> buf;
};

void write_planar_frame(uint8_t* psrc, int* order, vf* dstf)
{
	size_t np = dstf->plane.size();
	for (int rp = 0; rp < np; rp++)
	{
		int p = order[rp];
		vfplane& dstp = dstf->plane[p];
		for (size_t i = 0; i < dstp.h; i++)
		{
			std::copy_n(psrc, dstp.rowsize, dstp.p + i * dstp.stride);
			psrc += dstp.rowsize;
		}
	}
}

filter::filter(vmap&& map, VSCore* core)
{
#define PACK(a,b,c,d) (a|b<<2|c<<4|d<<6)
	const struct
	{
		const char* format_name;
		uint8_t num_planes;
		uint8_t order;
		VSColorFamily cf;
		VSSampleType st;
		uint8_t bps;
		uint8_t ssw;
		uint8_t ssh;
	} table[] = {
		{ "gray", 1, 0, cfGray, stInteger, 8, 0, 0 },
		{ "gray10", 1, 0, cfGray, stInteger, 10, 0, 0 },
		{ "gray16", 1, 0, cfGray, stInteger, 16, 0, 0 },
		{ "grayh", 1, 0, cfGray, stFloat, 16, 0, 0 },
		{ "grays", 1, 0, cfGray, stFloat, 32, 0, 0 },
		{ "yuv420p8", 3, PACK(0, 1, 2, 0), cfYUV, stInteger, 8, 1, 1 },
		{ "yuv420p10", 3, PACK(0, 1, 2, 0), cfYUV, stInteger, 10, 1, 1 },
		{ "yuv420p16", 3, PACK(0, 1, 2, 0), cfYUV, stInteger, 16, 1, 1 },
		{ "yuv444p8", 3, PACK(0, 1, 2, 0), cfYUV, stInteger, 8, 0, 0 },
		{ "yuv444p10", 3, PACK(0, 1, 2, 0), cfYUV, stInteger, 10, 0, 0 },
		{ "yuv444p16", 3, PACK(0, 1, 2, 0), cfYUV, stInteger, 16, 0, 0 },
		{ "yuv444ps", 3, PACK(0, 1, 2, 0), cfYUV, stFloat, 32, 0, 0 },
		{ "yv12", 3, PACK(0, 2, 1, 0), cfYUV, stInteger, 8, 1, 1 },
		{ "yv24", 3, PACK(0, 2, 1, 0), cfYUV, stInteger, 8, 0, 0 },
	};
	int h = map.get<int>("height");
	int w = map.get<int>("width");
	int np, ssh, ssw, bps;
	int nf;
	if (std::string format_name; !map.get(format_name, "format"))
	{
		int id = map.get<int>("format_id");
		if (!vsapi->getVideoFormatByID(&vi.format, id, core))
			throw "invalid format_id"s;
		bps = vi.format.bitsPerSample;
		np = vi.format.numPlanes;
		ssh = vi.format.subSamplingH;
		ssw = vi.format.subSamplingW;
		for (int i = 0; i < 4; i++)
			order[i] = i;
	}
	else
	{
		auto end = table + std::size(table), it = std::find_if(table, end,
			[&](auto& x) -> bool { return !strcasecmp(x.format_name, format_name.c_str()); });
		if (it == end)
			throw "unsupported format"s;
		bps = it->bps, np = it->num_planes, ssh = it->ssh, ssw = it->ssw;
		vsapi->queryVideoFormat(&vi.format, it->cf, it->st, bps, ssw, ssh, core);
		for (int i = 0, p = it->order; i < 4; i++, p >>= 2)
			order[i] = p & 0b11;
	}
	if (w >> ssw << ssw != w)
		throw "invalid width"s;
	if (h >> ssh << ssh != h)
		throw "invalid height"s;
	frame_size = 0;
	for (int p = 0; p < np; p++)
		frame_size += h * w * uf::ceil2_n(bps, 8) >> 3 >> (p ? ssh + ssw : 0);
	pptr = (uint8_t**)map.get<intptr_t>("pptr", 0);
	if (!pptr)
	{
		ptr = (uint8_t*)map.get<intptr_t>("ptr", 0);
		if (ptr)
			pptr = &ptr;
	}
	if (!pptr)
	{
		std::string fn;
		if (!map.get(fn, "file"))
			throw ""s;
		fh.open(fn.c_str());
		nf = (int)(fh.size / frame_size);
		buf.reset(new uint8_t[frame_size]);
	}
	else if (!map.get(nf, "num_frames"))
		nf = 1;
	vi.height = h;
	vi.width = w;
	vi.numFrames = nf;
	vi.fpsNum = 25;
	vi.fpsDen = 1;
	index.resize(nf);
	size_t off = 0;
	for (int i = 0; i < nf; i++)
	{
		index[i] = off;
		off += frame_size;
	}
}

void filter::proc(int n, vf* dstf)
{
	uint8_t* ptr;
	if (pptr)
		ptr = *pptr + index[n];
	else
	{
		ptr = buf.get();
		fh.read(ptr, index[n], frame_size);
	}
	write_planar_frame(ptr, order, dstf);
}

const VSFrame* get(int n, int activationReason,
	void* instanceData, void**, VSFrameContext*, VSCore* core, const VSAPI*)
{
	filter* d = static_cast<filter*>(instanceData);
	if (activationReason != arInitial)
		return 0;
	vvf dstf;
	dstf.create(&d->vi, 0, core);
	d->proc(n, &dstf);
	return dstf.release();
}

void create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi)
{
	if (!::vsapi)
		::vsapi = vsapi;
	try
	{
		filter* d = new filter(in, core);
		vsapi->createVideoFilter(out, filter::name, &d->vi, get, uf::free<filter>, fmUnordered, 0, 0, d, core);
	}
	catch (const std::string& e)
	{
		vsapi->mapSetError(out, e.c_str());
	}
}

const char filter::name[] = "raws";
const char filter::args[] =
"width:int;" "height:int;" "format:data:opt;" "file:data:opt;" "format_id:int:opt;" "ptr:int:opt;" "pptr:int:opt;" "num_frames:int:opt;";

void reg_f(VSPlugin* plugin, const VSPLUGINAPI* vspapi)
{
	vspapi->registerFunction(filter::name, filter::args, "clip:vnode;", create, 0, plugin);
}

PUSH_REG_F(reg_f);
}
#endif
