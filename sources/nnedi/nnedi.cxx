#include "config.h"
#if defined(BUILD_NNEDI)

#include <boost/compute.hpp>
#include <lz4.h>

#include "filter.h"

typedef void(*edi_t)(size_t, void*, const void*, const void*);
typedef std::remove_pointer_t<edi_t> edi_f;

extern "C" {
	edi_f nnedi_4x12x16_ivb;
}

namespace nnedi
{

struct impl
{
	uf::copy_pad_t copy_pad;
	uf::copy_t copy;
	std::stringstream log;
	virtual ~impl() = default;
	virtual void process(size_t, size_t, void*, ptrdiff_t, const void*, ptrdiff_t, void*, int) = 0;
	size_t tmp_size = 0;
};

struct context
{
	vf_t srcf;
	uf::aligned_unique<> tmp;
	context(const impl* d)
	{
		if (d->tmp_size)
			tmp = uf::make_aligned_unique(d->tmp_size, 64);
	}
};

struct filter
{
	static const char name[], args[];
	VSVideoInfo vi;
	VSFilterDependency deps[1];
	std::mutex mutex;
	std::vector<std::unique_ptr<context>> cc;
	filter(vmap&&);
	~filter()
	{
	}
	void proc(context*, vf*);
	std::unique_ptr<context> get();
	vvfcache_t srcc;
	std::unique_ptr<impl> nnedi;
	int field;
	int gpu;
};

#include "nnedi_q1_cl.h"

struct nnedi_q1 : public impl
{
	static constexpr size_t h_pad = 4, w_pad = 8;
	boost::compute::device device;
	boost::compute::context context;
	boost::compute::command_queue queue;
	std::vector<boost::compute::buffer> buf, filter_buf, bias_buf;
	std::vector<boost::compute::kernel> kernel;
	int read_weights(const std::string& path, uint8_t* data, size_t n);
	nnedi_q1(size_t h, size_t w, sample_type_e st, const std::string& weights_path, int gpu) try
	{
		switch (st)
		{
		case sample_type_e::BYTE:
			copy = uf::copy<float, uint8_t>;
			copy_pad = uf::copy_pad_replicate<uint8_t, float>;
			break;
		case sample_type_e::WORD:
			copy = uf::copy<float, uint16_t>;
			copy_pad = uf::copy_pad_replicate<uint16_t, float>;
			break;
		case sample_type_e::FLOAT:
			copy = uf::copy<float, float>;
			copy_pad = uf::copy_pad_replicate<float, float>;
			break;
		default:
			throw "unsupported format"s;
		}
		constexpr size_t weights_n = 575232, n = 32;
		std::vector<uint8_t> weights(weights_n * sizeof(float));
		if (read_weights(weights_path, weights.data(), weights_n))
			throw "reading weights error"s;
		size_t padded_h = h + 2, padded_w = w + 2;
		std::vector<boost::compute::device> devices = boost::compute::system::devices();
		if (devices.empty())
			throw "opencl device is not available"s;
		device = boost::compute::system::default_device();
		cl_device_id default_id = device.id();
		for (int i = 0; i < devices.size();)
		{
			boost::compute::device& d = devices[i++];
			bool enabled = i == gpu;
			if (enabled)
				device = d;
			log << (enabled ? "+" : "-") << "device[" << i << "]"
				<< (d.id() == default_id ? " (default)" : "") << ": " << d.get_info<CL_DEVICE_NAME>() << "\n";
			log << "device vendor: " << d.get_info<CL_DEVICE_VENDOR>() << "\n";
			log << "compute units: " << d.get_info<CL_DEVICE_MAX_COMPUTE_UNITS>() << "\n";
			log << "clock frequency: " << d.get_info<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << "\n";
			log << "global memory size: " << d.get_info<CL_DEVICE_GLOBAL_MEM_SIZE>() / 1048576 << " MB\n";
			log << "max memory alloc size = " << d.get_info<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() / 1048576 << " MB\n";
			log << "global memory cache size: " << d.get_info<CL_DEVICE_GLOBAL_MEM_CACHE_SIZE>() / 1024 << " KB\n";
			log << "local memory size: " << d.get_info<CL_DEVICE_LOCAL_MEM_SIZE>() / 1024 << " KB\n";
			log << "constant memory size: " << d.get_info<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>() / 1024 << " KB\n";
			log << "max work group size: " << d.get_info<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << "\n";
			log << "opencl version: " << d.get_info<CL_DEVICE_OPENCL_C_VERSION>() << "\n";
			log << "extensions: " << d.get_info<CL_DEVICE_EXTENSIONS>() << "\n\n";
		}
		context = boost::compute::context(device);
		queue = boost::compute::command_queue(context, device);
		std::vector<char> source(65536);
		if (LZ4_decompress_safe((const char*)nnedi_q1_cl, source.data(), sizeof(nnedi_q1_cl), source.size() - 1) < 0)
			throw "program decompress error"s;
		boost::compute::program program = boost::compute::program::build_with_source(source.data(), context);
		log << "program build log:\n" << program.build_log() << "\n";
		size_t padded_size = padded_h * padded_w * n * sizeof(float), dst_size = h * w * sizeof(float);
		buf.reserve(6);
		log << "create io buffers size " << padded_size * 5 + dst_size << " bytes\n";
		for (int i = 0; i < 5; i++)
			buf.emplace_back(context, padded_size, CL_MEM_READ_WRITE);
		buf.emplace_back(context, dst_size, CL_MEM_READ_WRITE);
		size_t f_size = 3 * 3, f0_size = f_size * n, f1_size = f0_size * n;
		uint8_t* ptr = weights.data();
		cl_mem_flags flags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR | CL_MEM_HOST_NO_ACCESS;
		filter_buf.reserve(68);
		log << "create kernel buffers size " << weights.size() << " bytes\n";
		for (int i = 0; i < 2; i++)
		{
			size_t size = 4 * 8 * n * sizeof(float);
			filter_buf.emplace_back(context, size, flags, ptr), ptr += size;
		}
		for (int i = 0; i < 62; i++)
		{
			size_t size = f1_size * sizeof(float);
			filter_buf.emplace_back(context, size, flags, ptr), ptr += size;
		}
		for (int i = 62; i < 66; i++)
		{
			size_t size = f0_size * sizeof(float);
			filter_buf.emplace_back(context, size, flags, ptr), ptr += size;
		}
		bias_buf.reserve(20);
		for (int i = 0; i < 20; i++)
		{
			size_t size = n * sizeof(float);
			bias_buf.emplace_back(context, size, flags, ptr), ptr += size;
		}
		kernel.reserve(9);
		kernel.emplace_back(program, "conv2d_relu");
		kernel.emplace_back(program, "conv2d_add_output_relu");
		kernel.emplace_back(program, "conv2d_add_input");
		kernel.emplace_back(program, "conv2d");
		kernel.emplace_back(program, "conv2d_add_output");
		kernel.emplace_back(program, "conv2d_add_output_reduce");
		kernel.emplace_back(program, "conv2d_add_input_reduce");
		kernel.emplace_back(program, "pad");
		kernel.emplace_back(program, "nnedi");
		kernel[8].set_arg(0, buf[1]);
		kernel[8].set_arg(1, buf[0]);
		kernel[8].set_arg(2, buf[5]);
		kernel[8].set_arg(3, filter_buf[0]);
		kernel[8].set_arg(4, filter_buf[1]);
	}
	catch (const boost::compute::program_build_failure& e)
	{
		log << e.what() << "\n" << e.build_log();
		throw log.str();
	}
	catch (const boost::compute::opencl_error& e)
	{
		log << e.what();
		throw log.str();
	}
	~nnedi_q1()
	{
#if 0
		FILE* fp = fopen("nnedi.log", "wb");
		std::string str = log.str();
		fwrite(str.data(), 1, str.size(), fp);
		fclose(fp);
#endif
	}
	void process(size_t h, size_t w,void* dst_p, ptrdiff_t dst_stride1,
		const void* src_p, ptrdiff_t src_stride1, void*, int field) try
	{
		size_t global_work_size[2], local_work_size[2], pad_global_work_size[1];
		global_work_size[0] = w;
		global_work_size[1] = h;
		local_work_size[0] = std::gcd(w, 8);
		local_work_size[1] = std::gcd(h, 8);
		pad_global_work_size[0] = std::max(h + 1, w + 1);
		const uint8_t rnds[] = {
			0, 0, 1,
			3, 0, 2, 1, 1, 2,
			3, 0, 4, 4, 1, 4, 1, 2, 4,
			3, 0, 3, 4, 1, 3, 4, 2, 3, 1, 4, 3,
			3, 0, 1, 4, 3, 1, 4, 2, 1, 1, 4, 1,
			3, 0, 2, 4, 3, 2, 4, 1, 2, 1, 4, 2,
			3, 2, 4, 4, 3, 4, 4, 1, 4, 2, 0, 4,
			0, 4, 1,
			3, 4, 2, 1, 1, 2,
			3, 4, 0, 4, 1, 0, 1, 2, 0,
			3, 4, 3, 4, 1, 3, 4, 2, 3, 1, 0, 3,
			3, 4, 1, 4, 3, 1, 4, 2, 1, 1, 0, 1,
			3, 4, 2, 4, 3, 2, 4, 1, 2, 1, 0, 2,
			3, 2, 0, 4, 3, 0, 4, 1, 0, 2, 4, 0,
			0, 0, 1,
			3, 0, 2, 1, 1, 2,
			3, 0, 4, 4, 1, 4, 1, 2, 4,
			3, 0, 3, 4, 1, 3, 4, 2, 3, 1, 4, 3,
			3, 0, 1, 4, 3, 1, 4, 2, 1, 1, 4, 1,
			3, 0, 2, 4, 3, 2, 4, 1, 2, 1, 4, 2,
			5, 2, 5, 5, 3, 5, 5, 1, 5, 6, 0, 5,
		}, * rnds_p = rnds;
		kernel[7].set_arg(1, static_cast<int>(h + 1));
		kernel[7].set_arg(2, static_cast<int>(w + 1));
		kernel[7].set_arg(3, static_cast<int>(w + 2 << 3)); // 32/float4
		void* ptr = queue.enqueue_map_buffer(buf[1], CL_MAP_WRITE, 0, buf[1].size());
		copy_pad(h + h_pad, w + w_pad, 0, ptr, h, w, src_stride1, src_p, h_pad / 2 - field, w_pad / 2 - 1);
		queue.enqueue_unmap_buffer(buf[1], ptr);
		queue.enqueue_nd_range_kernel(kernel[8], 2, 0, global_work_size, local_work_size);
		boost::compute::buffer* fbuf = filter_buf.data() + 2, * bbuf = bias_buf.data();
		for (int fi = 2, bi = 0, i = 0; i < 66; i++)
		{
			size_t j = *rnds_p++, a = *rnds_p++, b = *rnds_p++;
			if (j < 1)
			{
				kernel[7].set_arg(0, buf[a]);
				queue.enqueue_nd_range_kernel(kernel[7], 1, 0, pad_global_work_size, 0);
			}
			boost::compute::kernel& k = kernel[j];
			k.set_arg(0, buf[a]);
			k.set_arg(1, buf[b]);
			k.set_arg(2, *fbuf++);
			if (j < 3)
				k.set_arg(3, *bbuf++);
			queue.enqueue_nd_range_kernel(k, 2, 0, global_work_size, 0);
			if (j < 2)
			{
				kernel[7].set_arg(0, buf[b]);
				queue.enqueue_nd_range_kernel(kernel[7], 1, 0, pad_global_work_size, 0);
			}
			if (j < 3)
				queue.finish();
		}
		ptr = queue.enqueue_map_buffer(buf[5], CL_MAP_READ, 0, buf[5].size());
		copy(h, w, dst_stride1, dst_p, h, w, 0, ptr);
		queue.enqueue_unmap_buffer(buf[5], ptr);
	}
	catch (const boost::compute::opencl_error& e)
	{
		throw e.error_string();
	}
};

int nnedi_q1::read_weights(const std::string& path, uint8_t* data, size_t n)
{
	uf::file_h f;
	if (!f.open(path.c_str()))
		throw "can't open weights file"s;
	size_t size = f.size;
	if (size == n * sizeof(float))
	{
		if (f.read(data, size) != size)
			return -1;
	}
	else if (size == n * sizeof(uf::sf16))
	{
		std::vector<uf::sf16> tmp(n);
		if (f.read(tmp.data(), size) != size)
			return -1;
		uf::conv<uf::sf16, float>::transform_n(tmp.data(), tmp.size(), reinterpret_cast<float*>(data));
	}
	else
		throw "invalid weights file"s;
	return 0;
}

struct nnedi_q0 : public impl
{
	static constexpr size_t h_pad = 4, w_pad = 12;
	alignas(64) float nnedi_data[1616];
	nnedi_q0(size_t h, size_t w, sample_type_e st, const std::string& weights_path, float pscr)
	{
		switch (st)
		{
		case sample_type_e::BYTE:
			copy = uf::copy<float, uint8_t>;
			copy_pad = uf::copy_pad_replicate<uint8_t, float>;
			break;
		case sample_type_e::WORD:
			copy = uf::copy<float, uint16_t>;
			copy_pad = uf::copy_pad_replicate<uint16_t, float>;
			break;
		case sample_type_e::FLOAT:
			copy = uf::copy<float, float>;
			copy_pad = uf::copy_pad_replicate<float, float>;
			break;
		default:
			throw "unsupported format"s;
		}
		tmp_size = (h + h_pad) * (w + w_pad) * sizeof(float);
		if (read_weights(weights_path))
			throw "reading weights error"s;
		*nnedi_data *= pscr;
	}
	int read_weights(const std::string& path)
	{
		uf::file_h f;
		if (!f.open(path.c_str()))
			throw "can't open weights file"s;
		size_t size = f.size;
		if (size == sizeof(nnedi_data))
		{
			if (f.read(nnedi_data, size) != size)
				return -1;
		}
		else
			throw "invalid weights file"s;
		return 0;
	}
	void process(size_t h, size_t w, void* dst_p, ptrdiff_t dst_stride1,
		const void* src_p, ptrdiff_t src_stride1, void* tmp_p, int field)
	{
		float* tmp = static_cast<float*>(tmp_p);
		float* dst = static_cast<float*>(dst_p);
		ptrdiff_t dst_stride = dst_stride1 / sizeof(float);
		ptrdiff_t padded_stride = w + w_pad;
		copy_pad(h + h_pad, w + w_pad, 0, tmp_p, h, w, src_stride1, src_p, h_pad / 2 - field, w_pad / 2 - 1);
		for (unsigned i = 0; i < h; i++)
			nnedi_4x12x16_ivb(w * sizeof(float), tmp + i * w, tmp + i * padded_stride, nnedi_data);
		copy(h, w, dst_stride1, dst_p, h, w, 0, tmp_p);
	}
};

filter::filter(vmap&& map)
{
	VSNode* node = map.get<VSNode*>("clip");
	vi = *vsapi->getVideoInfo(node);
	srcc = std::make_unique<vvfcache>(node);
	field = map.get("field");
	int q = map.get("q", 0);
	if (q < 0 || q > 1)
		throw "q: 0, 1"s;
	std::string weights;
	if (!map.get(weights, "weights") || weights.compare(0, 2, "./"))
		weights = get_plugins_path() + (weights.empty() ? (q == 0 ? "nnedi_q0.bin" : "nnedi_q1.bin") : weights);
	float pscr = map.get("pscr", 1.f);
	gpu = map.get("gpu", -1); // [-1 auto 1.. mask 0 disabled]
	int h = srcc->h(), w = srcc->w();
	if (field < 0 || field > 2)
		throw "field: 0 = pred even lines, 1 = odd, 2 = both. [3 = auto, ? = upscale]"s;
	else
	{
		if (q == 0)
			nnedi = std::make_unique<nnedi_q0>(h / 2, w, srcc->st, weights, pscr);
		if (q == 1)
			nnedi = std::make_unique<nnedi_q1>(h / 2, w, srcc->st, weights, gpu);
	}
	deps[0] = { node, rpStrictSpatial };
}

std::unique_ptr<context> filter::get()
{
	return std::make_unique<context>(nnedi.get());
}

void filter::proc(context* ctx, vf* dstf)
{
	for (int p = 0; p < srcc->np; p++)
	{
		vfplane& dstp = dstf->plane[p], & srcp = ctx->srcf->plane[p];
		size_t src_stride = srcp.stride, dst_stride = dstp.stride;
		uint8_t* src_p = srcp.p, * dst_p = dstp.p;
		int h = srcp.h, w = srcp.w;
		if (field + 1 & 1)
			nnedi->process(h / 2, w, dst_p, dst_stride << 1,
				src_p + src_stride, src_stride << 1, ctx->tmp.get(), 0);
		else
			srcp.copy_to(dstp, 0);
		if (field + 1 & 2)
			nnedi->process(h / 2, w, dst_p + dst_stride, dst_stride << 1,
				src_p, src_stride << 1, ctx->tmp.get(), 1);
		else
			srcp.copy_to(dstp, 1);
	}
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
	vvf dstf;
	try
	{
		if (!ctx)
			ctx = d->get();
		d->srcc->get(ctx->srcf, n, frameCtx);
		dstf.create(*ctx->srcf, core);
		d->proc(ctx.get(), &dstf);
	}
	catch (const std::string& e)
	{
		vsapi->setFilterError(e.c_str(), frameCtx);
		dstf.reset();
	}
	if (ctx)
	{
		lock.lock();
		d->cc.emplace_back(std::move(ctx));
		lock.unlock();
	}
	return dstf.release();
}

void create(const VSMap* in, VSMap* out, void*, VSCore* core, const VSAPI* vsapi)
{
	if (!::vsapi)
		::vsapi = vsapi;
	try
	{
		filter* d = new filter(in);
		vsapi->createVideoFilter(out, filter::name, &d->vi, get, uf::free<filter>, d->gpu ? fmParallelRequests : fmParallel,
			d->deps, std::size(d->deps), d, core);
	}
	catch (const std::string& e)
	{
		vsapi->mapSetError(out, e.c_str());
	}
}

const char filter::name[] = "nnedi";
const char filter::args[] =
"clip:vnode;"
"field:int;"
"q:int:opt;"
"weights:data:opt;"
"pscr:float:opt;"
"gpu:int:opt;";

void reg_f(VSPlugin* plugin, const VSPLUGINAPI* vspapi)
{
	vspapi->registerFunction(filter::name, filter::args, "clip:vnode;", create, 0, plugin);
}

PUSH_REG_F(reg_f);
}
#endif