#pragma once

namespace uf {

enum cpuflags
{
	NONE,
	rdrand = 1 << 0,
	rdseed = 1 << 1,
	avx = 1 << 2,
	f16c = 1 << 3,
	fma = 1 << 4,
	avx2 = 1 << 5,
	aes = 1 << 6,
	vaes = 1 << 7,
	sha = 1 << 8,
	avx512f = 1 << 9,
	avx512cd = 1 << 10,
	avx512bw = 1 << 11,
	avx512dq = 1 << 12,
	avx512vl = 1 << 13,
};

enum class cpuclass
{
	unknown,
	ivb,
};

struct cpu
{
	cpuflags flags;
	cpu()
		: flags()
	{
	}
	cpu(int cpuclass_i)
	{
		switch (static_cast<cpuclass>(cpuclass_i))
		{
		case cpuclass::ivb:
			break;
		default:
			flags = cpuflags::NONE;
		}
	}
	void set_flag(cpuflags f)
	{
		flags = static_cast<cpuflags>(flags | f);
	}
};

cpu query_x86_capabilities();

}