#pragma once
namespace uf {

enum class cpuclass
{
	unknown,
	snb,
	ivb,
	hsw,
	bdw,
	skl,
	zn2
};

struct cpuflags
{
};

struct cpu
{
	cpuflags flags;
	cpuclass alias;
};

cpu query_x86_capabilities();
}