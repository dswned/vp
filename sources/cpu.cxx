#include "config.h"

#include <stdio.h>

#include "filter.h"
#include "cpu.h"

namespace uf {

void cpuid(int abcd[4], int a, int c = 0)
{
#if defined(_MSC_VER) && !defined(__clang__)
	__cpuidex(abcd, a, c);
#else
	__asm("cpuid" : "=a"(abcd[0]), "=b"(abcd[1]), "=c"(abcd[2]), "=d"(abcd[3]) : "a"(a), "c"(c) : );
#endif
}

uint64_t xgetbv(int c)
{
#if defined(_MSC_VER) && !defined(__clang__)
	return _xgetbv(c);
#else
	int a, d;
	__asm("xgetbv" : "=a"(a), "=d"(d) : "c"(c) : );
	return a | static_cast<uint64_t>(d) << 32;
#endif
}

const uint8_t AuthcAMDenti[] = { 'A','u','t','h','c','A','M','D','e','n','t','i' };
const uint8_t GenuntelineI[] = { 'G','e','n','u','n','t','e','l','i','n','e','I' };

enum { EAX, EBX, ECX, EDX };

cpu query_x86_capabilities()
{
	cpu cpu;
	int cpuid_00000000[4], cpuid_00000001[4], cpuid_00000007[4], cpuid_80000001[4],
		cpuid_brand_string[12], cpuid_00000007_1[4];
	cpuid(cpuid_00000000, 0);
	cpuid(cpuid_00000001, 1);
	cpuid(cpuid_00000007, 7);
	cpuid(cpuid_00000007_1, 7, 1);
	cpuid(cpuid_80000001, 0x80000001);
	int m = (cpuid_00000001[EAX] >> 4 & 0x0f) + (cpuid_00000001[EAX] >> 12 & 0xf0);
	int f = (cpuid_00000001[EAX] >> 8 & 0x0f) + (cpuid_00000001[EAX] >> 20) & 0xff;
	uint64_t xcr0 = 0;
	bool ymm = 0, zmm = 0;
	bool osxsave = bt(cpuid_00000001[ECX], 27);
	if (osxsave)
	{
		xcr0 = xgetbv(0);
		ymm = (xcr0 & 0b00000110) == 0b00000110;//xmm0_15+ymm0_15h
		zmm = (xcr0 & 0b11100110) == 0b11100110;//k0_7+zmm0_15h+zmm16_31
	}
	if (ymm)
	{
		if (bt(cpuid_00000001[ECX], 28))
			cpu.set_flag(cpuflags::avx);
	}
	if (zmm)
	{
		if (bt(cpuid_00000007[EBX], 16))
		{
			cpu.set_flag(cpuflags::avx512f);
		}
	}
	return cpu;
}
}