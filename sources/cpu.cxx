#include "config.h"

#include "cpu.h"
#include "uf.h"

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

cpu query_x86_capabilities()
{
	cpu cpu = {};
	return cpu;
}
}