#include "config.h"

#include <vector>
#include <vapoursynth.h>
#if defined(_MSC_VER)
#include "vp.h"
#endif

const VSAPI* vsapi = 0;
std::vector<void(*)(VSRegisterFunction, VSPlugin*)> vregf;

VS_EXTERNAL_API(void) VapourSynthPluginInit(VSConfigPlugin configFunc, VSRegisterFunction registerFunc, VSPlugin* plugin)
{
	configFunc("xxx.xyz.vp", "vp", "VapourSynth Filter", VAPOURSYNTH_API_VERSION, 1, plugin);
	for (auto regf : vregf)
		regf(registerFunc, plugin);
}