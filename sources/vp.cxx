#include "config.h"
#include <vapoursynth4.h>
#include <vector>

#if defined _MSC_VER
#define VP_CXX
#include "vp.h"
#endif

const VSPlugin* plugin = 0;
const VSAPI* vsapi = 0;
std::vector<void(*)(VSPlugin*, const VSPLUGINAPI*)> v_reg_f;

VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin* plugin, const VSPLUGINAPI* vspapi)
{
	::plugin = plugin;
	vspapi->configPlugin("xxx.xyz.vp", "vp", "VapourSynth Filter", 0, VAPOURSYNTH_API_VERSION, 0, plugin);
	for (auto reg_f : v_reg_f)
		reg_f(plugin, vspapi);
}