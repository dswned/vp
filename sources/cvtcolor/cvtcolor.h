#include "config.h"

#define LAB_CBRT_TAB_SIZE 1024

#ifdef __cplusplus
#include <type_traits>
typedef void(*cvt_t)(size_t, void*, void*, void*, void*, void*, void*);
typedef std::remove_pointer_t<cvt_t> cvt_f;
extern "C" {
	cvt_f lab2yuv_f32_ivb, yuv2lab_f32_ivb,
		lab2yuv_i16_avx, yuv2lab_i16_avx,
		lab2yuv_f32_avx512;
	void init_lab_cbrt_tab();
	extern const float lab_yuv_i16_data[];
	extern float lab_cbrt_tab[];
}
constexpr float
k_sa = +2.54485449026360e+00f,
k_sb = +9.27166736516954e-01f,
k_sai = +3.92949775252737e-01f,
k_sbi = -1.07855465539743e+00f,
k_116_100 = 116 / 100.f,
k_16_100 = 16 / 100.f,
k_100_116 = 100 / 116.f,
k_16_116 = 16 / 116.f,
k_512_65535 = 512 / 65535.f,
k_1_65535 = 1 / 65535.f,
k_1_65536 = 1 / 65536.f,
k_65535 = 65535.f,
k_65536 = 65536.f,
k_32768 = 32768.f,
k_half = .5f,
k_2700_24389 = 2700 / 24389.f,
k_008 = 0.08f,
k_6_29 = 6 / 29.f,
k_sai_3132_24389 = k_sai * 3132 / 24389,
k_sbi_3132_24389 = k_sbi * 3132 / 24389;
#endif