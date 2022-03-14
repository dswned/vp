#define mad(x,y,z) x*y+z
__kernel void conv2d(const __global float4* input, __global float4* output, __constant float4* f)
{
	int k, x, y, stride = get_global_size(0) * 8 + 8 * 2,
		off = mad24((int)get_global_id(1), stride, (int)get_global_id(0) * 8);
	float4 sum[8] = { 0.f }, t;
	input += off, output += off + stride + 8;
	for (y = 0; y < 3; y++)
	{
		for (x = 0; x < 8 * 3; x++)
		{
			t = input[x];
			for (k = 0; k < 8; k++)
			{
				sum[k] = mad(t.x, f[8 * 0 + k], sum[k]);
				sum[k] = mad(t.y, f[8 * 1 + k], sum[k]);
				sum[k] = mad(t.z, f[8 * 2 + k], sum[k]);
				sum[k] = mad(t.w, f[8 * 3 + k], sum[k]);
			}
			f += 8 * 4;
		}
		input += stride;
	}
	for (k = 0; k < 8; k++)
		output[k] = sum[k];
}
__kernel void conv2d_relu(const __global float4* input, __global float4* output,
	__constant float4* f, __constant float4* b)
{
	int k, x, y, stride = get_global_size(0) * 8 + 8 * 2,
		off = mad24((int)get_global_id(1), stride, (int)get_global_id(0) * 8);
	float4 sum[8] = { 0.f }, t;
	input += off, output += off + stride + 8;
	for (y = 0; y < 3; y++)
	{
		for (x = 0; x < 8 * 3; x++)
		{
			t = input[x];
			for (k = 0; k < 8; k++)
			{
				sum[k] = mad(t.x, f[8 * 0 + k], sum[k]);
				sum[k] = mad(t.y, f[8 * 1 + k], sum[k]);
				sum[k] = mad(t.z, f[8 * 2 + k], sum[k]);
				sum[k] = mad(t.w, f[8 * 3 + k], sum[k]);
			}
			f += 8 * 4;
		}
		input += stride;
	}
	for (k = 0; k < 8; k++)
		output[k] = max(sum[k] + b[k], 0.f);
}
__kernel void conv2d_add_output(const __global float4* input, __global float4* output, __constant float4* f)
{
	int k, x, y, stride = get_global_size(0) * 8 + 8 * 2,
		off = mad24((int)get_global_id(1), stride, (int)get_global_id(0) * 8);
	float4 sum[8] = { 0.f }, t;
	input += off, output += off + stride + 8;
	for (y = 0; y < 3; y++)
	{
		for (x = 0; x < 8 * 3; x++)
		{
			t = input[x];
			for (k = 0; k < 8; k++)
			{
				sum[k] = mad(t.x, f[8 * 0 + k], sum[k]);
				sum[k] = mad(t.y, f[8 * 1 + k], sum[k]);
				sum[k] = mad(t.z, f[8 * 2 + k], sum[k]);
				sum[k] = mad(t.w, f[8 * 3 + k], sum[k]);
			}
			f += 8 * 4;
		}
		input += stride;
	}
	for (k = 0; k < 8; k++)
		output[k] = output[k] + sum[k];
}
__kernel void conv2d_add_output_relu(const __global float4* input, __global float4* output,
	__constant float4* f, __constant float4* b)
{
	int k, x, y, stride = get_global_size(0) * 8 + 8 * 2,
		off = mad24((int)get_global_id(1), stride, (int)get_global_id(0) * 8);
	float4 sum[8] = { 0.f }, t;
	input += off, output += off + stride + 8;
	for (y = 0; y < 3; y++)
	{
		for (x = 0; x < 8 * 3; x++)
		{
			t = input[x];
			for (k = 0; k < 8; k++)
			{
				sum[k] = mad(t.x, f[8 * 0 + k], sum[k]);
				sum[k] = mad(t.y, f[8 * 1 + k], sum[k]);
				sum[k] = mad(t.z, f[8 * 2 + k], sum[k]);
				sum[k] = mad(t.w, f[8 * 3 + k], sum[k]);
			}
			f += 8 * 4;
		}
		input += stride;
	}
	for (k = 0; k < 8; k++)
		output[k] = max(output[k] + sum[k] + b[k], 0.f);
}
__kernel void conv2d_add_input(const __global float4* input, __global float4* output,
	__constant float4* f, __constant float4* b)
{
	int k, y, x, stride = get_global_size(0) * 8 + 8 * 2,
		off = mad24((int)get_global_id(1) + 1, stride, (int)get_global_id(0) * 8 + 8);
	float4 sum[8] = { 0.f }, t;
	for (y = -1; y < 2; y++)
	{
		for (x = -1 * 8; x < 8 * 2; x++)
		{
			t = input[y * stride + x + off];
			for (k = 0; k < 8; k++)
			{
				sum[k] = mad(t.x, f[8 * 0 + k], sum[k]);
				sum[k] = mad(t.y, f[8 * 1 + k], sum[k]);
				sum[k] = mad(t.z, f[8 * 2 + k], sum[k]);
				sum[k] = mad(t.w, f[8 * 3 + k], sum[k]);
			}
			f += 8 * 4;
		}
	}
	for (k = 0; k < 8; k++)
		output[k + off] = output[k + off] + sum[k] + b[k] + input[k + off];
}
__kernel void conv2d_add_input_reduce(const __global float4* input, __global float* dst, __constant float4* f)
{
	int gx = get_global_id(0), gy = get_global_id(1), gsx = get_global_size(0);
	int k, y, x, stride = gsx * 8 + 8 * 2, off = mad24(gy + 1, stride, gx * 8 + 8);
	float4 sum = 0.f;
	dst += gy * gsx + gx;
	for (y = -1; y < 2; y++)
	{
		for (x = -1 * 8; x < 8 * 2; x++)
			sum = mad(input[y * stride + x + off], *f++, sum);
	}
	for (k = 0; k < 8; k++)
		sum += input[k + off];
	*dst += sum.x + sum.y + sum.z + sum.w;
}
__kernel void conv2d_add_output_reduce(const __global float4* input, __global float* dst, __constant float4* f)
{
	int x, y, gx = get_global_id(0), gy = get_global_id(1), gsx = get_global_size(0),
		stride = gsx * 8 + 8 * 2, off = mad24(gy, stride, gx * 8);
	float4 sum = 0.f;
	input += off;
	for (y = 0; y < 3; y++)
	{
		for (x = 0; x < 8 * 3; x++)
			sum = mad(input[x], *f++, sum);
		input += stride;
	}
	off = mad24(gy, gsx, gx);
	dst[off] += sum.x + sum.y + sum.z + sum.w;
}
__kernel void pad(__global float4* p, int h, int w, int stride)
{
	int x = get_global_id(0);
	__global float4* d;
	if (x != 0)
	{
		if (x < h)
		{
			d = p + x * stride;
			for (int i = 0; i < 8; i++)
				d[i] = d[i + 8];
			d = d + stride - 8;
			for (int i = 0; i < 8; i++)
				d[i] = d[i - 8];
		}
		if (x < w)
		{
			d = p + x * 8;
			for (int i = 0; i < 8; i++)
				d[i] = d[i + stride];
			d = d + stride * h;
			for (int i = 0; i < 8; i++)
				d[i] = d[i - stride];
		}
	}
	else
	{
		d = p;
		for (int i = 0; i < 8; i++)
			d[i] = d[i + stride + 8];
		d = d + stride - 8;
		for (int i = 0; i < 8; i++)
			d[i] = d[i + stride - 8];
		d = p + stride * h;
		for (int i = 0; i < 8; i++)
			d[i] = d[i - stride + 8];
		d = d + stride - 8;
		for (int i = 0; i < 8; i++)
			d[i] = d[i - stride - 8];
	}
}
__kernel void nnedi(const __global float* input, __global float4* output, __global float* dst,
	__constant float4* e, __constant float4* s)
{
	__local float cache[19][24];
	int lx = get_local_id(0), ly = get_local_id(1), gx = get_global_id(0), gy = get_global_id(1);
	int sx = get_local_size(0), sy = get_local_size(1), gsx = get_global_size(0);
	int x, y, n = 0, k, stride = gsx + 8;
	float4 xe[8] = { 0.f }, xs[8] = { 0.f }, sum = 0.f, t;
	input += mad24((int)get_group_id(1) * sy, stride, (int)get_group_id(0) * sx);
	output += mad24(gy + 1, gsx + 2, gx + 1) * 8;
	dst += mad24(gy, gsx, gx);
	for (y = ly; y < sy + 3; y += sy)
	{
		for (x = lx; x < sx + 7; x += sx)
			cache[y][x] = input[y * stride + x];
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	for (y = 0; y < 4; y++)
	{
		t = vload4(0, &cache[ly + y][lx]);
		sum += t;
		for (k = 0; k < 8; k++)
		{
			xe[k] += t.x * e[n + 8 * 0 + k];
			xs[k] += t.x * s[n + 8 * 0 + k];
			xe[k] += t.y * e[n + 8 * 1 + k];
			xs[k] += t.y * s[n + 8 * 1 + k];
			xe[k] += t.z * e[n + 8 * 2 + k];
			xs[k] += t.z * s[n + 8 * 2 + k];
			xe[k] += t.w * e[n + 8 * 3 + k];
			xs[k] += t.w * s[n + 8 * 3 + k];
		}
		t = vload4(1, &cache[ly + y][lx]);
		sum += t;
		for (k = 0; k < 8; k++)
		{
			xe[k] += t.x * e[n + 8 * 4 + k];
			xs[k] += t.x * s[n + 8 * 4 + k];
			xe[k] += t.y * e[n + 8 * 5 + k];
			xs[k] += t.y * s[n + 8 * 5 + k];
			xe[k] += t.z * e[n + 8 * 6 + k];
			xs[k] += t.z * s[n + 8 * 6 + k];
			xe[k] += t.w * e[n + 8 * 7 + k];
			xs[k] += t.w * s[n + 8 * 7 + k];
		}
		n += 8 * 8;
	}
	float mean = (sum.x + sum.y + sum.z + sum.w) * 0.03125f;
	for (sum = 0.f, k = 0; k < 8; k++)
	{
		xs[k] *= xs[k];
		sum += xs[k];
	}
	float xs_sum = max(sum.x + sum.y + sum.z + sum.w, 1e-20f);
	for (sum = 0.f, k = 0; k < 8; k++)
	{
		t = xe[k] * xs[k] / xs_sum;
		output[k] = t, sum += t;
	}
	*dst = mean + sum.x + sum.y + sum.z + sum.w;
}