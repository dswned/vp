#include "uf.h"

#include <fcntl.h>
#if _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

namespace uf {

template<>
A_F16C void conv<sf16, float>::transform_n(const sf16* s, size_t n, float* d)
{
	size_t i = 0;
#if defined VP_ENABLE_INTRIN
	for (; i < n; i += 8)
	{
		if (i > n - 8)
		{
			if (i == 0)
				break;
			i = n - 8;
		}
		_mm256_storeu_ps(d + i, _mm256_cvtph_ps(_mm_castps_si128(_mm_loadu_ps(reinterpret_cast<const float*>(s + i)))));
	}
#endif
	for (; i < n; i++)
		d[i] = s[i];
}

#if _WIN32

bool file_h::open(const char* path)
{
	h = ::_open(path, _O_BINARY | _O_RDONLY);
	if (h < 0)
		return false;
	size = ::_lseeki64(h, 0, SEEK_END);
	if (size < 0 || ::_lseeki64(h, 0, SEEK_SET) < 0)
		return false;
	return true;
}

int file_h::read(void* buf, size_t nbyte)
{
	return ::_read(h, buf, (unsigned)nbyte);
}

int file_h::read(void* buf, intmax_t off, size_t nbyte)
{
	::_lseeki64(h, off, SEEK_SET);
	return ::_read(h, buf, (unsigned)nbyte);
}

file_h::~file_h()
{
	if (h >= 0)
		::_close(h);
}
#else

bool file_h::open(const char* fn)
{
	h = ::open(fn, O_RDONLY);
	if (h < 0)
		return false;
	size = ::lseek64(h, 0, SEEK_END);
	if (size < 0 || ::lseek64(h, 0, SEEK_SET) < 0)
		return false;
	return true;
}

int file_h::read(void* buf, size_t nbyte)
{
	return ::read(h, buf, (unsigned)nbyte);
}

int file_h::read(void* buf, intmax_t off, size_t nbyte)
{
	::lseek64(h, off, SEEK_SET);
	return ::read(h, buf, (unsigned)nbyte);
}

file_h::~file_h()
{
	if (h >= 0)
		::close(h);
}
#endif

}
