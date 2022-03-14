# vp

vapoursynth plugin

## Dependencies

- [vapoursynth](https://github.com/vapoursynth/vapoursynth)
- [fftw3](https://github.com/fftw/fftw3)
- [xxhash](https://github.com/cyan4973/xxhash)
- [lz4](https://github.com/lz4/lz4)

## Building

```
cmake -S <path-to-source> -B <path-to-build> -DCMAKE_BUILD_TYPE=Release -DINCLUDE_PATH= -DLIBRARY_PATH=
cmake --build <path-to-build>
```
