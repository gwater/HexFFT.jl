# HexFFT

This package implements (Birdsong and Rummelt's 2016 algorithm)[http://ieeexplore.ieee.org/document/7532670/] for fast Fourier transforms on hexagonal lattices.

The module `HexFFT` currently exports two functions: `hfft2()` and `ihfft2()` which compute the FFT and its inverse, respectively.

The current implementation is ca. 10 times slower than `Base.fft()` for comparable rectangular grids. However significant optimization should be possible by preallocating and reusing temporary arrays.
