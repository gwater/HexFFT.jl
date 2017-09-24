module HexFFT

import Base: fft, ifft

export hfft2, ihfft2, fft, ifft, OffsetHexData

# HELPER FUNCTIONS
pad_zeros(arr) = hcat(arr, zeros(arr))
decimate(data, offset) = view(data, :, offset + 1:2:size(data, 2))
div2(n::Int) = n >> 1
function fold_half(data, op)
    n = size(data, 1)
    op.(view(data, 1:div2(n), :), view(data, div2(n) + 1:n, :))
end

# FORWARD TRANSFORM
function dft_nst1(data)
    out = fft(pad_zeros(data), 2)
    return decimate(out, 0), decimate(out, 1)
end

nst2(data) = repmat(fft(fold_half(data, (+)), 1), 2, 1)

coeff_inner(i, n) = exp(-2pi * im * i / n)
nst3_coeffs(n, m) = reshape(coeff_inner.(0:div2(n) - 1, n), div2(n), 1)
nst3(data) = repmat(fft(nst3_coeffs(size(data)...) .*
                        fold_half(data, (-)), 1), 2, 1)

w_inner(b, s, d, n, m) = exp(-pi * im * ((b + 2d) / 2m + (b + 2s) / n))
w(b, n, m) = [w_inner(b, s, d, n, m) for s in 0:n - 1, d in 0:m - 1]

function hfft2(data0, data1)
    g00, g01 = dft_nst1(data0)
    g10, g11 = dft_nst1(data1)
    return nst2(g00) .+ w(0, size(g10)...) .* nst2(g10),
           nst3(g01) .+ w(1, size(g11)...) .* nst3(g11)
end

# INVERSE TRANSFORM
function idft_inst1(data)
    out = 2ifft(pad_zeros(data), 2)
    return decimate(out, 0), decimate(out, 1)
end

inst2(data) = repmat(0.5ifft(fold_half(data, (+)), 1), 2, 1)

icoeff_inner(i, n) = exp(2pi * im * i / n)
inst3_coeffs(n, m) = reshape(icoeff_inner.(0:div2(n) - 1, n), div2(n), 1)
inst3(data) = repmat(0.5ifft(inst3_coeffs(size(data)...) .*
                             fold_half(data, (-)), 1), 2, 1)

iw_inner(a, r, c, n, m) = exp(pi * im * ((a + 2c) / 2m + (a + 2r) / n))
iw(a, n, m) = [iw_inner(a, r, c, n, m) for r in 0:n - 1, c in 0:m - 1]

function ihfft2(data0, data1)
    g00, g01 = idft_inst1(data0)
    g10, g11 = idft_inst1(data1)
    return 0.5(inst2(g00) .+ iw(0, size(g10)...) .* inst2(g10)),
           0.5(inst3(g01) .+ iw(1, size(g11)...) .* inst3(g11))
end

# Hex Array type
immutable OffsetHexData{T, S <: AbstractArray{T, 2}}
    odd_rows::S
    even_rows::S
end
function OffsetHexData{S}(odd_rows::S, even_rows::S)
    @assert size(odd_rows) == size(even_rows)
    OffsetHexData{eltype(odd_rows), S}(odd_rows, even_rows)
end

fft(H::OffsetHexData) = OffsetHexData(hfft2(H.odd_rows, H.even_rows)...)
ifft(H::OffsetHexData) = OffsetHexData(ihfft2(H.odd_rows, H.even_rows)...)

end # module
