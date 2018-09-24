using NFFT

struct container
fourier_transform::NFFTPlan
end

function crash(x::Array{Complex{Float64},2}, cont::container)
result = nfft(cont.fourier_transform, x);
return result
end

samples = rand(Float64, 2, 1000)
cont = container(NFFTPlan(samples, (64,64)))
x = rand(Complex{Float64}, 64,64)
crash(x, cont)
