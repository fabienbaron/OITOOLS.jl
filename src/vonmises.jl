using SpecialFunctions, LinearAlgebra

function evaluate_polynomial(p, z)
  res = 0
  for i=1:length(p)
    res += p[i]*z^(i-1)
  end
  return res
end

function log1pexp(x::Real)
  x < 18.0 ? log1p(exp(x)) : x < 33.3 ? x + exp(-x) : oftype(exp(-x), x)
end



function logbesselI0(z)
  if  z < 7.75
    p=[1.00000000000000000e+00, 2.49999999999999909e-01,
    2.77777777777782257e-02, 1.73611111111023792e-03,
    6.94444444453352521e-05, 1.92901234513219920e-06,
    3.93675991102510739e-08, 6.15118672704439289e-10,
    7.59407002058973446e-12, 7.59389793369836367e-14,
    6.27767773636292611e-16, 4.34709704153272287e-18,
    2.63417742690109154e-20, 1.13943037744822825e-22,
    9.07926920085624812e-25]
    return log1pexp(2 * log.(z) - log(4.0) + log.(evaluate_polynomial(p, 0.25 * z.^2)))

  elseif z < 500
    p = [ 3.98942280401425088e-01,  4.98677850604961985e-02,
    2.80506233928312623e-02,  2.92211225166047873e-02,
    4.44207299493659561e-02,  1.30970574605856719e-01,
    -3.35052280231727022e+00, 2.33025711583514727e+02,
    -1.13366350697172355e+04, 4.24057674317867331e+05,
    -1.23157028595698731e+07, 2.80231938155267516e+08,
    -5.01883999713777929e+09, 7.08029243015109113e+10,
    -7.84261082124811106e+11, 6.76825737854096565e+12,
    -4.49034849696138065e+13, 2.24155239966958995e+14,
    -8.13426467865659318e+14, 2.02391097391687777e+15,
    -3.08675715295370878e+15, 2.17587543863819074e+15]
    return z + log(evaluate_polynomial(p,1.0/z)) - 0.5 * log.(z);
  else # for z>500
    p = [3.98942280401432905e-01, 4.98677850491434560e-02, 2.80506308916506102e-02, 2.92179096853915176e-02, 4.53371208762579442e-02];
    return z + log(evaluate_polynomial(p,1.0/z)) - 0.5 * log.(z);
  end
end

function logbesselI1(z) # for z>500
  p = [3.989422804014314820e-01, -1.496033551467584157e-01,-4.675105322571775911e-02, -4.090421597376992892e-02,-5.843630344778927582e-02];
  return z + log(p[1]+p[2]/z+p[3]/z^2+p[4]/z^3+p[5]/z^4) - 0.5 * log(z);
end


function AAA(k)
  # The problem is that besseli.(1,k) -> infty when k is high
  # k > 700 will overfloat, but this corresponds to only σ <~ 1 degree
  # So this is a BIG issue
  # besseli.(1,k)./besseli.(0,k) is bounded by [sqrt(1+1/k^2)-1/k , sqrt(1+1/(2*k)^2)-1/(2*k)]  for large k
  # and nearly equal to the upper bound so we use this approximation
  overf = findall(k.>700)
  ok = findall(k.<=700)
  res = zeros(size(k))
  res[ok] = besseli.(1,k[ok])./besseli.(0,k[ok]);
  #res[overf] = sqrt.(1 .+ 1.0 ./ (2*k[overf]).^2) - 1.0 ./(2*k[overf]);
  res[overf] = exp.(logbesselI1.(k[overf])-logbesselI0.(k[overf])) # <- this is another possibility
  return res
end

function gaussianwrapped_to_vonmises(σ)
  # Input: σ, a vector of the standard deviations of a Gaussian wrapped distributions
  # σ should be in radian
  # Output: k, a vector of the corresponding von Mises concentration parameters
  # Method: simple Newton-Raphson algorithm
  # This will not accept values 0 <= phi or var_phi > 180 degrees
  # Note: don't confuse circular variance of an angle with its variance, which was one of PAINTER's mistakes
  #f(k)= 1 .- AAA(k) .- (1 .- exp.(-0.5*σ.^2));
  f(k)= exp.(-0.5*σ.^2) .- AAA(k);
  g(k) = -1 .+ AAA(k)./k + AAA(k).^2;
  k = ones(size(σ)); #init value that works ok in practice
  for it=1:10  # use norm f(k)./g(k) ?
    k .-= f(k)./g(k)
  end
  return k
end

function gaussianwrapped_to_vonmises_fast(σ)
    #Best & Fisher approximation - usually better than fiddling with the exact expression 
    x = exp.(-0.5*σ.^2)
    if x<0.53
        return 2*x+x^3+5*x^5/6
    elseif x < 0.85
        return -0.4+1.39x+0.43/(1-x)
    else
        return 1/(x^3-4x^2+3x)
    end
end
