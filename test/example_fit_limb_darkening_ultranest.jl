using OITOOLS, UltraNest
oifitsfile = "../demos/data/AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.

visfunc = visibility_ud
init_param = [8.0]
weights=[1.0,0,0]
nparams=length(init_param)
lbounds = 0.0
hbounds = 10.0

# setup priors
function prior_transform(u::AbstractVector{<:Real})
        Δx = hbounds - lbounds
        u .* Δx .+ lbounds
end

prior_transform_vectorized = let trafo = prior_transform
    (U::AbstractMatrix{<:Real}) -> reduce(vcat, (u -> trafo(u)').(eachrow(U)))
end

loglikelihood=param->model_to_chi2(data, visfunc, param, weights=weights);

loglikelihood_vectorized = let loglikelihood = loglikelihood
    # UltraNest has variate in rows:
    (X::AbstractMatrix{<:Real}) -> loglikelihood.(eachrow(X))
end

paramnames = ["UD"]
smplr = ultranest.ReactiveNestedSampler(paramnames, loglikelihood_vectorized, transform = prior_transform_vectorized, vectorized = true)
result = smplr.run(min_num_live_points = 1000, cluster_num_live_points = 400)
