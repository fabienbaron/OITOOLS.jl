using OITOOLS, OptimPackNextGen, LinearAlgebra

cd("/home/baron/SOFTWARE/OITOOLS.jl/demos/")
include("../../dpi.jl/src/DPI.jl");
prior_file = "../../dpi.jl/demos/trained_networks/Glow32_AdaBelief/Glow32_YSOs.jld2"
Network = load_prior(prior_file);

oifitsfile = "./data/2019_v1295Aql.WL_SMOOTH.A.oifits"
pixsize = 0.25 #mas/pixel
nx = 32 # pixels
data = readoifits(oifitsfile, filter_bad_data=true)[1,1];
ft = setup_nfft(data, nx, pixsize);
params_start = [0.444, 0, 0, 4, 1.6e-6] #V2:12.36 T3A: 4.14 T3P: 1.61


function crit_sparco_dpi_fg(z, g, ft, data, params, Network; verb = verb, weights=weights)
    #weights=[1.0,0,1.0]; verb=true;
    x = Network.inverse(Float32.(z))
    flux = sum(x);
    nparams = length(params)
    ndata=data.nv2+data.nt3phi
    nx = size(x,1)
    α = 1.0f0;
    likelihood_g = zeros(nx*nx+nparams);
    likelihood_f = chi2_sparco_fg([params;vec(x)],  likelihood_g, ft, data, nparams, weights=weights, verb=verb)/ndata; 
    likelihood_g = (likelihood_g[nparams+1:end]  .- sum(vec(x).*likelihood_g[nparams+1:end]) / flux ) / flux/ndata; # gradient correction to take into account the non-normalized image
    #likelihood_g = likelihood_g[nparams+1:end]
    Δx = Float32.(reshape(likelihood_g, nx,nx, 1,1))
    crit_f = likelihood_f + α*norm(z)^2/(nx*nx); 
    #crit_f = likelihood_f + α*(norm(z)^2-nx*nx)^2/(2*nx*nx);
    Δz, _ = Network.backward_inv(Δx, x) # sets the data gradients
    g[:] = Δz + 2.0f0*α*z/(nx*nx);
    #g[:] = Δz + 2α*z*(norm(z)^2-nx*nx)/(nx*nx)
    return crit_f
end

function reconstruct_sparco_gray_dpi(x_start::Array{Float64,1}, params::Array{Float64,1}, data::OIdata, ft, Network; printcolor = :black, verb = true, maxiter = 100, weights=[1.0,1.0,1.0],ftol= (0,1e-8), xtol=(0,1e-8), gtol=(0,1e-8)) #grey environment
        # verb = true; maxiter = 300; weights=[1.0,0.0,1.0];ftol= (0,1e-12); xtol=(0,1e-12); gtol=(0,1e-12);
    nx = Int(sqrt(length(x_start)))
    z_start = randn(Float32,nx,nx,1,1);#Network.forward(reshape(x_start, nx, nx, 1, 1))
    crit = (z,g)->crit_sparco_dpi_fg(z, g, ft, data, params, Network, verb = false, weights=weights);
    z_opt = OptimPackNextGen.vmlmb(crit, z_start, verb=verb, maxiter=maxiter, blmvm=false, xtol = xtol, ftol = ftol, gtol=gtol);
    x_opt = Network.inverse(z_opt)
    return x_opt, z_opt
end

x_start = vec(rand(nx,nx))
x_start /= sum(x_start)
chi2_sparco_f(x_start, ft, data, params_start, weights=[1.0,0.0,1.0])
minchi2, params,ret = optimize_sparco_parameters(params_start, x_start, ft, data, weights = [1.0,0.0,1.0] )
maxiter = 200
params=copy(params_start)
x_opt, z_opt = reconstruct_sparco_gray_dpi(x_start, params, data, ft, Network, verb = true, maxiter = maxiter, weights=[1.0,0.0,1.0], xtol = (0,1e-8), ftol =(0,1e-8), gtol= (0,1e-8)) 
imdisp(vec(x_opt), pixsize = pixsize, figtitle="v1295 Aql")

chi2_sparco_f(Float64.(vec(x_opt)), ft, data, params, weights=[1.0,0.0,1.0])
