using Pkg
if in("NLopt",keys(Pkg.installed()))
    @eval using NLopt

function fit_model_v2(data::OIdata, visfunc, init_param::Array{Float64,1}; fitter=:LN_NELDERMEAD, verbose = true, calculate_vis = true)
    nparams=length(init_param)
    nv2 = length(data.v2);
    chisq=(param,g)->norm( (abs2.(visfunc(param,data.v2_baseline))-data.v2)./data.v2_err)^2/nv2;
    opt = Opt(fitter, nparams);
    min_objective!(opt, chisq)
    xtol_rel!(opt,1e-5)
    bounds = zeros(size(init_param)); # note: this enforces positivity on radius and coefficients
    # maybe not desirable for quadratic law ?
    lower_bounds!(opt, bounds);
    (minf,minx,ret) = optimize(opt, init_param);
    if verbose == true
        println("Chi2: $minf \t parameters:$minx \t \t $ret")
    end
    cvis_model = [];
    if calculate_vis == true
        cvis_model = visfunc(minx,sqrt.(data.uv[1,:].^2+data.uv[2,:].^2))
    end
    return (minf,minx,cvis_model)
end

using Statistics

function bootstrap_v2_fit(nbootstraps, data::OIdata, visfunc, init_param::Array{Float64,1}; fitter=:LN_NELDERMEAD)
    npars=length(init_param)
    println("Finding mode...")
    f_chi2, params_mode, cvis_model = fit_model_v2(data, visfunc, init_param);#diameter, ld1, ld2 coeffs
    params = zeros(Float64, npars , nbootstraps)
    println("Now boostraping to estimate errors...")
    for k=1:nbootstraps
    if (k% Int.(ceil.(nbootstraps/100)) == 0)
        println("Boostrap $(k) out of $(nbootstraps)");
    end
     f_chi2, paropt, ~ = fit_model_v2(resample_v2_data(data), visfunc, params_mode, fitter= fitter, verbose = false, calculate_vis = false);#diameter, ld1, ld2 coeffs
     params[:, k]= paropt;
    end

    for i=1:npars
    fig  = figure("Histogram $(i)",figsize=(5,5));
    clf();
    plt[:hist](params[i,:],50);
    title("Bootstrap for parameter $(i)");
    xlabel("Value of parameter $(i)");
    ylabel("Boostrap samples");
    end
    params_mean = mean(params, dims=2)
    params_err = std(params, dims=2, corrected=false)
    println("Mode (maximum likelihood from original data): $(params_mode)")
    println("Boostrap mean: $(params_mean)");
    println("Bootstrap standard deviation: $(params_err)");
    return params_mode, params_mean,params_err
end

else
    @warn("NLopt not installed: related functions will not be loaded")
end

#

if in("MultiNest",keys(Pkg.installed()))
    @eval using MultiNest

function fit_model_v2_nest(data::OIdata, visfunc, bounds::Array{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},1})
npar = size(bounds,1);
v2_context = [npar*1.0; minimum.(bounds); maximum.(bounds); data.v2_baseline;data.v2;data.v2_err]; # should be const ?
minx = NaN*zeros(npar);
function loglike_v2(cube::Vector{Cdouble}, context::Array{Float64, 1})
   npar = Int.(context[1])
   nv2 = Int.((length(context)-2*npar-1)/3)
   sampledparam = cube[1:npar].*(context[2+npar:1+2*npar]-context[2:1+npar])+context[2:1+npar]
   cube[1:npar] = sampledparam[:]
   #v2_baseline = context[2:nv2+1]; v2 = context[nv2+2:2*nv2+1]; v2_err = context[2*nv2+2:end];
   return -0.5 * nv2 * log(2*pi) - 0.5*sum(log.(context[2*nv2+2*npar+2:end]))-0.5 * norm( (abs2.(visfunc(cube,context[2*npar+2:nv2+2*npar+1]))-context[nv2+2*npar+2:2*nv2+2*npar+1])./context[2*nv2+2*npar+2:end])^2;
end

function dumper(physlive::Array{Cdouble, 2},posterior::Array{Cdouble, 2}, paramconstr::Array{Cdouble, 2},
                                                    maxloglike::Cdouble,logz::Cdouble, inslogz::Cdouble, logzerr::Cdouble, context::Array{Cdouble,1})
                                                    println(size(posterior, 1), " samples and logZ= ",logz, " +/- ", logzerr);
                                                    println("Mean parameters: ", paramconstr[:,1], " +/- ", paramconstr[:,2]);
                                                    println("Best parameters: ", paramconstr[:,3]);
                                                    println("MAP parameters:  ", paramconstr[:,4]);
                                                    minx = paramconstr[:,1];
end

nest = nested(loglike_v2, npar, "chains/eggbox_context_jl-",
ins = true,  # do Nested Importance Sampling?
mmodal = false, # do mode separation?
ceff = false,  # run in constant efficiency mode?
nlive = 1000,  # number of live points
efr = 0.8,  # set the required efficiency
tol = .1,  # tol, defines the stopping criteria
npar = npar,   # total no. of parameters including free & derived parameters
nclspar = npar, # no. of parameters to do mode separation on
updint = 100, # after how many iterations feedback is required & the output files should be updated
# note: posterior files are updated & dumper routine is called after every updInt*10 iterations
ztol = -1E90, # all the modes with logZ < Ztol are ignored
maxmodes = 10, # expected max no. of modes (used only for memory allocation)
wrap = false, # which parameters to have periodic boundary conditions?
seed = -1, # random no. generator seed, if < 0 then take the seed from system clock
fb = false, # need feedback on standard output?
resume = false, # resume from a previous job?
outfile = true, # write output files?
initmpi = true, # initialize MPI routines?, relevant only if compiling with MPI
logzero = nextfloat(-Inf), # points with loglike < logZero will be ignored by MultiNest
maxiter = 0,
dumper = dumper,
context = v2_context);
# run MultiNest
run(nest);
minf =norm( (abs2.(visfunc(minx,data.v2_baseline))-data.v2)./data.v2_err)^2;
cvis_model = visfunc(minx,sqrt.(data.uv[1,:].^2+data.uv[2,:].^2));
return (minf,minx,cvis_model)
end

else
   @warn("MultiNest is not installed: related functions will not be loaded")
end



#
# BOOSTRAPING
#


function resample_v2_data(data_input)
    # returns a new sample
data_out = deepcopy(data_input);
indx_resampling = Int.(ceil.(length(data_input.v2)*rand(length(data_input.v2))));
data_out.v2 = data_input.v2[indx_resampling]
data_out.v2_err = data_input.v2_err[indx_resampling]
data_out.v2_baseline = data_input.v2_baseline[indx_resampling]
data_out.nv2 = length(data_input.v2)
data_out.v2_mjd  = data_input.v2_mjd[indx_resampling]
data_out.v2_lam  = data_input.v2_lam[indx_resampling]
data_out.v2_dlam = data_input.v2_dlam[indx_resampling]
data_out.v2_flag = data_input.v2_flag[indx_resampling]
return data_out
end
