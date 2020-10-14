#
# Model fitting
# TODO: - go beyond V2+T3 fits by adding complex vis and/or diffvis
#       - currently we may feed data.uv or data.uv_baseline to visfunc to speed up calculation - however the polychromatic implementation is slightly different
#       - update bootstrap for non equal number of T3amp and T3phi
#       - update MultiNest fitter

using Statistics, LinearAlgebra, NLopt

function model_to_chi2(data::OIdata, visfunc, params::Array{Float64,1}; weights=[1.0,1.0,1.0])
    cvis_model = visfunc(params, data.uv)
    chi2_v2 =0.0; chi2_t3amp =0.0; chi2_t3phi=0.0;
    if (data.nv2>0) && (weights[1]>0.0)
        v2_model = cvis_to_v2(cvis_model, data.indx_v2);
        chi2_v2 = sum( ((v2_model - data.v2)./data.v2_err).^2)/data.nv2;
    else
        weights[1]=0.0
    end

    if (data.nt3amp>0 || data.nt3phi>0)  && (weights[2]>0 || weights[3]>0)
        t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
        if (data.nt3amp>0) && (weights[2]>0.0)
        chi2_t3amp = sum( ((t3amp_model - data.t3amp)./data.t3amp_err).^2)/data.nt3amp;
        else
            weights[2]=0.0
        end
        if (data.nt3phi>0) && (weights[3]>0.0)
        chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi)./data.t3phi_err).^2)/data.nt3phi;
        else
            weights[3] = 0.0;
        end
    else
        weights[2] = 0.0;
        weights[3] = 0.0;
    end

    chi2 = (weights'*[chi2_v2, chi2_t3amp, chi2_t3phi])[1]/sum(weights)
end

function model_to_chi2(data::Array{OIdata,1}, visfunc, params::Array{Float64,1}; chromatic_vector=[], weights=[1.0,1.0,1.0]) # polychromatic data case
    nwavs = length(data)
    chi2 = zeros(Float64, nwavs)
    for i=1:nwavs
        cvis_model = []
        if chromatic_vector ==[]
            cvis_model = visfunc(params, data[i].uv, data[i].uv_baseline)
        else
            cvis_model = visfunc(params, data[i].uv, data[i].uv_baseline, chromatic_vector[i]) # sometimes we need a chromatic constant term
        end
        
        chi2_v2 =0.0; chi2_t3amp =0.0; chi2_t3phi=0.0;
        if (data[i].nv2>0) && (weights[1]>0.0)
            v2_model = cvis_to_v2(cvis_model, data[i].indx_v2);
            chi2_v2 = sum( ((v2_model - data[i].v2)./data[i].v2_err).^2)/data[i].nv2;
        else
            weights[1]=0.0
        end
        if (data[i].nt3amp>0 || data[i].nt3phi>0)  && (weights[2]>0 || weights[3]>0)
            t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data[i].indx_t3_1, data[i].indx_t3_2 ,data[i].indx_t3_3);
            if (data[i].nt3amp>0) && (weights[2]>0.0)
            chi2_t3amp = sum( ((t3amp_model - data[i].t3amp)./data[i].t3amp_err).^2)/data[i].nt3amp;
            else
            weights[2]=0.0
            end
            if (data[i].nt3phi>0) && (weights[3]>0.0)
            chi2_t3phi = sum( (mod360(t3phi_model - data[i].t3phi)./data[i].t3phi_err).^2)/data[i].nt3phi;
            else
            weights[3] = 0.0;
            end
        else
            weights[2] = 0.0;
            weights[3] = 0.0;
        end
        chi2[i] = (weights'*[chi2_v2, chi2_t3amp, chi2_t3phi])[1]/sum(weights)
    end
    return sum(chi2)/nwavs;
end

function fit_model(data::OIdata, visfunc, init_param::Array{Float64,1}; fitter=:LN_NELDERMEAD, lbounds = [], hbounds = [], verbose = true, calculate_vis = true, weights=[1.0,1.0,1.0])
    nparams=length(init_param)
    chisq=(param,g)->model_to_chi2(data, visfunc, param, weights=weights)
    opt = Opt(fitter, nparams);
    min_objective!(opt, chisq)
    xtol_rel!(opt,1e-5)
    if lbounds == []
        lbounds,~ = init_bounds(visfunc)
    end
    if hbounds == []
        ~,hbounds = init_bounds(visfunc)
    end
    # maybe not desirable for quadratic law ?
    lower_bounds!(opt, lbounds);
    upper_bounds!(opt, hbounds);
    (minf,minx,ret) = optimize(opt, init_param);
    if verbose == true
        println("Chi2: $minf \t parameters:$minx \t \t $ret")
    end
    cvis_model = [];
    if calculate_vis == true
        cvis_model = visfunc(minx, data.uv)
    end
    return (minf,minx,cvis_model)
end

#
# BOOSTRAPING
#

function resample_data(data_input; weights=[1.0,1.0,1.0]) # weights=0 are used to disable resampling
data_out = deepcopy(data_input);

# V2
if weights[1]>0
indx_resampling = Int.(ceil.(data_input.nv2*rand(data_input.nv2)));
data_out.v2 = data_input.v2[indx_resampling]
data_out.v2_err = data_input.v2_err[indx_resampling]
data_out.v2_baseline = data_input.v2_baseline[indx_resampling]
data_out.nv2 = length(data_input.v2)
data_out.v2_mjd  = data_input.v2_mjd[indx_resampling]
data_out.v2_lam  = data_input.v2_lam[indx_resampling]
data_out.v2_dlam = data_input.v2_dlam[indx_resampling]
data_out.v2_flag = data_input.v2_flag[indx_resampling]
data_out.v2_sta_index = data_input.v2_sta_index[:,indx_resampling]
data_out.indx_v2= data_input.indx_v2[indx_resampling]
end
# T3
if weights[2]>0 || weights[3]>0
    indx_resampling = Int.(ceil.(data_input.nt3phi*rand(data_input.nt3phi))); # needs updating if nt3amp =/= nt3phi
    data_out.t3amp = data_input.t3amp[indx_resampling]
    data_out.t3amp_err = data_input.t3amp_err[indx_resampling]
    data_out.nt3amp = data_input.nt3amp
    data_out.t3phi = data_input.t3phi[indx_resampling]
    data_out.t3phi_err = data_input.t3phi_err[indx_resampling]
    data_out.nt3phi = data_input.nt3phi
    data_out.t3_baseline = data_input.t3_baseline[indx_resampling]
    data_out.t3_mjd  = data_input.t3_mjd[indx_resampling]
    data_out.t3_lam  = data_input.t3_lam[indx_resampling]
    data_out.t3_dlam = data_input.t3_dlam[indx_resampling]
    data_out.t3_flag = data_input.t3_flag[indx_resampling]
    data_out.t3_sta_index = data_input.t3_sta_index[:,indx_resampling]
    data_out.indx_t3_1= data_input.indx_t3_1[indx_resampling]
    data_out.indx_t3_2= data_input.indx_t3_2[indx_resampling]
    data_out.indx_t3_3= data_input.indx_t3_3[indx_resampling]
end
return data_out
end

function bootstrap_fit(nbootstraps, data::OIdata, visfunc, init_param::Array{Float64,1}; fitter=:LN_NELDERMEAD, weights=[1.0,1.0,1.0])
    npars=length(init_param)
    println("Finding mode...")
    f_chi2, params_mode, cvis_model = fit_model(data, visfunc, init_param, weights=weights);#diameter, ld1, ld2 coeffs
    params = zeros(Float64, npars , nbootstraps)
    if data.nt3phi != data.nt3amp
        @warn("This function needs updating to be used with nt3amp =/= nt3phi. ")
    end
    println("Now boostraping to estimate errors...")
    for k=1:nbootstraps
    if (k% Int.(ceil.(nbootstraps/100)) == 0)
        println("Boostrap $(k) out of $(nbootstraps)");
    end
     f_chi2, paropt, ~ = fit_model(resample_data(data), visfunc, params_mode, fitter= fitter, verbose = false, calculate_vis = false, weights=weights);#diameter, ld1, ld2 coeffs
     params[:, k]= paropt;
    end

    for i=1:npars
    fig  = figure("Histogram $(i)",figsize=(5,5));
    clf();
    plt.hist(params[i,:],50);
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

#
# function fit_model_v2_nest(data::OIdata, visfunc, bounds::Array{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},1})
# npar = size(bounds,1);
# v2_context = [npar*1.0; minimum.(bounds); maximum.(bounds); data.v2_baseline;data.v2;data.v2_err]; # should be const ?
# minx = NaN*zeros(npar);
# function loglike_v2(cube::Vector{Cdouble}, context::Array{Float64, 1})
#    npar = Int.(context[1])
#    nv2 = Int.((length(context)-2*npar-1)/3)
#    sampledparam = cube[1:npar].*(context[2+npar:1+2*npar]-context[2:1+npar])+context[2:1+npar]
#    cube[1:npar] = sampledparam[:]
#    #v2_baseline = context[2:nv2+1]; v2 = context[nv2+2:2*nv2+1]; v2_err = context[2*nv2+2:end];
#    return -0.5 * nv2 * log(2*pi) - 0.5*sum(log.(context[2*nv2+2*npar+2:end]))-0.5 * norm( (abs2.(visfunc(cube,context[2*npar+2:nv2+2*npar+1]))-context[nv2+2*npar+2:2*nv2+2*npar+1])./context[2*nv2+2*npar+2:end])^2;
# end
#
# function dumper(physlive::Array{Cdouble, 2},posterior::Array{Cdouble, 2}, paramconstr::Array{Cdouble, 2},
#                                                     maxloglike::Cdouble,logz::Cdouble, inslogz::Cdouble, logzerr::Cdouble, context::Array{Cdouble,1})
#                                                     println(size(posterior, 1), " samples and logZ= ",logz, " +/- ", logzerr);
#                                                     println("Mean parameters: ", paramconstr[:,1], " +/- ", paramconstr[:,2]);
#                                                     println("Best parameters: ", paramconstr[:,3]);
#                                                     println("MAP parameters:  ", paramconstr[:,4]);
#                                                     minx = paramconstr[:,1];
# end
#
# nest = nested(loglike_v2, npar, "chains/eggbox_context_jl-",
# ins = true,  # do Nested Importance Sampling?
# mmodal = false, # do mode separation?
# ceff = false,  # run in constant efficiency mode?
# nlive = 1000,  # number of live points
# efr = 0.8,  # set the required efficiency
# tol = .1,  # tol, defines the stopping criteria
# npar = npar,   # total no. of parameters including free & derived parameters
# nclspar = npar, # no. of parameters to do mode separation on
# updint = 100, # after how many iterations feedback is required & the output files should be updated
# # note: posterior files are updated & dumper routine is called after every updInt*10 iterations
# ztol = -1E90, # all the modes with logZ < Ztol are ignored
# maxmodes = 10, # expected max no. of modes (used only for memory allocation)
# wrap = false, # which parameters to have periodic boundary conditions?
# seed = -1, # random no. generator seed, if < 0 then take the seed from system clock
# fb = false, # need feedback on standard output?
# resume = false, # resume from a previous job?
# outfile = true, # write output files?
# initmpi = true, # initialize MPI routines?, relevant only if compiling with MPI
# logzero = nextfloat(-Inf), # points with loglike < logZero will be ignored by MultiNest
# maxiter = 0,
# dumper = dumper,
# context = v2_context);
# # run MultiNest
# run(nest);
# minf =norm( (abs2.(visfunc(minx,data.v2_baseline))-data.v2)./data.v2_err)^2;
# cvis_model = visfunc(minx,sqrt.(data.uv[1,:].^2+data.uv[2,:].^2));
# return (minf,minx,cvis_model)
# end
