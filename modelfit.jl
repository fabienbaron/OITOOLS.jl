if Pkg.installed("NLopt") != nothing
      using NLopt

      function fit_model_v2(data::OIdata, visfunc, init_param::Array{Float64,1})
            nparams=length(init_param)
            indx= data.indx_v2
            nv2 = length(data.v2[indx]);
            r=data.v2_baseline
            chisq=(param,g)->norm( (abs2.(visfunc(param,r))-data.v2[indx])./data.v2_err[indx])^2/nv2;
            opt = Opt(:LN_NELDERMEAD, nparams);
            min_objective!(opt, chisq)
            xtol_rel!(opt,1e-5)
            bounds = zeros(size(init_param)); # note: this enforces positivity on radius and coefficients
            # maybe not desirable for quadratic law ?
            lower_bounds!(opt, bounds);
            (minf,minx,ret) = optimize(opt, init_param);
            println("Chi2: $minf \t parameters:$minx \t \t $ret")
            cvis_model = visfunc(minx,sqrt.(data.uv[1,:].^2+data.uv[2,:].^2))
            return (minf,minx,cvis_model)
      end

end
#
 if Pkg.installed("MultiNest") != nothing
       using MultiNest

 function fit_model_v2_nest(data::OIdata, visfunc, bounds::Array{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},1})
             npar = size(bounds,1);
             v2_context = [npar*1.0; minimum.(bounds); maximum.(bounds); data.v2_baseline;data.v2;data.v2_err]; # should be const ?
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
                   #	 paramConstr(4*nPar):
                   #   paramConstr(1) to paramConstr(nPar)	     	= mean values of the parameters
                   #   paramConstr(nPar+1) to paramConstr(2*nPar)    	= standard deviation of the parameters
                   #   paramConstr(nPar*2+1) to paramConstr(3*nPar)  	= best-fit (maxlike) parameters
                   #   paramConstr(nPar*4+1) to paramConstr(4*nPar)  	= MAP (maximum-a-posteriori) parameters
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
             run(nest)
             return (0.0,0.0,[])
       end
 end
