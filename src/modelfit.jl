# Model fitting
# TODO: - add differential visibilities
#       - merge orbit code
#       - expression evaluation for custom models
#       - more parameter distributions than uniform (e.g. Gausian, von Mises or lognormal)
#       - black body spectral law
#       - spectral lines (Gaussian, Voigt)
#       - custom r/μ profiles
using Statistics, LinearAlgebra, Parameters, PyCall, UltraNest, LsqFit, NLopt, Printf

@with_kw mutable struct OIparam
           name::String = "" # optional name of the compoment (e.g. "primary", "central source")
           val::Float64 = 0
           minval::Float64 = val
           maxval::Float64 = val
           step::Float64 = 0.01
           free::Bool = true
end

function pos_fixed(pos_params::Array{OIparam,1})
    return (pos_params[1].val, pos_params[2].val)
end

function spectrum_powerlaw(spectrum_params::Array{OIparam,1}, data::OIdata )
    # (λ/λ0)^d  where λ0=spectrum[2] and d=spectrum[3]
    return spectrum_params[1].val.*(data.uv_lam/spectrum_params[2].val).^spectrum_params[3].val
end


function spectrum_powerlaw(spectrum_params::Array{OIparam,1}, λ::Array{Float64,1})
    # (λ/λ0)^d  where λ0=spectrum[2] and d=spectrum[3]
    return spectrum_params[1].val.*(λ/spectrum_params[2].val).^spectrum_params[3].val
end

function spectrum_gray(spectrum_params::Array{OIparam,1}, λ::Array{Float64,1} )
    return spectrum_params[1].val #*ones(Float64,length(data.uv_lam))
end


function spectrum_gray(spectrum_params::Array{OIparam,1}, data::OIdata )
    return spectrum_params[1].val #*ones(Float64,length(data.uv_lam))
end




@with_kw mutable struct OIcomponent
           type::String  # Type of component ("UD","LDLIN", "Ring")
           name::String = "Component1" # optional name of the compoment (e.g. "primary", "central source")
           vis_function # compute visibility for this component
           vis_params::Array{OIparam,1} = [] # visibility function parameters
           pos_function = pos_fixed
           pos_params::Array{OIparam,1} = [OIparam(name="ra", val=0.0,free=false), OIparam(name="dec", val=0.0, free=false)]  # positional parameters
           spectrum_function = spectrum_gray
           spectrum_params::Array{OIparam,1} = [OIparam(name="flux", val=1.0, free=false, minval=0.0, maxval=1.0)] # spectral law parameters for spectral law
end

function Base.display(comp::OIcomponent) # extend display for Components
                pos_type ="Fixed location"
                spec_type = "Unknown"
                if comp.spectrum_function == spectrum_gray
                    spec_type = "Spectrally Gray"
                elseif comp.spectrum_function == spectrum_powerlaw
                    spec_type = "Polychromatic (power law)"
                end
                println();
                @printf("Component: %s -- Model: %s - %s - %s \n", comp.name, comp.type, pos_type, spec_type);
                println();
                println("Visibility Parameters\t\tCurrent Value\tMinimum\t\tMaximum\t\tStepsize\tFree");
                println("---------------------------------------------------------------------");
                color=:black
                for i=1:length(comp.vis_params)
                    marksym=Char(0x02716)
                    color=:black
                    if comp.vis_params[i].free==true
                        color=:red
                        marksym=Char(0x02714)
                    end
                    printstyled(@sprintf("%16s\t\t%f\t%f\t%f\t%f\t%s\n", comp.vis_params[i].name, comp.vis_params[i].val, comp.vis_params[i].minval, comp.vis_params[i].maxval,comp.vis_params[i].step,marksym),color=color)
                end
                println("---------------------------------------------------------------------");
                println("Position Parameters\t\tCurrent Value\tMinimum\t\tMaximum\t\tStepsize\tFree");
                println("---------------------------------------------------------------------");
                for i=1:length(comp.pos_params)
                    marksym=Char(0x02716)
                    color=:black
                    if comp.pos_params[i].free==true
                        color=:red
                        marksym=Char(0x02714)
                    end
                    printstyled(@sprintf("%16s\t\t%f\t%f\t%f\t%f\t%s\n", comp.pos_params[i].name, comp.pos_params[i].val, comp.pos_params[i].minval, comp.pos_params[i].maxval,comp.pos_params[i].step,marksym),color=color)
                end
                println("---------------------------------------------------------------------");
                println("Spectral Parameters\t\tCurrent Value\tMinimum\t\tMaximum\t\tStepsize\tFree");
                println("---------------------------------------------------------------------");
                for i=1:length(comp.spectrum_params)
                    marksym=Char(0x02716)
                    color=:black
                    if comp.spectrum_params[i].free==true
                        color=:red
                        marksym=Char(0x02714)
                    end
                    printstyled(@sprintf("%16s\t\t%f\t%f\t%f\t%f\t%s\n", comp.spectrum_params[i].name, comp.spectrum_params[i].val, comp.spectrum_params[i].minval, comp.spectrum_params[i].maxval,comp.spectrum_params[i].step,marksym),color=color)
                end

end


@with_kw mutable struct OImodel
    components::Array{OIcomponent,1}
    param_map
end

function Base.display(model::OImodel)
    for i=1:length(model.components)
        display(model.components[i])
    end
    println();
    println("Parameter map:");
    display(model.param_map)
    println("Total number of free parameters: $(length(model.param_map))");
    # Sanity check here:
    # minval >= maxval for free parameters
    # or val < minval or val > maxval
    for n=1:length(model.param_map)
        i,j,k = model.param_map[n]
                if j==1
            if (model.components[i].vis_params[k].minval >= model.components[i].vis_params[k].maxval)
                @warn("Wrong bounds values for $(model.components[i].name) $(model.components[i].vis_params[k].name)")
            end
            if ((model.components[i].vis_params[k].minval > model.components[i].vis_params[k].val) || (model.components[i].vis_params[k].maxval < model.components[i].vis_params[k].val))
                @warn("Value not within bounds for $(model.components[i].name) $(model.components[i].vis_params[k].name)")
            end

        elseif j==2
            if (model.components[i].pos_params[k].minval >= model.components[i].pos_params[k].maxval)
                @warn("Wrong bounds values for $(model.components[i].name) $(model.components[i].pos_params[k].name)")
            end
            if ((model.components[i].pos_params[k].minval > model.components[i].pos_params[k].val) || (model.components[i].pos_params[k].maxval < model.components[i].pos_params[k].val))
                @warn("Value not within bounds for $(model.components[i].name) $(model.components[i].vis_params[k].name)")
            end
        elseif j==3
            if (model.components[i].spectrum_params[k].minval >= model.components[i].spectrum_params[k].maxval)
                @warn("Wrong bounds values for $(model.components[i].name) $(model.components[i].spectrum_params[k].name)")
            end
            if ((model.components[i].spectrum_params[k].minval > model.components[i].spectrum_params[k].val) || (model.components[i].spectrum_params[k].maxval < model.components[i].spectrum_params[k].val))
                @warn("Value not within bounds for $(model.components[i].name) $(model.components[i].vis_params[k].name)")
            end

        end
    end
end


function create_component(;type::String=[], name::String="") # this autofill defaults for arguments
if type=="ud"
    model = OIcomponent(type="ud", name=name,
                   vis_function=visibility_ud,
                   vis_params= [OIparam(name="diameter", val=1.0, minval=0.0, maxval = 40.0)],
                   pos_function = pos_fixed,
                   pos_params = [OIparam(name="ra", val=0.0, free=false), OIparam(name="dec", val=0.0, free=false)],  # positional parameters
                   spectrum_function = spectrum_gray,
                   spectrum_params = [OIparam(name="flux", val=1.0, minval=0.0, maxval=1.0, free=false)])
    return model
elseif type=="Gaussian"
    model = OIcomponent(type="Gaussian", name=name,
                       vis_function=visibility_Gaussian,
                       vis_params= [OIparam(name="FWHM", val=1.0, minval=0.0, maxval = 400.0, free=true), OIparam(name="Inclination", val=0.0, minval=-180.0, maxval = 180.0, free=false), OIparam(name="Position angle", val=0.0, minval=-90.0, maxval = 90.0, free=false)  ],
                       pos_function = pos_fixed,
                       pos_params = [OIparam(name="ra", val=0.0, free=false), OIparam(name="dec", val=0.0, free=false)],  # positional parameters
                       spectrum_function = spectrum_gray,
                       spectrum_params = [OIparam(name="flux", val=1.0, minval=0.0, maxval=1.0, free=false)])
    return model
elseif type=="ud-polychromatic"
        model = OIcomponent(type="ud", name=name,
                       vis_function=visibility_ud,
                       vis_params= [OIparam(name="diameter", val=1.0, minval=0.0, maxval = 40.0)],
                       pos_function = pos_fixed,
                       pos_params = [OIparam(name="ra", val=0.0, free=false), OIparam(name="dec", val=0.0, free=false)],  # positional parameters
                       spectrum_function = spectrum_powerlaw,
                       spectrum_params = [OIparam(name="flux", val=0.5, minval=0.0, maxval=1.0, free=false), OIparam(name="λ0", val=1.6e-6, free=false), OIparam(name="Spectral index", val=-4.0, minval=-5.0, maxval=2.0, free=false)])
        return model
elseif type == "ldlin"
    model = OIcomponent(type="ldlin", name=name,
                   vis_function=visibility_ldlin,
                   vis_params= [OIparam(name="diameter", val=1.0, minval=0.0, maxval = 40.0), OIparam(name="ld1", val=0.2, minval=0.0, maxval=1.0)],
                   pos_function = pos_fixed,
                   pos_params = [OIparam(name="ra", val=0.0, free=false), OIparam(name="dec", val=0.0, free=false)],  # positional parameters
                   spectrum_function = spectrum_gray,
                   spectrum_params = [OIparam(name="flux", val=1.0, minval=0.0, maxval=1.0, free=false)])

    return model
elseif type == "ldquad"
        model = OIcomponent(type="ldquad", name=name,
                       vis_function=visibility_ldquad,
                       vis_params= [OIparam(name="diameter", val=1.0, minval=0.0, maxval = 40.0), OIparam(name="ld1", val=0.2, minval=0.0, maxval=1.0), OIparam(name="ld2", val=0.2, minval=0.0, maxval=1.0)],
                       pos_function = pos_fixed,
                       pos_params = [OIparam(name="ra", val=0.0), OIparam(name="dec", val=0.0)],  # positional parameters
                       spectrum_function = spectrum_gray,
                       spectrum_params = [OIparam(name="flux", val=1.0, minval=0.0, maxval=1.0, free=false)])
    return model
elseif type == "ldpow"
    model = OIcomponent(type="ldpow", name=name,
                   vis_function=visibility_ldpow,
                   vis_params= [OIparam(name="diameter", val=1.0, minval=0.0, maxval = 40.0), OIparam(name="ld1", val=0.2, minval=0.0, maxval=1.0)],
                   pos_function = pos_fixed,
                   pos_params = [OIparam(name="ra", val=0.0, free=false), OIparam(name="dec", val=0.0, free=false)],  # positional parameters
                   spectrum_function = spectrum_gray,
                   spectrum_params = [OIparam(name="flux", val=1.0, minval=0.0, maxval=1.0, free=false)])
    return model
elseif type == "ldsqrt"
    model = OIcomponent(type="ldsqrt", name=name,
                   vis_function=visibility_ldsquareroot,
                   vis_params= [OIparam(name="diameter", val=1.0, minval=0.0, maxval = 40.0), OIparam(name="ld1", val=0.2, minval=0.0, maxval=1.0), OIparam(name="ld2", val=0.2, minval=0.0, maxval=1.0)],
                   pos_function = pos_fixed,
                   pos_params = [OIparam(name="ra", val=0.0, free=false), OIparam(name="dec", val=0.0, free=false)],  # positional parameters
                   spectrum_function = spectrum_gray,
                   spectrum_params = [OIparam(name="flux", val=1.0, minval=0.0, maxval=1.0, free=false)])
    return model
elseif type == "ring"
    model = OIcomponent(type="ring", name=name,
                   vis_function=visibility_GaussianLorentzian_ring_az,
                   vis_params= [OIparam(name="Ring radius", val=1.0, minval=0.0, maxval = 40.0), OIparam(name="Position Angle", minval=0.0, maxval=180.0), OIparam(name="Inclination", minval=0.0, maxval=180.0), OIparam(name="Ring FWHM/Radius", minval=0.0, maxval=1.0),
                    OIparam(name="Az1", minval=-1.0, maxval=1.0),OIparam(name="Az2", minval=-1.0, maxval=1.0),OIparam(name="Az3", minval=-1.0, maxval=1.0),OIparam(name="Az4", minval=-1.0, maxval=1.0),OIparam(name="Gauss/Lorentz ratio", val=0, minval=0, maxval=1.0)],
                   pos_function = pos_fixed,
                   pos_params = [OIparam(name="ra", val=0.0, free=false), OIparam(name="dec", val=0.0, free=false)],  # positional parameters
                   spectrum_function = spectrum_gray,
                   spectrum_params = [OIparam(name="flux", val=1.0, minval=0.0, maxval=1.0, free=false)])
    return model
elseif type == "ring-polychromatic"
    model = OIcomponent(type="ring", name=name,
                   vis_function=visibility_GaussianLorentzian_ring_az,
                   vis_params= [OIparam(name="Ring radius", val=1.0, minval=0.0, maxval = 40.0), OIparam(name="Position Angle", minval=0.0, maxval=180.0), OIparam(name="Inclination", minval=0.0, maxval=180.0), OIparam(name="Ring FWHM/Radius", minval=0.0, maxval=1.0),
                    OIparam(name="Az1", minval=-1.0, maxval=1.0),OIparam(name="Az2", minval=-1.0, maxval=1.0),OIparam(name="Az3", minval=-1.0, maxval=1.0),OIparam(name="Az4", minval=-1.0, maxval=1.0),OIparam(name="Gauss/Lorentz ratio", val=0, minval=0, maxval=1.0)],
                   pos_function = pos_fixed,
                   pos_params = [OIparam(name="ra", val=0.0, free=false), OIparam(name="dec", val=0.0, free=false)],  # positional parameters
                   spectrum_function = spectrum_powerlaw,
                   spectrum_params = [OIparam(name="flux", val=0.5, minval=0.0, maxval=1.0, free=false), OIparam(name="λ0", val=1.6e-6, free=false), OIparam(name="Spectral index", val=1.0, minval=-5.0, maxval=2.0, free=true)])
    return model
else
    @warn("Trying to call undefined component type");
end
end


function create_model(args...)
    components=OIcomponent[]
    for i=1:length(args)
        if typeof(args[i])!=OIcomponent
            @warn("Warning: one of the arguments of create_model was not a component")
        end
        push!(components, args[i])
    end
    param_map = []
    for i=1:length(components)
        for j=1:length(components[i].vis_params)
            if(components[i].vis_params[j].free)
                push!(param_map, [i,1,j])
            end
        end
        for j=1:length(components[i].pos_params)
            if(components[i].pos_params[j].free)
                push!(param_map, [i,2,j])
            end
        end
        for j=1:length(components[i].spectrum_params)
            if(components[i].spectrum_params[j].free)
                push!(param_map, [i,3,j])
            end
        end
    end
    return OImodel(components, param_map);
end



function model_to_cvis(model::OImodel, data::OIdata)
    V=zeros(Complex{Float64},size(data.uv,2))
    flux = zeros(Float64, size(data.uv,2)); # normalization
    for i=1:length(model.components)
        # Estimate the flux polychromatic behavior
        f = model.components[i].spectrum_function(model.components[i].spectrum_params, data)
        x,y = model.components[i].pos_function(model.components[i].pos_params)
        # Visibility calculation
        visparams = [model.components[i].vis_params[j].val for j=1:length(model.components[i].vis_params)]  # slow step... any way to speed this up?
        V += f.*model.components[i].vis_function(visparams, data.uv).*cis.(2*pi/206264806.2*(data.uv[1,:]*x - data.uv[2,:]*y))
        flux .+= f
    end
return V./flux
end


function model_to_cvis(model::OImodel, uv::Array{Float64,2}, λ::Array{Float64,1})
    V=zeros(Complex{Float64},size(uv,2))
    flux = zeros(Float64, size(uv,2)); # normalization
    for i=1:length(model.components)
        # Estimate the flux polychromatic behavior
        f = model.components[i].spectrum_function(model.components[i].spectrum_params, λ)
        x,y = model.components[i].pos_function(model.components[i].pos_params)
        # Visibility calculation
        visparams = [model.components[i].vis_params[j].val for j=1:length(model.components[i].vis_params)]  # slow step... any way to speed this up?
        V += f.*model.components[i].vis_function(visparams, uv).*cis.(2*pi/206264806.2*(uv[1,:]*x - uv[2,:]*y))
        flux .+= f
    end
return V./flux
end

function dispatch_params(params::AbstractVector{<:Real}, model::OImodel)
  # we should have length(params) == length(model.param_map)
  for n=1:length(model.param_map)
      i,j,k = model.param_map[n]
      if j==1
          model.components[i].vis_params[k].val = params[n]
      elseif j==2
          model.components[i].pos_params[k].val = params[n]
      elseif j==3
          model.components[i].spectrum_params[k].val = params[n]
      end
  end
end

function model_to_chi2(data::OIdata, model::OImodel, params::AbstractVector{<:Real}; chi2_weights=[1.0,1.0,1.0])
    # Dispatch params to model
    dispatch_params(params, model);
    #Compute vis
    cvis_model = model_to_cvis(model, data);
    chi2_v2 =0.0; chi2_t3amp =0.0; chi2_t3phi=0.0;
    if (data.nv2>0) && (chi2_weights[1]>0.0)
        v2_model = cvis_to_v2(cvis_model, data.indx_v2);
        chi2_v2 = sum( ((v2_model - data.v2)./data.v2_err).^2)/data.nv2;
    else
        chi2_weights[1]=0.0
    end

    if (data.nt3amp>0 || data.nt3phi>0)  && (chi2_weights[2]>0 || chi2_weights[3]>0)
        t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
        if (data.nt3amp>0) && (chi2_weights[2]>0.0)
        chi2_t3amp = sum( ((t3amp_model - data.t3amp)./data.t3amp_err).^2)/data.nt3amp;
        else
            chi2_weights[2]=0.0
        end
        if (data.nt3phi>0) && (chi2_weights[3]>0.0)
        chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi)./data.t3phi_err).^2)/data.nt3phi;
        else
            chi2_weights[3] = 0.0;
        end
    else
        chi2_weights[2] = 0.0;
        chi2_weights[3] = 0.0;
    end

    chi2 = (chi2_weights'*[chi2_v2, chi2_t3amp, chi2_t3phi])[1]/sum(chi2_weights)
end


function get_model_bounds(model::OImodel)
    # Setup bounds
    lbounds = Float64[]
    hbounds = Float64[]
    for i=1:length(model.components)
        for j=1:length(model.components[i].vis_params)
            if(model.components[i].vis_params[j].free)
                push!(lbounds, model.components[i].vis_params[j].minval)
                push!(hbounds, model.components[i].vis_params[j].maxval)
            end
        end
        for j=1:length(model.components[i].pos_params)
            if(model.components[i].pos_params[j].free)
                push!(lbounds, model.components[i].pos_params[j].minval)
                push!(hbounds, model.components[i].pos_params[j].maxval)
            end
        end
        for j=1:length(model.components[i].spectrum_params)
            if(model.components[i].spectrum_params[j].free)
                push!(lbounds, model.components[i].spectrum_params[j].minval)
                push!(hbounds, model.components[i].spectrum_params[j].maxval)
            end
        end
    end
    return lbounds, hbounds
end


function get_model_params(model::OImodel)
    # Setup bounds
    params = Float64[]
    for i=1:length(model.components)
        for j=1:length(model.components[i].vis_params)
            if(model.components[i].vis_params[j].free)
                push!(params, model.components[i].vis_params[j].val)
            end
        end
        for j=1:length(model.components[i].pos_params)
            if(model.components[i].pos_params[j].free)
                push!(params, model.components[i].pos_params[j].val)
            end
        end
        for j=1:length(model.components[i].spectrum_params)
            if(model.components[i].spectrum_params[j].free)
                push!(params, model.components[i].spectrum_params[j].val)
            end
        end
    end
    return params
end


function get_model_pnames(model::OImodel)
    param_names = String[]
    for n=1:length(model.param_map)
        i,j,k = model.param_map[n]
        if j==1
            push!(param_names, string(model.components[i].name, " - ", model.components[i].vis_params[k].name))
        elseif j==2
            push!(param_names, string(model.components[i].name, " - ", model.components[i].pos_params[k].name))
        elseif j==3
            push!(param_names, string(model.components[i].name, " - ", model.components[i].spectrum_params[k].name))
        end
    end
    return param_names
end



function fit_model_ultranest(data::OIdata, model::OImodel; lbounds = Float64[], hbounds = Float64[],
    verbose = true, calculate_vis = true, cornerplot = true, chi2_weights=[1.0,1.0,1.0], min_num_live_points = 400, cluster_num_live_points = 100)

    lbounds, hbounds = get_model_bounds(model);

    function prior_transform(u::AbstractVector{<:Real}) # To be modified to accept other distributions via distributions.jl?
            Δx = hbounds - lbounds
            u .* Δx .+ lbounds
    end

    prior_transform_vectorized = let trafo = prior_transform
        (U::AbstractMatrix{<:Real}) -> reduce(vcat, (u -> trafo(u)').(eachrow(U)))
    end

    loglikelihood=param::AbstractVector{<:Real}->-0.5*model_to_chi2(data, model, param, chi2_weights=chi2_weights);


    loglikelihood_vectorized = let loglikelihood = loglikelihood
        # UltraNest has variate in rows:
        (X::AbstractMatrix{<:Real}) -> loglikelihood.(eachrow(X))
    end

    param_names = get_model_pnames(model);

    smplr = ultranest.ReactiveNestedSampler(param_names, loglikelihood_vectorized, transform = prior_transform_vectorized, vectorized = true)
    result = smplr.run(min_num_live_points = min_num_live_points, cluster_num_live_points = cluster_num_live_points)

    minx = result["maximum_likelihood"]["point"]
    minf = model_to_chi2(data, model, minx, chi2_weights=chi2_weights);

    if verbose == true
        printstyled("Log Z: $(result["logz_single"]) Chi2: $minf \t parameters:$minx ",color=:red)
    end

    if cornerplot == true
        PyDict(pyimport("matplotlib")."rcParams")["font.size"]=[10];
        pyimport("ultranest.plot").cornerplot(result);
    end
    cvis_model = [];
    if calculate_vis == true
        dispatch_params(minx, model);
        cvis_model = model_to_cvis(model, data);
    end
    return (minf,minx,cvis_model, result);
end


function lsqmodelobs(params, model::OImodel, data::OIdata; chi2_weights=[1.0,1.0,1.0]) # LsqFit observables
# Dispatch params to model
dispatch_params(params, model);
#Compute vis
cvis_model = model_to_cvis(model, data);
# Compute observables
obs = Float64[]
if (chi2_weights[1]>0) && (data.nv2>0)
    append!(obs, cvis_to_v2(cvis_model, data.indx_v2))
end
if ((chi2_weights[2]>0) && (data.nt3amp>0))||(((chi2_weights[3]>0) && (data.t3phi>0)))

    t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
    if ((chi2_weights[2]>0) && (data.nt3amp>0))
        append!(obs, t3amp_model)
    end
    if ((chi2_weights[3]>0) && (data.nt3phi>0))
        append!(obs, t3phi_model)
    end
end
return obs
end

function fit_model_levenberg(data::OIdata, model::OImodel; verbose = true, calculate_vis = true, chi2_weights=[1.0,1.0,1.0])

    println("OITOOLS Warning: LSQFIT doesn't support mod360() on residuals");

    # Setup chi2_weights and data for weighted least squares
    wt = Float64[]
    if ((chi2_weights[1]>0) && (data.nv2>0))
        append!(wt, chi2_weights[1]./data.v2_err.^2)
    end
    if ((chi2_weights[2]>0) && (data.nt3amp>0))
        append!(wt, chi2_weights[2]./data.t3amp_err.^2)
    end
    if ((chi2_weights[3]>0) && (data.nt3phi>0))
        append!(wt, chi2_weights[3]./data.t3phi_err.^2)
    end

    ydata = Float64[]
    if (chi2_weights[1]>0) && (data.nv2>0)
        append!(ydata, data.v2)
    end
    if (chi2_weights[2]>0) && (data.nt3amp>0)
        append!(ydata, data.t3amp)
    end
    if (chi2_weights[3]>0) && (data.nt3phi>0)
        append!(ydata, data.t3phi)
    end


    lbounds, hbounds = get_model_bounds(model);
    pinit = get_model_params(model);
    m = (x,p)->lsqmodelobs(p, model, data, chi2_weights=chi2_weights);
    fit = curve_fit(m, [], ydata, wt, pinit, lower=lbounds, upper=hbounds, show_trace=true); # todo: add lower/upper bounds
    minx = fit.param
    minf = model_to_chi2(data, model, minx, chi2_weights=chi2_weights);

    if fit.converged == true
        println("Levenberg-Marquardt fit converged to chi2 = $(minf) for p=$(minx)\n")
    end
    sigma = stderror(fit)
    covar = estimate_covar(fit)
    if verbose==true
        println("Name       \t\tMinimum\t\tMaximum\t\tInit\t\tConverged ± Error");
        pnames = get_model_pnames(model);
        for i=1:length(pinit)
            @printf("%s \t%f\t%f\t%f\t%f ± %f\n", pnames[i], lbounds[i], hbounds[i], pinit[i], minx[i], sigma[i]);
        end
        println("\nCovariance matrix:");
        display("text/plain", covar);
    end
    cvis_model = [];
    if calculate_vis == true
        dispatch_params(minx, model);
        cvis_model = model_to_cvis(model, data);
    end
    return (minf, minx, cvis_model, fit)
end

function fit_model_nlopt(data::OIdata, model::OImodel; fitter=:LN_NELDERMEAD, verbose = true, calculate_vis = true, chi2_weights=[1.0,1.0,1.0])
    if verbose == true
        println("NLopt optimization with ", NLopt.algorithm_name(fitter))
    end
    pinit = get_model_params(model);
    nparams=length(pinit)
    chisq=(param,g)->model_to_chi2(data, model, param, chi2_weights=chi2_weights);
    opt = Opt(fitter, nparams);
    min_objective!(opt, chisq)
    xtol_rel!(opt,1e-5)
    lbounds, hbounds = get_model_bounds(model);
    lower_bounds!(opt, lbounds);
    upper_bounds!(opt, hbounds);
    (minf,minx,ret) = optimize(opt, pinit);
    if verbose == true
        numevals = opt.numevals # the number of function evaluations
        println("NLopt found a minimum of chi2=$minf at $minx after $numevals iterations (returned $ret)")
        println("Name       \t\tMinimum\t\tMaximum\t\tInit\t\tConverged");
        pnames = get_model_pnames(model);
        for i=1:length(pinit)
            @printf("%s \t%f\t%f\t%f\t%f\n", pnames[i], lbounds[i], hbounds[i], pinit[i], minx[i]);
        end
    end
    cvis_model = [];
    if calculate_vis == true
        dispatch_params(minx, model);
        cvis_model = model_to_cvis(model, data);
    end
    return (minf,minx,cvis_model, ret)
end

#
#
# #
# # OLD model interface -- by visibility function
# #
function visfunc_to_chi2(data::OIdata, visfunc, params::AbstractVector{<:Real}; chi2_weights=[1.0,1.0,1.0])
    cvis_model = visfunc(params, data.uv)
    chi2_v2 =0.0; chi2_t3amp =0.0; chi2_t3phi=0.0;
    if (data.nv2>0) && (chi2_weights[1]>0.0)
        v2_model = cvis_to_v2(cvis_model, data.indx_v2);
        chi2_v2 = sum( ((v2_model - data.v2)./data.v2_err).^2)/data.nv2;
    else
        chi2_weights[1]=0.0
    end

    if (data.nt3amp>0 || data.nt3phi>0)  && (chi2_weights[2]>0 || chi2_weights[3]>0)
        t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
        if (data.nt3amp>0) && (chi2_weights[2]>0.0)
        chi2_t3amp = sum( ((t3amp_model - data.t3amp)./data.t3amp_err).^2)/data.nt3amp;
        else
            chi2_weights[2]=0.0
        end
        if (data.nt3phi>0) && (chi2_weights[3]>0.0)
        chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi)./data.t3phi_err).^2)/data.nt3phi;
        else
            chi2_weights[3] = 0.0;
        end
    else
        chi2_weights[2] = 0.0;
        chi2_weights[3] = 0.0;
    end

    chi2 = (chi2_weights'*[chi2_v2, chi2_t3amp, chi2_t3phi])[1]/sum(chi2_weights)
end
#
# function model_to_chi2(data::Array{OIdata,1}, visfunc, params::Array{Float64,1}; chromatic_vector=[], chi2_weights=[1.0,1.0,1.0]) # polychromatic data case
#     nwavs = length(data)
#     chi2 = zeros(Float64, nwavs)
#     for i=1:nwavs
#         cvis_model = []
#         if chromatic_vector ==[]
#             cvis_model = visfunc(params, data[i].uv, data[i].uv_baseline)
#         else
#             cvis_model = visfunc(params, data[i].uv, data[i].uv_baseline, chromatic_vector[i]) # sometimes we need a chromatic constant term
#         end
#         chi2_v2 =0.0; chi2_t3amp =0.0; chi2_t3phi=0.0;
#         if (data[i].nv2>0) && (chi2_weights[1]>0.0)
#             v2_model = cvis_to_v2(cvis_model, data[i].indx_v2);
#             chi2_v2 = sum( ((v2_model - data[i].v2)./data[i].v2_err).^2)/data[i].nv2;
#         else
#             chi2_weights[1]=0.0
#         end
#         if (data[i].nt3amp>0 || data[i].nt3phi>0)  && (chi2_weights[2]>0 || chi2_weights[3]>0)
#             t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data[i].indx_t3_1, data[i].indx_t3_2 ,data[i].indx_t3_3);
#             if (data[i].nt3amp>0) && (chi2_weights[2]>0.0)
#             chi2_t3amp = sum( ((t3amp_model - data[i].t3amp)./data[i].t3amp_err).^2)/data[i].nt3amp;
#             else
#             chi2_weights[2]=0.0
#             end
#             if (data[i].nt3phi>0) && (chi2_weights[3]>0.0)
#             chi2_t3phi = sum( (mod360(t3phi_model - data[i].t3phi)./data[i].t3phi_err).^2)/data[i].nt3phi;
#             else
#             chi2_weights[3] = 0.0;
#             end
#         else
#             chi2_weights[2] = 0.0;
#             chi2_weights[3] = 0.0;
#         end
#         chi2[i] = (chi2_weights'*[chi2_v2, chi2_t3amp, chi2_t3phi])[1]/sum(chi2_weights)
#     end
#     return sum(chi2)/nwavs;
# end

#
# BOOSTRAPING
#
function resample_data(data_input; chi2_weights=[1.0,1.0,1.0]) # chi2_weights=0 are used to disable resampling
data_out = deepcopy(data_input);

# V2
if chi2_weights[1]>0
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
if chi2_weights[2]>0 || chi2_weights[3]>0
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

function bootstrap_fit(nbootstraps, data::OIdata, model::OImodel; fitter=:LN_NELDERMEAD, chi2_weights=[1.0,1.0,1.0])
    println("WARNING: WORK IN PROGRESS ON MODEL-FITTING - THIS FUNCTION DOES NOT WORK AT THE MOMENT");
    pinit = get_model_params(model);
    nparams=length(pinit)
    println("Finding mode...")
    f_chi2, params_mode, cvis_model = fit_model(data, model, pinit, chi2_weights=chi2_weights);#diameter, ld1, ld2 coeffs
    params = zeros(Float64, nparams, nbootstraps)
    if data.nt3phi != data.nt3amp
        @warn("This function needs updating to be used with nt3amp =/= nt3phi. ")
    end
    println("Now boostraping to estimate errors...")
    for k=1:nbootstraps
    if (k% Int.(ceil.(nbootstraps/100)) == 0)
        println("Boostrap $(k) out of $(nbootstraps)");
    end
     f_chi2, par_opt, ~ = fit_model(resample_data(data), model, params_mode, fitter= fitter, verbose = false, calculate_vis = false, chi2_weights=chi2_weights);#diameter, ld1, ld2 coeffs
     params[:, k]= par_opt;
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
