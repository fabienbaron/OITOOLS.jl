oifitsfiles = ["./DATA/2011Sep02.lam_And_prepped.oifits", "./DATA/2011Sep06.lam_And_prepped.oifits",
"./DATA/2011Sep10.lam_And_prepped.oifits","./DATA/2011Sep14.lam_And_prepped.oifits",
"./DATA/2011Sep19.lam_And_prepped.oifits","./DATA/2011Sep24.lam_And_prepped.oifits"];

nepochs, tepochs, data = readoifits_multiepochs(oifitsfiles);

# Setup Fourier transform (polygons -> complex visibilities)
polyflux, polyft = setup_polygon_ft(data, star_epoch_geom);
# Define epoch weights and function to minimize
epochs_weights = ones(Float64, nepochs)/nepochs;
# Define regularization weight
mu = 1e3;

crit_imaging = (x,g)->chi2_allepochs_fg(x, g, epochs_weights, polyflux, polyft, data);

f = chi2_allepochs_f(temperature_map, epochs_weights, polyflux_new, polyft_new, data);


function chi2_nfft_fg(x::Array{Float64,1}, g::Array{Float64,1}, fftplan::FFTPLAN, data::OIdata; verb = false)
  flux = sum(x);
# Likelihood
  cvis_model = image_to_cvis_nfft(x, fftplan.fftplan_uv);
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  chi2_v2 = vecnorm((v2_model - data.v2)./data.v2_err)^2;
  chi2_t3amp = vecnorm((t3amp_model - data.t3amp)./data.t3amp_err)^2;
  chi2_t3phi = vecnorm(mod360(t3phi_model - data.t3phi)./data.t3phi_err)^2;
  g_v2 = real(nfft_adjoint(fftplan.fftplan_v2, (4*((v2_model-data.v2)./data.v2_err.^2).*cvis_model[data.indx_v2])));
  g_t3amp = real(nfft_adjoint(fftplan.fftplan_t3_1, (2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*cvis_model[data.indx_t3_1]./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3]) ))) + real(nfft_adjoint(fftplan.fftplan_t3_2,(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*cvis_model[data.indx_t3_2]./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3]) ))) +real(nfft_adjoint(fftplan.fftplan_t3_3,(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*cvis_model[data.indx_t3_3]./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]) )))
  g_t3phi = -360./pi*imag(nfft_adjoint(fftplan.fftplan_t3_1, ((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*conj(cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3]).*t3_model)
                        +nfft_adjoint(fftplan.fftplan_t3_2, ((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*conj(cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3]).*t3_model)
                        +nfft_adjoint(fftplan.fftplan_t3_3, ((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*conj(cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2]).*t3_model))
  g[:] = vec(g_v2 + g_t3amp + g_t3phi)
  # g[:] = (g - sum(vec(x).*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
 if verb==true
   print("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi," Flux: ", flux)
 end
   return chi2_v2 + chi2_t3amp + chi2_t3phi
end

function chi2_nfft_allepochs_fg(x,g,epoch_weights,fftplan::FTTPLAN,data::OITdata)
    f=0;
    g[:]=0;
    npix=size(x);
    singleepoch_g=zeros(Float64,npix);
    for i=1:nepochs
        f+=epochs_weights[i]*chi2_nfft_fg(x,singleepoch_g,fftplan[i],data[i],true);
        g[:]+=epochs_weights[i]*singleepoch_g;
    end
    println("All epochs, weighted chi2: ", f, "\n");
    return f;
end

function crit_nfft_fg_allepochs(x,g, fftplan, data;regularizers=[], verb = true)
  chi2_f = chi2_nfft_fg(x, g, fftplan, data, verb = verb);
  reg_f = regularization(x,g, regularizers=regularizers, verb = verb);
  flux = sum(x)
  g[:] = (g - sum(vec(x).*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
  return chi2_f + reg_f;
end


function reconstruct(x_start::Array{Float64,1}, ft; verb = true, maxiter = 100, regularizers =[])
  x_sol = []
  if typeof(ft) == FFTPLAN
    crit = (x,g)->crit_nfft_fg(x, g, ft, data, regularizers=regularizers, verb = verb)
    x_sol = OptimPack.vmlmb(crit, x_start, verb=verb, lower=0, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
  else
    crit = (x,g)->crit_dft_fg(x, g, ft, data, regularizers=regularizers, verb = verb)
    x_sol = OptimPack.vmlmb(crit, x_start, verb=verb, lower=0, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
  end
return x_sol
end


function chi2_allepochs_fg(x, g, epochs_weights, polyflux, polyft, data)
f = 0;
g[:] = 0;
npix = size(x);
singleepoch_g = zeros(Float64, npix);
for i=1:nepochs # weighted sum -- should probably do the computation in parallel
  f += epochs_weights[i]*chi2_fg(x, singleepoch_g, polyflux[i], polyft[i], data[i], true);
  g[:] += epochs_weights[i]*singleepoch_g;
end
println("All epochs, weighted chi2: ", f, "\n");
return f;
end

function chi2_allepochs_f(x, epochs_weights, polyflux, polyft, data)
f = 0;
npix = size(x);
for i=1:nepochs # weighted sum -- should probably do the computation in parallel
  f += epochs_weights[i]*chi2_f(x, polyflux[i], polyft[i], data[i], true);
end
println("All epochs, weighted chi2: ", f, "\n");
return f;
end
