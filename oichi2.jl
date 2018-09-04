function setup_dft(uv, nx, pixsize)
scale_rad = pixsize * (pi / 180.0) / 3600000.0;
nuv = size(uv,2)
dft = zeros(Complex{Float64}, nuv, nx*nx);
xvals = -2 * pi * scale_rad * ([(mod(i-1,nx)+1) for i=1:nx*nx] - (nx+1)/2);
yvals=  -2 * pi * scale_rad * ([(div(i-1,nx)+1) for i=1:nx*nx] - (nx+1)/2);
for uu=1:nuv
    dft[uu,:] = cis.( (uv[1,uu] * xvals + uv[2,uu] * yvals));
end
return dft
end

mutable struct FFTPLAN
fftplan_uv
fftplan_v2
fftplan_t3_1
fftplan_t3_2
fftplan_t3_3
end

using NFFT
using SpecialFunctions
FFTW.set_num_threads(8);
function setup_nfft(data, nx, pixsize)
  scale_rad = pixsize * (pi / 180.0) / 3600000.0;
  fftplan_uv  = NFFTPlan(scale_rad*data.uv, (nx,nx), 4, 2.0);
  fftplan_v2   = NFFTPlan(scale_rad*data.uv[:, data.indx_v2], (nx,nx), 4, 2.0);
  fftplan_t3_1 = NFFTPlan(scale_rad*data.uv[:, data.indx_t3_1], (nx,nx), 4, 2.0);
  fftplan_t3_2 = NFFTPlan(scale_rad*data.uv[:, data.indx_t3_2], (nx,nx), 4, 2.0);
  fftplan_t3_3 = NFFTPlan(scale_rad*data.uv[:, data.indx_t3_3], (nx,nx), 4, 2.0);
  return FFTPLAN(fftplan_uv,fftplan_v2,fftplan_t3_1,fftplan_t3_2,fftplan_t3_3)
end



function setup_nfft_multiepochs(data, nx, pixsize)
  nepochs = size(data,1);
  scale_rad = pixsize * (pi / 180.0) / 3600000.0;
  fftplan_multi = Array{FFTPLAN}(nepochs);
  for i=1:nepochs
    fftplan_uv  = NFFTPlan(scale_rad*data[i].uv, (nx,nx), 4, 2.0);
    fftplan_v2   = NFFTPlan(scale_rad*data[i].uv[:, data[i].indx_v2], (nx,nx), 4, 2.0);
    fftplan_t3_1 = NFFTPlan(scale_rad*data[i].uv[:, data[i].indx_t3_1], (nx,nx), 4, 2.0);
    fftplan_t3_2 = NFFTPlan(scale_rad*data[i].uv[:, data[i].indx_t3_2], (nx,nx), 4, 2.0);
    fftplan_t3_3 = NFFTPlan(scale_rad*data[i].uv[:, data[i].indx_t3_3], (nx,nx), 4, 2.0);
    fftplan_multi[i]=FFTPLAN(fftplan_uv,fftplan_v2,fftplan_t3_1,fftplan_t3_2,fftplan_t3_3)
  end
return fftplan_multi
end



function mod360(x)
  mod.(mod.(x+180,360.)+360., 360.) - 180.
end

function cvis_to_v2(cvis, indx)
  v2_model = abs2.(cvis[indx]);
end

function cvis_to_t3(cvis, indx1, indx2, indx3)
  t3 = cvis[indx1].*cvis[indx2].*cvis[indx3];
  t3amp = abs.(t3);
  t3phi = angle.(t3)*180./pi;
  return t3, t3amp, t3phi
end


function image_to_cvis_dft(x, dft)
  cvis_model = Array{Complex{Float64}}(size(dft, 1));
  cvis_model = dft * vec(x) / sum(x);
end

function image_to_cvis_nfft(x, nfft_plan)
  flux = sum(x)
  if (ndims(x) == 1)
    nx = Int64(sqrt(length(x)))
    cvis_model = nfft(nfft_plan, reshape(x,(nx,nx)) + 0. * im) / flux;
  else
    cvis_model = nfft(nfft_plan, x + 0. * im) / flux;
  end
end

function chi2_dft_f(x, dft, data, verbose = true)
  cvis_model = image_to_cvis_dft(x, dft);
  # compute observables from all cvis
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  chi2_v2 = sum( ((v2_model - data.v2)./data.v2_err).^2);
  chi2_t3amp = sum( ((t3amp_model - data.t3amp)./data.t3amp_err).^2);
  chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi)./data.t3phi_err).^2);
  if verbose == true
    flux = sum(x);
    println("Chi2  -  Total: ", chi2_v2 + chi2_t3amp + chi2_t3phi, " V2: ", chi2_v2, " T3A: ", chi2_t3amp, " T3P: ", chi2_t3phi," Flux: ", flux)
    println("Chi2r -  Total:", (chi2_v2 + chi2_t3amp + chi2_t3phi)/(data.nv2+ data.nt3amp+ data.nt3phi), " V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi," Flux: ", flux)
  end
  return chi2_v2 + chi2_t3amp + chi2_t3phi
end

function chi2_nfft_f(x::Array{Float64,1}, fftplan::FFTPLAN, data::OIdata ) # criterion function for nfft
  cvis_model = image_to_cvis_nfft(x, fftplan.fftplan_uv);
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  chi2_v2 = vecnorm((v2_model - data.v2)./data.v2_err)^2;
  chi2_t3amp = vecnorm((t3amp_model - data.t3amp)./data.t3amp_err)^2;
  chi2_t3phi = vecnorm(mod360(t3phi_model - data.t3phi)./data.t3phi_err)^2;
  println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi," Flux: ", sum(x))
  return chi2_v2 + chi2_t3amp + chi2_t3phi
end

function chi2_vis_nfft_f(x, fftplan::FFTPLAN, data ) # criterion function plus its gradient w/r x
  cvis_model = image_to_cvis_nfft(x, fftplan.fftplan_uv);
  # compute observables from all cvis
  visamp_model = abs.(cvis_model);
  visphi_model = angle.(cvis_model)*(180./pi);
  chi2_visamp = vecnorm((visamp_model - data.visamp)./data.visamp_err)^2;
  chi2_visphi = vecnorm(mod360(visphi_model - data.visphi)./data.visphi_err)^2;
  println("VISAMP: ", chi2_visamp/data.nvisamp, " VISPHI: ", chi2_visphi/data.nvisphi)
  return chi2_visamp+chi2_visphi
end


function chi2_vis_dft_fg(x, g, dft, data ) # criterion function plus its gradient w/r x
  cvis_model = image_to_cvis_dft(x, dft);
  # compute observables from all cvis
  visamp_model = abs.(cvis_model);
  visphi_model = angle.(cvis_model)*(180./pi);
  chi2_visamp = vecnorm((visamp_model - data.visamp)./data.visamp_err)^2;
  chi2_visphi = vecnorm(mod360(visphi_model - data.visphi)./data.visphi_err)^2;
  # Original formulas
  # g_visamp = 2*sum(((visamp_model-data.visamp)./data.visamp_err.^2).*real( conj(cvis_model./visamp_model).*dft),1);
  # g_visphi = 360./pi*sum(((mod360(visphi_model-data.visphi)./data.visphi_err.^2)./abs2.(cvis_model)).*(-imag(cvis_model).*real(dft)+real(cvis_model).*imag(dft)),1);
  # Improved formulas
  g_visamp = 2.0*real(transpose(dft)*(conj(cvis_model./visamp_model).*(visamp_model-data.visamp)./data.visamp_err.^2))
  g_visphi = 360./pi*imag(transpose(dft)*((mod360(visphi_model-data.visphi)./data.visphi_err.^2)./visamp_model.^2.*conj(cvis_model)));
  g[:] = g_visamp + g_visphi;
  flux = sum(x);
  g[:] = (g - sum(x.*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
#  imdisp(x);
  println("VISAMP: ", chi2_visamp/data.nvisamp, " VISPHI: ", chi2_visphi/data.nvisphi, " Flux: ", flux)
  return chi2_visamp + chi2_visphi
end

function chi2_vis_nfft_fg(x, g, fftplan::FFTPLAN, data; mu = 1e7 ) # criterion function plus its gradient w/r x
  # Total squared variation
  # Total variation
    ϵ=1e-10
    nx = Int(sqrt(length(x)))
    y = reshape(x, nx, nx);
    yx = circshift(y,(0,-1));
    yy = circshift(y,(-1,0));
    tv_f = sum(sqrt.((y-yx).^2+(y-yy).^2+ϵ))
    tv_g =  vec((2*y-yx-yy)./(sqrt.((y-yx).^2+(y-yy).^2+ϵ)));


  cvis_model = image_to_cvis_nfft(x, fftplan.fftplan_uv);
  # compute observables from all cvis
  visamp_model = abs.(cvis_model);
  visphi_model = angle.(cvis_model)*(180./pi);
  chi2_visamp = vecnorm((visamp_model - data.visamp)./data.visamp_err)^2;
  chi2_visphi = vecnorm(mod360(visphi_model - data.visphi)./data.visphi_err)^2;
  g_visamp = 2.0*real(nfft_adjoint(fftplan.fftplan_uv,(cvis_model./visamp_model.*(visamp_model-data.visamp)./data.visamp_err.^2)));
  g_visphi = 360./pi*-imag(nfft_adjoint(fftplan.fftplan_uv,cvis_model.*((mod360(visphi_model-data.visphi)./data.visphi_err.^2)./visamp_model.^2)));
  g[:] = vec(g_visamp + g_visphi) +  mu * tv_g;
  flux = sum(x);
  g[:] = (g - sum(x.*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
  println("VISAMP: ", chi2_visamp/data.nvisamp, " VISPHI: ", chi2_visphi/data.nvisphi, " Flux: ", flux)
  return chi2_visamp + chi2_visphi +  mu * tv_f;
end

function gaussian2d(n,m,sigma)
  g2d = [exp(-((X-(m/2)).^2+(Y-n/2).^2)/(2*sigma.^2)) for X=1:m, Y=1:n]
  return g2d
end

# function create_start_image(type, nx, pixsize, param)
#   if type="gaussian"
#   gaussian2d(nx,nx,sigma)
#   elseif type="disc"
#     x = 

#   end
# return 
# end


function cdg(x) #note: this takes a 2D array
  xvals=[i for i=1:size(x,1)]
  return [sum(xvals'*x) sum(x*xvals)]/sum(x)
end

function reg_centering(x,g; verb = false) # takes a 1D array
  nx = Int(sqrt(length(x)))
  flux= sum(x)
  c = cdg(reshape(x,nx,nx))
  f = (c[1]-(nx+1)/2)^2+(c[2]-(nx+1)/2)^2
  xx = [mod(i-1,nx)+1 for i=1:nx*nx]
  yy = [div(i-1,nx)+1 for i=1:nx*nx]
  g[:] = 2*(c[1]-(nx+1)/2)*xx + 2*(c[2]-(nx+1)/2)*yy
  if verb == true
    print(" COG:", c, " REGC: ", f);
  end
    return f
end

function tvsq(x,tvsq_g; verb = false)
  # Total squared variation
  nx = Int(sqrt(length(x)))
  y = reshape(x, nx, nx);
  dx = circshift(y,(0,1))-y;
  dy = circshift(y,(1,0))-y;
  tvsq_f = vecnorm(dx)^2+vecnorm(dy)^2
  tvsq_g[:] = 2*vec((dx+dy))
  if verb == true
  print(" TVSQ:", tvsq_f);
  end
return tvsq_f
end

function tv(x,tv_g; verb = false)
ϵ=1e-10
nx = Int(sqrt(length(x)))
y = reshape(x, nx, nx);
yx = circshift(y,(0,-1));
yy = circshift(y,(-1,0));
tv_f = sum(sqrt.((y-yx).^2+(y-yy).^2+ϵ))
tv_g[:] =  vec((2*y-yx-yy)./(sqrt.((y-yx).^2+(y-yy).^2+ϵ)));
if verb == true
  print(" TV:", tv_f);
end
  return tv_f
end


function regularization(x, reg_g; printcolor = :black, regularizers=[], verb=true) # compound regularization
  reg_f = 0.0;
  for ireg =1:length(regularizers)
    temp_g = Array{Float64}(length(x))
    if regularizers[ireg][1] == "centering"
      reg_f += regularizers[ireg][2]*reg_centering(x, temp_g, verb = verb)
    elseif regularizers[ireg][1] == "tv"
      reg_f += regularizers[ireg][2]*tv(x,temp_g, verb = verb)
    elseif regularizers[ireg][1] == "tvsq"
      reg_f += regularizers[ireg][2]*tvsq(x,temp_g, verb = verb)
    elseif regularizers[ireg][1] == "EPLL"
      reg_f += regularizers[ireg][2]*EPLL_fg(x,temp_g, regularizers[ireg][3])
    end
    reg_g[:] += regularizers[ireg][2]*temp_g
  end
  if (verb==true)
   print("\n");
  end
  return reg_f
end

function crit_dft_fg(x, g, dft, data; regularizers=[], verb = true) # regularization + negloglikelihood
  chi2_f = chi2_dft_fg(x,  g, dft, data, verb = verb); # this resets g
  reg_f = regularization(x, g, regularizers=regularizers, verb = verb) ;# this adds to g
  flux = sum(x)
  g[:] = (g - sum(vec(x).*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
  return chi2_f + reg_f;
end

function crit_nfft_fg(x,g, fftplan, data;regularizers=[], verb = true)
  chi2_f = chi2_nfft_fg(x, g, fftplan, data, verb = verb);
  reg_f = regularization(x,g, regularizers=regularizers, verb = verb);
  flux = sum(x)
  g[:] = (g - sum(vec(x).*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
  return chi2_f + reg_f;
end


function chi2_dft_fg(x::Array{Float64,1}, g::Array{Float64,1}, dft::Array{Complex{Float64},2}, data::OIdata; printcolor =:black, verb=true)
  flux = sum(x);
  cvis_model = image_to_cvis_dft(x, dft);
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  # note: this is correct but slower
  chi2_v2 = vecnorm((v2_model - data.v2)./data.v2_err)^2;
  chi2_t3amp = vecnorm((t3amp_model - data.t3amp)./data.t3amp_err)^2;
  chi2_t3phi = vecnorm(mod360(t3phi_model - data.t3phi)./data.t3phi_err)^2;
  g_v2 = real(transpose(dft[data.indx_v2,:])*(4*((v2_model-data.v2)./data.v2_err.^2).*conj(cvis_model[data.indx_v2])));
  g_t3amp = real(transpose(dft[data.indx_t3_1,:])*(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_1])./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3]) ))+real(transpose(dft[data.indx_t3_2,:])*(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_2])./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3]) ))+real(transpose(dft[data.indx_t3_3,:])*(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_3])./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]) ))
  g_t3phi = 360./pi*imag(transpose(dft[data.indx_t3_1,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3].*conj(t3_model))
                         +transpose(dft[data.indx_t3_2,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3].*conj(t3_model))
                         +transpose(dft[data.indx_t3_3,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*conj(t3_model)))
  g[:] = g_v2 + g_t3amp + g_t3phi
 # g[:] = (g - sum(vec(x).*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
 if verb==true
  print_with_color(printcolor, "V2: $(chi2_v2/data.nv2) T3A: $(chi2_t3amp/data.nt3amp) T3P: $(chi2_t3phi/data.nt3phi)  Flux: $(flux) ");
 end
   return chi2_v2 + chi2_t3amp + chi2_t3phi
end

function chi2_nfft_fg(x::Array{Float64,1}, g::Array{Float64,1}, fftplan::FFTPLAN, data::OIdata; printcolor =:black,  verb = false)
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
   print_with_color(printcolor, "V2: $(chi2_v2/data.nv2) T3A: $(chi2_t3amp/data.nt3amp) T3P: $(chi2_t3phi/data.nt3phi)  Flux: $(flux) ")
 end
   return chi2_v2 + chi2_t3amp + chi2_t3phi 
end


function reconstruct(x_start::Array{Float64,1}, data, ft; printcolor = :black, verb = true, maxiter = 100, regularizers =[])
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


function crit_nfft_fg(x,g, fftplan, data; printcolor = :black, regularizers=[], verb = true)
  chi2_f = chi2_nfft_fg(x, g, fftplan, data, verb = verb);
  reg_f = regularization(x,g, regularizers=regularizers, printcolor = printcolor, verb = verb);
  flux = sum(x)
  g[:] = (g - sum(vec(x).*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
  return chi2_f + reg_f;
end

function reconstruct_multitemporal(x_start::Array{Float64,1}, data::Array{OIdata, 1}, ft; epochs_weights =[], printcolor= [], verb = true, maxiter = 100, regularizers =[])
  x_sol = []
  if typeof(ft) == Array{FFTPLAN,1}
    crit = (x,g)->crit_multitemporal_nfft_fg(x, g, ft, data, printcolor=printcolor, epochs_weights=epochs_weights, regularizers=regularizers, verb = verb)
    x_sol = OptimPack.vmlmb(crit, x_start, verb=verb, lower=0, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
  else 
    error("Sorry, polytemporal DFT methods not implemented yet");
    #crit = (x,g)->crit_dft_fg(x, g, ft, data, regularizers=regularizers, verb = verb)
    #x_sol = OptimPack.vmlmb(crit, x_start, verb=verb, lower=0, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
  end
return x_sol
end


function crit_multitemporal_nfft_fg(x::Array{Float64,1}, g::Array{Float64,1}, ft::Array{FFTPLAN,1}, data::Array{OIdata,1};printcolor= [], epochs_weights=[],regularizers=[], verb = false)
  nepochs = length(ft);
  if epochs_weights == []
    epochs_weights=ones(Float64, nepochs);
  end
  if printcolor == []
    printcolor=Array{Symbol}(nepochs);
    printcolor[:]=:black
  end
  npix = div(length(x),nepochs);
  f = 0.0;
  for i=1:nepochs # weighted sum -- should probably do the computation in parallel
    tslice = 1+(i-1)*npix:i*npix; # temporal slice
    subg = Array{Float64}(npix);
    print_with_color(printcolor[i], "Epoch $i ");
    f += epochs_weights[i]*crit_nfft_fg(x[tslice], subg, ft[i], data[i], regularizers=[], printcolor = printcolor[i], verb = verb);
    g[tslice] = epochs_weights[i]*subg
  end

  # cross temporal regularization
   if length(regularizers)>nepochs
     if (regularizers[nepochs+1][1][1] == "temporal_tvsq")  & (nepochs>1)
      y = reshape(x,(npix,nepochs))
      temporalf = sum( (y[:,2:end]-y[:,1:end-1]).^2 )
      tv_g = Array{Float64}(npix,nepochs)
      if nepochs>2
         tv_g[:,1] = 2*(y[:,1] - y[:,2])
         tv_g[:,2:end-1] = 4*y[:,2:end-1]-2*(y[:,1:end-2]+y[:,3:end])
         tv_g[:,end] = 2*(y[:,end] - y[:,end-1])
      else
         tv_g[:,1] = 2*(y[:,1]-y[:,2]);
         tv_g[:,2] = 2*(y[:,2]-y[:,1]);
      end
      f+= temporalf
      g[:] += regularizers[nepochs+1][1][2]*vec(tv_g);
      print_with_color(:yellow,"Temporal regularization: $temporalf\n")
     end
    end
   return f;
end


if Pkg.installed("Wavelets") !=nothing
  using Wavelets
  function W(mat; wavelet_bases=[WT.db1, WT.db2, WT.db3, WT.db4, WT.db5, WT.db6, WT.db7, WT.db8, WT.haar])
        nx = Int(sqrt(length(mat)))
      Wx = Array{Float64}(nx*nx,length(wavelet_bases));
        for i=1:length(wavelet_bases)
            Wx[:, i]=vec(dwt(reshape(mat,nx,nx), wavelet(wavelet_bases[i])));
        end
        return Wx;
  end
  
  function Wt(mat; wavelet_bases=[WT.db1, WT.db2, WT.db3, WT.db4, WT.db5, WT.db6, WT.db7, WT.db8, WT.haar])
    nx = size(mat,2)
    IWx = Array{Float64}(nx*nx,length(wavelet_bases));
     for i=1:length(wavelet_bases)
         IWx[:,i] = vec(idwt(reshape(mat[:,i],(nx,nx)), wavelet(wavelet_bases[i])));
      end
   return sum(IWx,2);#/length(wavelet_bases);
  end
  
  
  function regwav(x,wav_g; wavelet_bases=[WT.db1, WT.db2, WT.db3, WT.db4, WT.db5, WT.db6, WT.db7, WT.db8, WT.haar])
    wav_f = vecnorm(W(x,wavelet_bases))^2
    wav_g[:] = 2*length(wavelet_bases)*x
    return tv_f
  end
  end
  






