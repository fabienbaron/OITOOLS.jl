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

using NFFT
function setup_nfft(uv, nx, pixsize)
  scale_rad = pixsize * (pi / 180.0) / 3600000.0;
  nfft_plan = NFFTPlan(uv*scale_rad, (nx,nx), 4, 2.0);
  return nfft_plan
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

function chi2(x, dft, data, verbose = true)
  cvis_model = image_to_cvis_dft(x, dft);
  # compute observables from all cvis
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  chi2_v2 = sum( ((v2_model - data.v2)./data.v2_err).^2);
  chi2_t3amp = sum( ((t3amp_model - data.t3amp)./data.t3amp_err).^2);
  chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi)./data.t3phi_err).^2);
  if verbose == true
    flux = sum(x);
    println("V2: ", chi2_v2, " T3A: ", chi2_t3amp, " T3P: ", chi2_t3phi," Flux: ", flux)
    println("chi2V2: ", chi2_v2/data.nv2, " chi2T3A: ", chi2_t3amp/data.nt3amp, " chi2T3P: ", chi2_t3phi/data.nt3phi," Flux: ", flux)
  end
  return chi2_v2 + chi2_t3amp + chi2_t3phi
end

function chi2_fg(x, g, dft, data ) # criterion function plus its gradient w/r x
  #nx2 = length(x);
  cvis_model = image_to_cvis_dft(x, dft);
  # compute observables from all cvis
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  chi2_v2 = sum( ((v2_model - data.v2)./data.v2_err).^2);
  chi2_t3amp = sum( ((t3amp_model - data.t3amp)./data.t3amp_err).^2);
  chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi)./data.t3phi_err).^2);
  #tic();
  g_v2 = 4.0*sum(((v2_model-data.v2)./data.v2_err.^2).*real(conj(cvis_model[data.indx_v2]).*dft[data.indx_v2,:]),1);
  g_t3amp = 2.0*sum(((t3amp_model-data.t3amp)./data.t3amp_err.^2).*
  (real( conj(cvis_model[data.indx_t3_1]./abs.(cvis_model[data.indx_t3_1])).*dft[data.indx_t3_1,:]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3]) + real( conj(cvis_model[data.indx_t3_2]./abs.(cvis_model[data.indx_t3_2])).*dft[data.indx_t3_2,:]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3])+ real( conj(cvis_model[data.indx_t3_3]./abs.(cvis_model[data.indx_t3_3])).*dft[data.indx_t3_3,:]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2])),1);

  t3model_der = dft[data.indx_t3_1,:].*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3] + dft[data.indx_t3_2,:].*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3] + dft[data.indx_t3_3,:].*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2];
  g_t3phi =360./pi*sum(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*(-imag(t3_model).*real(t3model_der)+real(t3_model).*imag(t3model_der)),1);
  #toc();
  imdisp(x);
  g[:] = squeeze(g_v2 + g_t3amp + g_t3phi,1);
  flux = sum(x);
  g[:] = (g - sum(x.*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
  println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi," Flux: ", flux)
  return chi2_v2 + chi2_t3amp + chi2_t3phi
end

function gaussian2d(n,m,sigma)
  g2d = [exp(-((X-(m/2)).^2+(Y-n/2).^2)/(2*sigma.^2)) for X=1:m, Y=1:n]
  return g2d
end

function cdg(x) #note: this takes a 2D array
  xvals=[i for i=1:size(x,1)]
  return [sum(xvals'*x) sum(x*xvals)]/sum(x)
end

function reg_centering(x,g) # takes a 1D array
  nx = Int(sqrt(length(x)))
  flux= sum(x)
  c = cdg(reshape(x,nx,nx))
  f = (c[1]-(nx+1)/2)^2+(c[2]-(nx+1)/2)^2
  xx = [mod(i-1,nx)+1 for i=1:nx*nx]
  yy = [div(i-1,nx)+1 for i=1:nx*nx]
  g[1:nx*nx] = 2*(c[1]-(nx+1)/2)*xx + 2*(c[2]-(nx+1)/2)*yy
  return f
end

function chi2_centered_l2_fg(x, g, mu, dft, data )
flux = sum(x);
cvis_model = image_to_cvis_dft(x, dft);
# cvis_model = image_to_cvis_nfft(x, fft);
# compute observables from all cvis
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  chi2_v2 = sum( ((v2_model - data.v2)./data.v2_err).^2);
  chi2_t3amp = sum( ((t3amp_model - data.t3amp)./data.t3amp_err).^2);
  chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi)./data.t3phi_err).^2);

# centering
  rho = 1e4
  reg_der = zeros(size(x))
  reg = reg_centering(x, reg_der)
# L2
  l2 = sum(x.^2)
  l2_der = 2*x
  # note: this is correct but slower
  g_v2 = 4.0*sum(((v2_model-data.v2)./data.v2_err.^2).*real(conj(cvis_model[data.indx_v2]).*dft[data.indx_v2,:]),1);
  g_t3amp = 2.0*sum(((t3amp_model-data.t3amp)./data.t3amp_err.^2).*
  (   real( conj(cvis_model[data.indx_t3_1]./abs.(cvis_model[data.indx_t3_1])).*dft[data.indx_t3_1,:]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3])       + real( conj(cvis_model[data.indx_t3_2]./abs.(cvis_model[data.indx_t3_2])).*dft[data.indx_t3_2,:]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3])+ real( conj(cvis_model[data.indx_t3_3]./abs.(cvis_model[data.indx_t3_3])).*dft[data.indx_t3_3,:]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2])),1);
  t3model_der = dft[data.indx_t3_1,:].*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3] + dft[data.indx_t3_2,:].*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3] + dft[data.indx_t3_3,:].*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2];
  g_t3phi =360./pi*sum(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*(-imag(t3_model).*real(t3model_der)+real(t3_model).*imag(t3model_der)),1);
  imdisp(x)
  g[1:end] = vec(g_v2 + g_t3amp + g_t3phi) +  rho * reg_der + mu * l2_der;
  g[1:end] = (g - sum(vec(x).*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
  println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi," Flux: ", flux, " CENT: ", reg, " COG ", cdg(reshape(x,nx,nx)), " L2: ", mu*l2)
  return chi2_v2 + chi2_t3amp + chi2_t3phi + rho *reg + mu * l2
end



function chi2_centered_fg(x::Array{Float64,1}, g::Array{Float64,1}, dft::Array{Complex{Float64},2}, data::OIdata)
  flux = sum(x);
  cvis_model = image_to_cvis_dft(x, dft);
# cvis_model = image_to_cvis_nfft(x, fft);
  # compute observables from all cvis
  #tic();
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  chi2_v2 = sum( ((v2_model - data.v2)./data.v2_err).^2);
  chi2_t3amp = sum( ((t3amp_model - data.t3amp)./data.t3amp_err).^2);
  chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi)./data.t3phi_err).^2);
# centering
  rho = 1e4
  cent_g = zeros(size(x))
  cent_f = reg_centering(x, cent_g)
  # note: this is correct but slower
  g_v2 = 4.0*sum(((v2_model-data.v2)./data.v2_err.^2).*real(conj(cvis_model[data.indx_v2]).*dft[data.indx_v2,:]),1);
  g_t3amp = 2.0*sum(((t3amp_model-data.t3amp)./data.t3amp_err.^2).*
  (   real( conj(cvis_model[data.indx_t3_1]./abs.(cvis_model[data.indx_t3_1])).*dft[data.indx_t3_1,:]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3])       + real( conj(cvis_model[data.indx_t3_2]./abs.(cvis_model[data.indx_t3_2])).*dft[data.indx_t3_2,:]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3])+ real( conj(cvis_model[data.indx_t3_3]./abs.(cvis_model[data.indx_t3_3])).*dft[data.indx_t3_3,:]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2])),1);
  t3model_der = dft[data.indx_t3_1,:].*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3] + dft[data.indx_t3_2,:].*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3] + dft[data.indx_t3_3,:].*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2];
  g_t3phi =360./pi*sum(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*(-imag(t3_model).*real(t3model_der)+real(t3_model).*imag(t3model_der)),1);
  g[1:end] = vec(g_v2 + g_t3amp + g_t3phi) +  rho * cent_g;
  g[1:end] = (g - sum(vec(x).*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
  println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi," Flux: ", flux, " CENT: ", cent_f, " CDG ", cdg(reshape(x,nx,nx)))
  return chi2_v2 + chi2_t3amp + chi2_t3phi + rho * cent_f
end
