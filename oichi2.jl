function mod360(x)
mod(mod(x+180,360.)+360., 360.) - 180.
end

function cvis_to_v2(cvis, indx)
  v2_model = abs2(cvis[indx]);
end

function cvis_to_t3(cvis, indx1, indx2, indx3)
  t3 = cvis[indx1].*cvis[indx2].*cvis[indx3];
  t3amp = abs(t3);
  t3phi = angle(t3)*180./pi;
  return t3, t3amp, t3phi
end

#fig = figure("Image",figsize=(10,10));imshow(rotl90(image));PyPlot.draw();PyPlot.pause(1);

function image_to_cvis(x, dft)
  flux = sum(x);
  nuv = size(dft, 1)
  cvis_model = zeros(Complex{Float64},data.nuv);
  cvis_model = dft * x / flux;
end


function chi2(x, dft, data, verbose = true)
nx2 = length(x)
flux = sum(x);
cvis_model = zeros(Complex{Float64},div(data.nuv,data.nw),data.nw);
cvis_model[:,1] = dft * x / flux;
# compute observables from all cvis
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
chi2_v2 = sum( ((v2_model - data.v2_data)./data.v2_data_err).^2);
chi2_t3amp = sum( ((t3amp_model - data.t3amp_data)./data.t3amp_data_err).^2);
chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi_data)./data.t3phi_data_err).^2);
if verbose == true
  println("V2: ", chi2_v2, " T3A: ", chi2_t3amp, " T3P: ", chi2_t3phi," Flux: ", flux)
  println("chi2V2: ", chi2_v2/data.nv2, " chi2T3A: ", chi2_t3amp/data.nt3amp, " chi2T3P: ", chi2_t3phi/data.nt3phi," Flux: ", flux)
end
return chi2_v2 + chi2_t3amp + chi2_t3phi
end

function crit_fg(x, g, dft, data ) # criterion function plus its gradient w/r x
nx2 = length(x);
flux = sum(x);
cvis_model = zeros(Complex{Float64},data.nuv);
cvis_model = dft * x / flux;
# compute observables from all cvis
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
chi2_v2 = sum( ((v2_model - data.v2_data)./data.v2_data_err).^2);
chi2_t3amp = sum( ((t3amp_model - data.t3amp_data)./data.t3amp_data_err).^2);
chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi_data)./data.t3phi_data_err).^2);
g_v2 = Array(Float64, nx2);
g_t3amp = Array(Float64, nx2);
g_t3phi = Array(Float64, nx2);
#imdisp(x)
# note: this is correct but slower
#g = sum(4*((v2_model-v2_data)./v2_data_err.^2).*real(conj(cvis_model[indx_v2]).*dft[indx_v2,:]),1)
for ii = 1:nx2
g_v2[ii] = 4*sum(((v2_model-data.v2_data)./data.v2_data_err.^2).*real(conj(cvis_model[data.indx_v2]).*dft[data.indx_v2,ii]))
g_t3amp[ii] = 2*sum(((t3amp_model-data.t3amp_data)./data.t3amp_data_err.^2).*
                  (   real( conj(cvis_model[data.indx_t3_1]./abs(cvis_model[data.indx_t3_1])).*dft[data.indx_t3_1,ii]).*abs(cvis_model[data.indx_t3_2]).*abs(cvis_model[data.indx_t3_3])
                    + real( conj(cvis_model[data.indx_t3_2]./abs(cvis_model[data.indx_t3_2])).*dft[data.indx_t3_2,ii]).*abs(cvis_model[data.indx_t3_1]).*abs(cvis_model[data.indx_t3_3])
                    + real( conj(cvis_model[data.indx_t3_3]./abs(cvis_model[data.indx_t3_3])).*dft[data.indx_t3_3,ii]).*abs(cvis_model[data.indx_t3_1]).*abs(cvis_model[data.indx_t3_2]))
                   );
t3model_der = dft[data.indx_t3_1,ii].*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3] + dft[data.indx_t3_2,ii].*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3] + dft[data.indx_t3_3,ii].*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2];
g_t3phi[ii] = sum(360.0/pi*((mod360(t3phi_model-data.t3phi_data)./data.t3phi_data_err.^2)./abs2(t3_model)).*(-imag(t3_model).*real(t3model_der)+real(t3_model).*imag(t3model_der));
);
end
imdisp(x)
g[1:nx2] = g_v2 + g_t3amp + g_t3phi
g[1:nx2] = (g - sum(x.*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi," Flux: ", flux)
return chi2_v2 + chi2_t3amp + chi2_t3phi
end

function gaussian2d(n,m,sigma)
g2d = [exp(-((X-(m/2)).^2+(Y-n/2).^2)/(2*sigma.^2)) for X=1:m, Y=1:n]
return g2d
end

function cdg(x) #note: takes a 2D array
xvals=[i for i=1:size(x,1)]
return [sum(xvals'*x) sum(x_true*xvals)]/sum(x)
end

function crit_fgreg(x, g, dft, data, rho, x0) # criterion with regularization
nx2 = length(x)
flux = sum(x);
cvis_model = zeros(Complex{Float64},div(data.nuv,data.nw),data.nw);
cvis_model[:,1] = dft * x / flux;
# compute observables from all cvis
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
chi2_v2 = sum( ((v2_model - data.v2_data)./data.v2_data_err).^2);
chi2_t3amp = sum( ((t3amp_model - data.t3amp_data)./data.t3amp_data_err).^2);
chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi_data)./data.t3phi_data_err).^2);
g_v2 = Array(Float64, nx2);
g_t3amp = Array(Float64, nx2);
g_t3phi = Array(Float64, nx2);

# regularization
reg = mu*sum( (x-x0).^2);
reg_der = 2*mu*sum(x-x0);

#spatial_grad = [x_true[1:end-1]-x_true[2:end];0]
#reg_f = 1e6*sum(abs(spatial_grad));
#reg_g = 1e6*sum(convert(Array{Float64,1},spatial_grad.>0)-convert(Array{Float64,1},spatial_grad.<0));

# note: this is correct but slower
#g = sum(4*((v2_model-v2_data)./v2_data_err.^2).*real(conj(cvis_model[indx_v2]).*dft[indx_v2,:]),1)
for ii = 1:nx2
g_v2[ii] = 4*sum(((v2_model-data.v2_data)./data.v2_data_err.^2).*real(conj(cvis_model[data.indx_v2]).*dft[data.indx_v2,ii]))
g_t3amp[ii] = 2*sum(((t3amp_model-data.t3amp_data)./data.t3amp_data_err.^2).*
                  (   real( conj(cvis_model[data.indx_t3_1]./abs(cvis_model[data.indx_t3_1])).*dft[data.indx_t3_1,ii]).*abs(cvis_model[data.indx_t3_2]).*abs(cvis_model[data.indx_t3_3])
                    + real( conj(cvis_model[data.indx_t3_2]./abs(cvis_model[data.indx_t3_2])).*dft[data.indx_t3_2,ii]).*abs(cvis_model[data.indx_t3_1]).*abs(cvis_model[data.indx_t3_3])
                    + real( conj(cvis_model[data.indx_t3_3]./abs(cvis_model[data.indx_t3_3])).*dft[data.indx_t3_3,ii]).*abs(cvis_model[data.indx_t3_1]).*abs(cvis_model[data.indx_t3_2]))
                   );
t3model_der = dft[data.indx_t3_1,ii].*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3] + dft[data.indx_t3_2,ii].*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3] + dft[data.indx_t3_3,ii].*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2];
g_t3phi[ii] = sum(360./pi*((mod360(t3phi_model-data.t3phi_data)./data.t3phi_data_err.^2)./abs2(t3_model)).*(-imag(t3_model).*real(t3model_der)+real(t3_model).*imag(t3model_der));
);
end
g[1:end] = g_v2 + g_t3amp + g_t3phi +  rho * reg_der;
g[1:end] = (g - sum(x.*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi," Flux: ", flux)
return chi2_v2 + chi2_t3amp + chi2_t3phi + rho *reg
end



function proj_positivity(ztilde)
z = copy(ztilde)
z[ztilde.>0]=0
return z
end
