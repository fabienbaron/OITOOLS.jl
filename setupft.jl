#using NFFT
# NFFT
#plan = NFFTPlan(full_uv*scale_rad, (nx,nx));
#for n in 1:nw
#  cvis_model[:,n] = nfft(plan, image + 0. * im);
#end

#DFT -- warning, centering is different from NFFT
#xtransform = zeros(Complex{Float64},nuv, nx)
#ytransform = zeros(Complex{Float64},nuv, nx)
#for uu=1:nuv
# xtransform and ytransform use less memory than the full dft
# but the dft can take bandwidth smearing into account better
#    xtransform[uu, :] = exp(2 * pi * im * scale_rad * v2_u[uu] * [i for i=1:nx]); # could be  (ii - (nx+1) / 2) for centering
#    ytransform[uu, :] = exp(2 * pi * im * scale_rad * v2_v[uu] * [i for i=1:nx]);
#end

function setup_ft(data, nx)
dft = zeros(Complex{Float64}, data.nuv, nx*nx);
xvals = [((i-1)%nx+1) for i=1:nx*nx];
yvals=  [(div(i-1,nx)+1) for i=1:nx*nx];
for uu=1:data.nuv
    dft[uu,:] = exp(-2 * pi * im * scale_rad * (data.uv[1,uu] * xvals + data.uv[2,uu] * yvals));
end
return dft
end
