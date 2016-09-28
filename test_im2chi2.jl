using FITSIO
using OptimPack
using PyPlot,PyCall
PyPlot.show()
#@pyimport mpl_toolkits.axes_grid1 as axgrid

include("readoifits.jl")
include("setupft.jl")
include("chi2.jl")
include("oiplot.jl")
##########################################
##########################################
#
# Code actually starts
#
pixellation = 0.4; # in mas
fitsfile = "2004true64.fits";
oifitsfile = "2004-data1.oifits";
nw = 1;# monochromatic mode
#read model fits file
scale_rad = pixellation * (pi / 180.0) / 3600000.0;
x_true = (read((FITS(fitsfile))[1])); nx = (size(x_true))[1]; x_true=vec(x_true);

data = read_oifits(oifitsfile);
dft = setup_ft(data, nx);
#init required because of OptimPack way
#f_true = crit_fg(x_true, g, dft, data);


# Setup optimpack
mxvl = 20
mxtr = 10
stpmn = 1e-20
stpmx = 1e+20


# initial values: Z, rho, U, mu

#regularization param
mu = 1e5;
rho = 100.
u_old = ones(size(x_true))
z_old = rand(size(x_true));
z_old /= sum(z_old);

for t=1:10

  # X step
  xtilde = z_old - u_old / rho;
  x = proj_positivity(xtilde);
imdisp(x);
  # Z step
  ztilde = x + u_old / rho;
  crit_admm = (z,g)->crit_fg(z, g , dft,data, rho, ztilde);
  z = OptimPack.vmlm(crit_admm, z_old, 3, verb = true, grtol = 1e-6, gatol = 0, maxeval = mxvl, maxiter = mxtr, stpmin = stpmn, stpmax = stpmx, scaling = OptimPack.SCALING_OREN_SPEDICATO, lnsrch = OptimPack.MoreThuenteLineSearch(ftol = 1e-8, gtol = 0.95));

  # U step
  u = u_old + rho * (x-z);

  # end of iterations
  u_old = u;
  z_old = z;
  x_old = x;
end


#x_start = x_true
#x_start /= sum(x_start)
#g = zeros(size(x_true)); f_true = crit_fg(x_true, g, dft, data);

#x = OptimPack.nlcg(cost!, start_x, OptimPack.NLCG_HAGER_ZHANG)
#println(crit_fg(x, g, dft, data));
