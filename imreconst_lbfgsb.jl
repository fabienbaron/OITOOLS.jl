  using FITSIO
  using Lbfgsb
  include("readoifits.jl")
  include("setupft.jl")
  include("oichi2.jl")
  include("oiplot.jl")
  PyPlot.show()
  ##########################################
  ##########################################
  #
  # Code actually starts
  #
  pixellation = 0.2; # in mas
  fitsfile = "2004true64.fits";
  oifitsfile = "2004-data1.oifits";
  nw = 1;# monochromatic mode
  #read model fits file
  scale_rad = pixellation * (pi / 180.0) / 3600000.0;
  x_true = (read((FITS(fitsfile))[1])); nx = (size(x_true))[1]; x_true=vec(x_true)/sum(x_true);

  data = read_oifits(oifitsfile);
  dft = setup_ft(data, nx);
  #init required because of OptimPack way

  # initial values: Z, rho, U, mu
  x_start =  rand(size(x_true));
  x_start = vec(x_start)/sum(x_start);

  #regularization param
  crit = (x,g)->crit_fgreg(x, g , dft,data);
#  f, x, numCall, numIter, status = lbfgsb( crit, x_start, lb=zeros(size(x_start)), ub=ones(size(x_start)), iprint=1);
   f, x, numCall, numIter, status = lbfgsb( crit, x_start, lb=zeros(size(x_start)), ub=ones(size(x_start)), iprint=1);




#x_start = x_true
#x_start /= sum(x_start)
#g = zeros(size(x_true)); f_true = crit_fg(x_true, g, dft, data);

#x = OptimPack.nlcg(cost!, start_x, OptimPack.NLCG_HAGER_ZHANG)
#println(crit_fg(x, g, dft, data));
