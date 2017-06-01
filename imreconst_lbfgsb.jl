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

#  x_true = (read((FITS(fitsfile))[1])); nx = (size(x_true))[1]; x_true=vec(x_true)/sum(x_true);



nx = 64;
  data = read_oifits(oifitsfile);
  dft = setup_ft(data, nx);

  #initial image is a simple Gaussian
  x_start = Array(Float64, nx, nx)
  for i=1:nx
    for j=1:nx
      x_start[i,j] = exp(-((i-(nx+1)/2)^2+(j-(nx+1)/2)^2)/(2*(nx/6)^2))
    end
  end
  x_start = vec(x_start)/sum(x_start);




crit = (x,g)->crit_fgreg(x, g , dft, data);
f, x, numCall, numIter, status = lbfgsb( crit, x_start, lb=zeros(size(x_start)), ub=ones(size(x_start)), iprint=1);
