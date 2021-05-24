using OITOOLS
  # Test 1: match a gaussian wrapped distribution with Von Mises
  μ=20/180*pi
  σ=10/180*pi
  k=gaussianwrapped_to_vonmises_fast.(σ)
  p0(x) = 1/(sqrt(2*pi)*σ)*(  exp.(-0.5*(x.-μ).^2/σ^2) )
  p1(x) = 1/(sqrt(2*pi)*σ)*(  exp.(-0.5*(x.-μ).^2/σ^2) + exp.(-0.5*(x.-μ.-2*pi).^2/σ^2) + exp.(-0.5*(x.-μ.+2*pi).^2/σ^2)+ exp.(-0.5*(x.-μ.-4*pi).^2/σ^2) + exp.(-0.5*(x.-μ.+4*pi).^2/σ^2))
  p2(x) = exp.(k*cos.(x.-μ))./(2*pi*besseli.(0,k))
  x = collect(range(-720,720,length=2000));
  clf()
  plot(x,-log.(p1(x/180*pi)))
  plot(x,-log.(p2(x/180*pi)))
  # note: the distributions start to differ when σ > .4 rad = 23 degrees

  # Test 2: range of estimation and maxiμm error
  σ = collect(range(0.01, 360., length=3600))/180*pi
  kk=gaussianwrapped_to_vonmises_fast.(σ)
  function mod2pijdm(x)
    mod2pi.(mod2pi.(x.+pi).+2*pi) .- pi
  end
  n= 50
  x  = rand(n)
  μ = rand(n)
  σ = 0.4 .+ 0.0.*abs.(rand(n))*10.
  k=gaussianwrapped_to_vonmises_fast._fast.(σ)
  # -2 log p(x) = ...
  println(" -2 log p_gausswrap(x): ", norm(mod2pijdm(x.-μ)./σ)^2+n*log(2*pi)+2*sum(log.(σ)),
  " -2 log p_von(x): ",  -2*sum(k.*cos.(x-μ)) + 2*n*log(2*pi) + 2*sum(logbesselI0.(k)))
  chi2r_gauss = norm(mod2pijdm(x-μ)./σ)^2/n
  chi2r_vonMises =  -2*sum(k.*cos.(x-μ))/n+ log(2*pi) + 2*sum(logbesselI0.(k)-log.(σ))/n
  #using OITOOLS
  # Test 3: BC2004 data
  oifitsfile = "./data/2004-data1.oifits";
  x0 = vec(readfits("./data/2004-64.fits"));
  x0 /= sum(x0)
  pixsize = 0.202;
  nx = 64;
  data = readoifits(oifitsfile)[1,1];
  data.t3phi_err = data.t3phi_err/180*pi;
  data.t3phi = data.t3phi/180*pi;
  H = setup_dft(data.uv, nx, pixsize);
  z = H*x0/sum(x0)
  t3phimodel = angle.(z[data.indx_t3_1].*z[data.indx_t3_2].*z[data.indx_t3_3]);
  lLike_gauss = norm(mod2pijdm.(data.t3phi - t3phimodel )./data.t3phi_err)^2 + data.nt3phi*log(2*pi) + 2*sum(log.(data.t3phi_err))
  chi2r_gauss =  norm(mod2pijdm.(data.t3phi - t3phimodel )./data.t3phi_err)^2/data.nt3phi
  kk = gaussianwrapped_to_vonmises_fast.(data.t3phi_err);
  lLike_vonMises = -2*sum(kk.*cos.(data.t3phi - t3phimodel))+ 2*data.nt3phi*log(2*pi) + 2*sum(logbesselI0.(kk))
  chi2r_vonMises =  -2*sum(kk.*cos.(data.t3phi - t3phimodel))/data.nt3phi+ log(2*pi) + 2*sum(logbesselI0.(kk)-log.(data.t3phi_err))/data.nt3phi#approximation of the equivalent chi2
  kk = gaussianwrapped_to_vonmises_fast.(data.t3phi_err);
  lLike_vonMises = sum(-2*kk.*cos.(data.t3phi - t3phimodel) .+ 2*log(2*pi) + 2*logbesselI0.(kk))
  chi2r_vonMises = sum(-2*kk.*cos.(data.t3phi - t3phimodel) .+ log(2*pi) + 2*(logbesselI0.(kk)-log.(data.t3phi_err)))/data.nt3phi#approximation of the equivalent chi2
