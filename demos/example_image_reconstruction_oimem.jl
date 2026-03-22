using OITOOLS
oifitsfile = "./data/2004-data1.oifits"
pixsize  = 0.101   # mas/pixel; set to 0 to use auto_pixsize
#oifitsfile = "./data/betlyr6t.oifits"
#pixsize = 0.05
nx       = 128
flux_err = 1e-5   # tight flux constraint
data = readoifits(oifitsfile);
ft = setup_ft(data, nx, pixsize);

# Prior: Gaussian centred on the image, FWHM = 1/5 of FoV
prior_fwhm = nx * pixsize / 5.0
prior = gaussian_prior(nx, pixsize; fwhm_mas = prior_fwhm);

x = reconstruct_bsmem(prior, data[1], ft[1];
                       regularizers = [["mem", prior]],
                       method       = [4, 1, 1, 2],
                       maxiter      = 100,
                       verbose      = true,
                       flux_err     = flux_err,
                       nrand        = 10,
                       #mackay_alpha = true,  # MacKay fixed-point α update
                       #ritz_alpha   = true,  # Ritz-value bisection α update
                       )
imdisp(x; pixsize)
outfile = replace(oifitsfile, r"\.oifits?$"i => "_oimem.fits")
writefits(x, outfile; pixsize)
println("Written: $outfile")
