#
# Basic use of "SPARCO"-type reconstruction
#
using OITOOLS
oifitsfile = "./data/MWC275_T4a.oifits"
pixsize = 0.075
data = readoifits(oifitsfile, filter_bad_data=true)[1,1];
nx = 128
ft = setup_nfft(data, nx, pixsize);
# Parameters
# 1: flux fraction of star at λ0
# 2: flux fraction of background at λ0
# 3: stellar angular diameter
# 4: spectral index of the environment+4
# 5: λ0
params_start=[0.46, 0., 0., 0., 1.6e-6] #V2:12.36 T3A: 4.14 T3P: 1.61
x_start = gaussian2d(nx,nx,nx/6);
regularizers = [["l1l2", 2e8, 1e-6]];
params, x = reconstruct_sparco_gray(x_start, params_start, data, ft, regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
imdisp(x, pixsize = pixsize)
minchi2, params,ret = optimize_sparco_parameters(params, x, ft, data; weights = [1.0,1.0,1.0] )
params, x = reconstruct_sparco_gray(x, params, data, ft, regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
imdisp(x, pixsize = pixsize)
minchi2, params,ret = optimize_sparco_parameters(params, x, ft, data; weights = [1.0,1.0,1.0] )
params, x = reconstruct_sparco_gray(x, params, data, ft, regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
imdisp(x, pixsize = pixsize)


#
# Write FITS image, and minimization information in the header
# 
using FITSIO
filename = "sparco_out"
f = FITS("$filename.fits", "w")
hdr_key = ["PIXSIZE","MINCHI2","FITMSG","STRFRAC","BGFRAC","STRDIAM","ENVSLOPE","REFWL",]
hdr_val = Any[pixsize, minchi2, string(ret), params...]
hdr_cmt = ["mas","","see OITOOLS","flux fraction of star at lamdba_0", "flux fraction of background at lambda_0",
"stellar angular diameter","spectral index of the environment plus 4 (gt 0)","lambda_0 reference wavelength",]
for (i, reg) in enumerate(regularizers)
    push!(hdr_key, "REG$i")
    push!(hdr_val, reg[2])
    push!(hdr_cmt, reg[1])
end
write(f, x; header=FITSHeader(hdr_key, hdr_val, hdr_cmt))
close(f)