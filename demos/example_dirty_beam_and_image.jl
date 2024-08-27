using OITOOLS,PyPlot,NFFT
set_oiplot_defaults()
oifitsfile = "./data/2004-data1.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.
nx=128
pixsize=.101
ft = setup_nfft(data, nx, pixsize);

# Dirty Beam
V= Complex.(visibility_ud([0.0], data.uv)) # point source
img = real.(adjoint(ft[1])*V); img = img.*(img .>0); imdisp(img, pixsize=pixsize, colormap="gist_earth");
text(6,5.8,"Dirty Beam",color="white",size="large")
savefig("dirty-beam.png")

#Direct Inversion
x = readfits("./data/2004true.fits")[66:193,66:193]
imdisp(recenter(x,mask=x.>maximum(x)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5.8,"Truth",color="white",size="large")
savefig("dirty-truth.png")

V = image_to_vis(x,ft[1])
img = real.(adjoint(ft[1])*V); img = img.*(img .>0);
imdisp(recenter(img,mask=img.>maximum(img)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5.8,"Dirty Map - No Noise",color="white",size="large")
savefig("dirty-map.png")

# V_noisy = abs.(V).*(1 .+0/100*randn(length(V))).*cis.(0*pi/180*(angle.(V)*180/pi.+30*randn(length(V))))
# img = real.(nfft_adjoint(ft[1], V_noisy)); img = img.*(img .>0);
# imdisp(recenter(img,mask=img.>maximum(img)/10),pixsize=pixsize,colormap="gist_earth");
# text(-6,5.8,"Direct Inversion - 1% Noise on Amplitude",color="white",size="x-large")
# savefig("dirty-inversion-1-amp.png")

# estimate V and 1/sigma_V^2 from V2 and V2_err using equation 3.98a in Data Analysis (Sivia/Skilling)
V = sqrt.(0.5*( data.v2 +  sqrt.(data.v2.^2+2*data.v2_err.^2)))
img = real.(adjoint(ft[3])*Complex.(V)); img = img.*(img .>0);
imdisp(recenter(img,mask=img.>maximum(img)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5.8,"Dirty Map - from |V|^2",color="white",size="large")
savefig("dirty-inversion-v2.png")
