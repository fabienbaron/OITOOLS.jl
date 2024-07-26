#
# Demonstrating different types of classic regularizations
#
using OITOOLS, PyPlot
set_oiplot_defaults()
oifitsfile = "./data/2004-data1.oifits"
pixsize = 0.1 # size of a pixel in milliarcseconds
nx = 137 # width of image (number of pixels)
data = readoifits(oifitsfile)[1,1];
ft = setup_nfft(data, nx, pixsize);
#initial image is a simple Gaussian
x_start = gaussian2d(nx,nx,nx/6);

x=copy(x_start)
imdisp(recenter(x,mask=x.>maximum(x)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5,"Starting Image",color="white",size="xx-large")
savefig("types-starting-image.png")

regularizers = [["centering", 1e3]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5,"Positivity only",color="white",size="xx-large")
savefig("types-positivity-only.png")


regularizers = [["centering", 1e3], ["l1hyp",1e3]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5,"l1 hyperbolic approx",color="white",size="xx-large")
savefig("types-l1hyp.png")


regularizers = [["centering", 1e3], ["tvsq",5e7]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5,"Squared Total Variation",color="white",size="xx-large")
savefig("types-tvsq.png")

regularizers = [["centering", 1e3], ["tv",1e2]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5,"Total Variation",color="white",size="x-large")
savefig("types-tv.png")

regularizers = [["centering", 1e3], ["l2sq",7e5]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5,"Tikhonov",color="white",size="xx-large")
savefig("types-l2sq.png")

regularizers = [["centering", 1e3], ["l1l2", 7e6, 1e-3]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5,"l1l2",color="white",size="xx-large")
savefig("types-l1l2.png")

regularizers = [["centering", 1e3], ["entropy", 7e2]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5,"Maximum Entropy",color="white",size="xx-large")
savefig("types-entropy.png")

regularizers = [["centering", 1e3], ["compactness", 7e5, 20.0]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5,"Compactness",color="white",size="xx-large")
savefig("types-compactness.png")

# Example of using weights
regularizers = [["centering", 1e3], ["l1l2w", 1e6]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5,"l1l2w - All data",color="white",size="xx-large")
savefig("types-l1l2w.png")

regularizers = [["centering", 1e3], ["l1l2w", 1e7]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, weights=[1.0,0.0,0.0], verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5,"l1l2w - V2 only",color="white",size="xx-large")
savefig("types-l1l2w-v2only.png")

regularizers = [["centering", 1e3], ["l1l2w", 1e6]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, weights=[0.0,1.0,1.0], verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5,"l1l2w - T3PHI + T3AMP",color="white",size="xx-large")
savefig("types-l1l2w-fullbispectrum.png")

regularizers = [["centering", 1e3], ["l1l2w", 2e7]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, weights=[0.0,0.0,1.0], verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixsize=pixsize,colormap="gist_earth");
text(6,5,"l1l2w - T3PHI only",color="white",size="xx-large")
savefig("types-l1l2w-closure-only.png")
