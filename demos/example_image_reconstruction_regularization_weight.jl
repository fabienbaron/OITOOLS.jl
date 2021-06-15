#
# Very Basic Image reconstruction code
#
using OITOOLS,PyPlot
oifitsfile = "./data/2004-data1.oifits"
pixsize = 0.101 # size of a pixel in milliarcseconds
nx = 128 # width of image (number of pixels)
data = readoifits(oifitsfile)[1,1];
ft = setup_nfft(data, nx, pixsize);
#initial image is a simple Gaussian
x_start = gaussian2d(nx,nx,nx/6);
x_start = vec(x_start)/sum(x_start);

x = readfits("./data/2004true.fits")[64:195,64:195]
imdisp(recenter(x,mask=x.>maximum(x)/10),pixscale=pixsize,colormap="gist_earth");
text(-6,6,"Truth",color="white",size="xx-large")
savefig("tvsqtruth.png")


img=zeros(128,128,8)

regularizers = [["centering", 1e3], ["tvsq",1e4]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixscale=pixsize,colormap="gist_earth");
text(-6,6,"μ=1e4",color="white",size="xx-large")
savefig("tvsq1e4.png")

regularizers = [["centering", 1e3], ["tvsq",1e5]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixscale=pixsize,colormap="gist_earth");
text(-6,6,"μ=1e5",color="white",size="xx-large")
savefig("tvsq1e5.png")

regularizers = [["centering", 1e3], ["tvsq",1e6]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixscale=pixsize,colormap="gist_earth");
text(-6,6,"μ=1e6",color="white",size="xx-large")
savefig("tvsq1e6.png")

regularizers = [["centering", 1e3], ["tvsq",1e7]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixscale=pixsize,colormap="gist_earth");
text(-6,6,"μ=1e7",color="white",size="xx-large")
savefig("tvsq1e7.png")

regularizers = [["centering", 1e3], ["tvsq",1e8]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixscale=pixsize,colormap="gist_earth");
text(-6,6,"μ=1e8",color="white",size="xx-large")
savefig("tvsq1e8.png")

regularizers = [["centering", 1e3], ["tvsq",1e9]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixscale=pixsize,colormap="gist_earth");
text(-6,6,"μ=1e9",color="white",size="xx-large")
savefig("tvsq1e9.png")

regularizers = [["centering", 1e3], ["tvsq",1e10]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = false, maxiter=500);
imdisp(recenter(x,mask=x.>maximum(x)/10),pixscale=pixsize,colormap="gist_earth");
text(-6,6,"μ=1e10",color="white",size="xx-large")
savefig("tvsq1e10.png")
