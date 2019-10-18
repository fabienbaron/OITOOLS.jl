#
# Examples of Model YSO images
#
using OITOOLS
using FITSIO

# Generate images - disc vs star

#image = readfits("./data/v1295_Aql.fits");
#imdisp(image.^.25)
#image[95:105, 95:105] .= 0;
#disk = image/sum(image)
#writefits(disk, "./data/disk.fits")
disk = readfits("./data/disk.fits")

lowK = 2.12e-6
hiK = 2.46e-6
λ0 = 0.5*(hiK+lowK) # central wavelength
Δλ = hiK-lowK # K band
#R = 30
#δλ = λ0 / R
#nwav = Δλ/δλ

nwav = 8 # spectral channels
λ = collect(range(lowK, stop=hiK, length=nwav))

disk_indx = 1.0
f0 = 0.8# fraction of total flux represented by the star at central wavelength
fstar = f0*(λ/λ0).^-4
fdisk = (1-f0)*(λ/λ0).^disk_indx
all_disks = (fdisk./(fstar+fdisk))*vec(disk)'
all_disks[:,201*101+101] = fstar./(fstar+fdisk)

# Generate images - YSO in several colors, From Claire Davies HD142666 paper
#using OITOOLS
include("../src/OITOOLS.jl"); using Main.OITOOLS
nwavs = 3#16
imageH = Float64.(readfits("./data/image_1p67.fits"))
#imageH ./=sum(imageH)
#imageH[75:76,75:76].=0
#imageH[75,75] = 0.2
imageK = Float64.(readfits("./data/image_2p13.fits"))
#imageK ./=sum(imageK)
#imageK[75:76,75:76].=0
#imageK[75,75] = 0.12

fade = collect(range(0,stop=1, length=nwavs))

dec=30.551191666666668 # AB Aur
facility_config_file="./data/CHARA.txt"
obsv_config_file="./data/default_obs.txt"
combiner_config_file="./data/MIRC.txt"

hour_angles=range(-1.70,stop=3.53,length=(108))
pixsize=0.1

low_wav = 1.51e-6
hi_wav  = 2.46e-6
wavelengths = collect(range(low_wav,stop=hi_wav, length=nwavs))
wavebands = [0.5*(wavelengths[2:end]-wavelengths[1:end-1]); 0.5*(wavelengths[end]-wavelengths[end-1])]

for i=1:nwavs
    image = (1-fade[i])*imageH + fade[i]*imageK
    image ./= sum(image)
    image_file="./data/image$i.fits"
    writefits(image, image_file)
    wave_config_file="./data/yso_wavelengths$i.txt"
    open(wave_config_file, "w") do f
    write(f, "combiner = \"MIRC\"\n")
    write(f, "mode = \"H_PRISM\"\n")
    write(f, "$(wavelengths[i]) $(wavebands[i])\n")
    end
    out_file="!./data/yso_sim$i.oifits"
    simulate_ha(facility_config_file,obsv_config_file,combiner_config_file,wave_config_file,hour_angles,dec,image_file,pixsize,out_file)
end

# Check ⁠- independent images
pixsize = 0.1; nx=150
for i=1:nwavs
data = readoifits("./data/yso_sim$i.oifits")[1,1];
x_true = readfits("./data/image$i.fits", normalize=true, vectorize=true)
ddft = setup_dft(data, nx, pixsize);
chi2_dft_f(x_true, ddft, data)

nx=150; ft = setup_nfft(data, nx, pixsize);
x_start = zeros(nx*nx); x_start[75+150*75]=1.0

regularizers = [   ["centering", 1e4] ];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = true, maxiter=100);
x = reconstruct(x, data, ft, regularizers = regularizers, verb = true, maxiter=100);
x = reconstruct(x, data, ft, regularizers = regularizers, verb = true, maxiter=300);
regularizers = [   ["centering", 1e4], ["tv", 2e4] ];
x = reconstruct(x, data, ft, regularizers = regularizers, verb = true, maxiter=300);
x = reconstruct(x, data, ft, regularizers = regularizers, verb = true, maxiter=300);
x = reconstruct(x, data, ft, regularizers = regularizers, verb = true, maxiter=300);
x = reconstruct(x, data, ft, regularizers = regularizers, verb = true, maxiter=300);
x = reconstruct(x, data, ft, regularizers = regularizers, verb = true, maxiter=300);
x = reconstruct(x, data, ft, regularizers = regularizers, verb = true, maxiter=300);
x = reconstruct(x, data, ft, regularizers = regularizers, verb = true, maxiter=300);

end

# CHARA realistic time observations
nwavs=3
oifitsfiles=["./data/yso_sim$i.oifits" for i=1:nwavs]
data = readoifits_multicolors(oifitsfiles)
printcolors = [:red, :green, :blue]
pixsize = 0.1;
nx = 150;
fftplan = setup_nfft_polychromatic(data, nx, pixsize);
x_start_single = zeros(nx*nx); x_start_single[75+150*75]=1.0
x_start = repeat(vec(x_start_single), nwavs)


x_true = [readfits("./data/image1.fits", normalize=true, vectorize=true);readfits("./data/image2.fits", normalize=true, vectorize=true);readfits("./data/image3.fits", normalize=true, vectorize=true)]

crit_polychromatic_nfft_fg(x_true, g, fftplan, data, verb=true)

#regularizers = [   [ ["centering", 1e4], ["tv", 2e4] ],  #wave 1
#                   [ ["tv", 2e4] ],                      #wave 2
#                   [ ["tv", 2e4] ],                      #wave 3
#                   [ ["spectral_tvsq", 3e8] ] ]; #transtemporal

regularizers = [  [ ["centering", 1e4] ],  #wave 1
                [    ["centering", 1e4]],                      #wave 2
                [    ["centering", 1e4]],                      #wave 3
                          [ ["spectral_tvsq", 1e5] ] ]; #transtemporal

x = reconstruct_polychromatic(x_start, data, fftplan, printcolor= printcolors, regularizers = regularizers, maxiter = 100);
x = reconstruct_polychromatic(x, data, fftplan, printcolor= printcolors, regularizers = regularizers, maxiter = 100);
x = reconstruct_polychromatic(x, data, fftplan, printcolor= printcolors, regularizers = regularizers, maxiter = 100);

imdisp_polychromatic(x.^.5, nwavs=3)

# Convolution of images with telescope of 300 m
