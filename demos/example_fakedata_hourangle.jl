using OITOOLS


#
#EXAMPLE 8
#Simulate an observation using an input image, given telescope parameters, and input hour angles

days=[30]
months=[10]
hours=[0]
minutes=[0]
secs=[0]
years=[2018]
longitude=118.0570313111
lstsarr=[]
hour_angles=[]
ra=[20,57,59.44]
for i=1:length(days)
    for j=1:length(months)
        for k=1:length(hours)
            for l=1:length(minutes)
                for m=1:length(secs)
                    for n=1:length(years)
                        lsts,hour_angle=hour_angle_calc(days[i],months[j],years[n],hours[k],minutes[l],secs[m],longitude,ra)
                        push!(lstsarr,lsts)
                        push!(hour_angles,hour_angle)
                    end
                end
            end
        end
    end
end
lst_hours=floor.(lstsarr)
lstmin=floor.((lst.-floor.(lstsarr)).*60)
lstsec=(((lst.-floor.(lstsarr)).*60)-floor.(lstmin))*60
facility_config_file="./data/example_facility_config.txt"
obsv_config_file="./data/example_obs_config.txt"
combiner_config_file="./data/example_combiner_config.txt"
wave_config_file="./data/example_wave_config.txt"

image_file="./data/2004true.fits"
pixsize=0.101
out_file="!./data/2004testsimulation.oifits"

simulate_ha(facility_config_file,obsv_config_file,combiner_config_file,wave_config_file,hour_angles,RA,image_file,pixsize,out_file)

#Compare simulated data to impate

#display the image
fitsfile = "./data/2004true.fits";
pixsize = 0.101; # in mas/pixel
x_true = readfits(fitsfile); nx = (size(x_true))[1]; x_true=vec(x_true);
imdisp(x_true, pixscale = pixsize, tickinterval = 1.0, beamsize = 1.0, beamlocation = [0.85, 0.85])

#read the data file

oifitsfile = "./data/2004testsimulation.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.
uvplot(data);
v2plot(data,logplot=true);# Alternatively, one can do v2plot(data.v2_baseline,data.v2,data.v2_err,logplot=true);
t3phiplot(data);

# Setup Fourier transform via DFT
dft = setup_dft(data, nx, pixsize);
# This computes the complete chi2
f_chi2 = chi2_dft_f(x_true, dft, data);
# Compute |V|^2 observables and plot
cvis_model = image_to_cvis_dft(x_true, dft);
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model);

# NFFT method
ft = setup_nfft(data, nx, pixsize);
# This computes the complete chi2
f_chi2 = chi2_nfft_f(x_true, ft, data);
# Compute |V|^2 observables and plot
cvis_model = image_to_cvis_nfft(x_true, ft);
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model);
