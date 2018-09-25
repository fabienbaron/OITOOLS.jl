using OITOOLS
using LinearAlgebra
include("satlas.jl");

oifitsfile = "./data/AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1];
satlas_model = read_satlas("./data/ld_satlas_surface.2t3300g-25m25");


# Compare visibiliy_satlas image and analytic methods
for i = 1:5
    cvis_model = visibility_satlas_img([6e0], data.v2_baseline, satlas_model);
end

for i = 1:5
    cvis_model = visibility_satlas_any([6e0], data.v2_baseline, satlas_model);
end

# FITS
#

#UNIFORM DISC
f_chi2, params, cvis_model = fit_model_v2(data, visibility_ud, [8.0]);# diameter is the parameter
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
chi2_v2 = norm((v2_model- data.v2)./data.v2_err)^2/data.nv2;
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model,logplot=true);

#SATLAS ANALYTIC (pure analytic function, relies on Struve from python which may be flaky)
#model_cvis = visibility__satlas_any(data[1].v2_baseline, 8, satlas_model)
f_chi2, params, cvis_model = fit_satlas_v2(data, visibility_satlas_any, [11.0], satlas_model);# diameter is the parameter
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
chi2_v2 = norm((v2_model- data.v2)./data.v2_err)^2/data.nv2;
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model,logplot=true);

#SATLAS IMAGE (computes an image then takes the DFT), slow
#model_cvis = visibility_satlas_img(data[1].v2_baseline, 6.0, satlas_model);
f_chi2, params, cvis_model = fit_satlas_v2(data, visibility_satlas_img, [8.0], satlas_model);# diameter is the parameter
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
chi2_v2 = norm((v2_model- data.v2)./data.v2_err)^2/data.nv2;

v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model,logplot=true);
