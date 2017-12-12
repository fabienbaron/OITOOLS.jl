include("oitools.jl");
include("satlas.jl")

oifitsfile = "AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1];
satlas_model = read_satlas("satlas_models/ld_satlas_surface.2t3300g-25m25");

#UNIFORM DISC
model_cvis=visibility_ud([6.0], data[1].v2_baseline)
model_v2 = abs2.(model_cvis[1:data[1].nv2]);
v2plot_modelvsdata(data[1].v2_baseline,data[1].v2,data[1].v2_err, model_v2)

#SATLAS IMAGE (computes an image then takes the DFT)
model_cvis = complex_vis_satlas_img(data[1].v2_baseline, 6.0, satlas_model);
model_v2 = abs2.(model_cvis[1:data[1].nv2]);
v2plot_modelvsdata(data[1].v2_baseline,data[1].v2,data[1].v2_err, model_v2)

#SATLAS ANALYTIC (pure analytic function, relies on Struve from python which may be flaky)
model_cvis = complex_vis_satlas_any(data[1].v2_baseline, 8, satlas_model)
model_v2 = abs2.(model_cvis[1:data[1].nv2]);
v2plot_modelvsdata(data[1].v2_baseline,data[1].v2,data[1].v2_err, model_v2)
