include("oitools.jl");
oifitsfile = "AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1];
fit_model_v2_nest(data, visibility_ud, [10.0]);
fit_model_v2_nest(data, visibility_ldpow, [10.0,0.5]);
fit_model_v2_nest(data, visibility_ldlin, [10.0,0.5]);
fit_model_v2_nest(data, visibility_ldquad, [10.0,0.5,0.1]);
