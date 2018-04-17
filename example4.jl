include("oitools.jl");
oifitsfile = "AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1];
fit_model_v2_nest(data, visibility_ud, [1.5:10.0]);
fit_model_v2_nest(data, visibility_ldpow, [1.0:10.0,0:1]);
fit_model_v2_nest(data, visibility_ldlin, [1.0:10.0,0:1]);
fit_model_v2_nest(data, visibility_ldquad, [1.0:10.0,0:1,0:1]);
