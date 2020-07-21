#
# BOOTSTRAP EXAMPLE
#
using OITOOLS
# VLTI example
oifitsfile = "./data/AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1];
# since we only want to use V2 to fit, we set weights=[1.0,0.0,0.0]
params_mode, params_mean, params_err = bootstrap_fit(1000, data, visibility_ud, [8.0], weights=[1.0,0.0,0.0]);
params_mode, params_mean, params_err = bootstrap_fit(1000, data, visibility_ldquad, [8.0,0.1,0.1], weights=[1.0,0.0,0.0]);
