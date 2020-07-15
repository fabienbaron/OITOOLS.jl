#
# BOOTSTRAP EXAMPLE
#

using OITOOLS
# VLTI example
oifitsfile = "./data/AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1];
params_mode, params_mean, params_err = bootstrap_v2_fit(1000, data, visibility_ud, [8.0]);
params_mode, params_mean, params_err = bootstrap_v2_fit(1000, data, visibility_ldquad, [8.0,0.1,0.1]);
