#
# BOOTSTRAP EXAMPLE
#
using OITOOLS
# VLTI example
oifitsfile = "./data/AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1];

model = create_model(create_component(type="ldlin", name="Model"));

# since we only want to use V2 to fit, we set weights=[1.0,0.0,0.0]
params_mode, params_mean, params_err = bootstrap_fit(1000, data, model, fitter=:LN_NELDERMEAD, weights=[1.0,0.0,0.0]);
params_mode, params_mean, params_err = bootstrap_fit(1000, data, model, fitter=:Levenberg, weights=[1.0,0.0,0.0]);
