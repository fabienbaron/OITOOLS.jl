# Model fittng example: multiple components
# Trying to model-fit the 2004 Beauty Contest data set
using OITOOLS

data= readoifits("./data/2004-data1.oifits")[1,1];

# Setup components for the 2004 Model
c1 = create_component(type="ud", name="star")
c1.spectrum_params[1].free = true; # we have two components, we just need to fit the flux of the first
c1.vis_params[1].free = false; # diameter in mas
c1.vis_params[1].val = 0.0; # maximum diameter in mas (we expect a barely reslved star)
c2 = create_component(type="ring", name="ring")
c2.spectrum_params[1].free = false; # we have two components, we just need to fit the flux of the first
c2.vis_params[5].free = false
c2.vis_params[6].free = false
c2.vis_params[7].free = false
c2.vis_params[8].free = false
c2.vis_params[9].free = false

# If one wants to fit the ring a little better, uncomment this
# More free parameters, so much much slower (with Ultranest)
# c2.vis_params[5].free = true
# c2.vis_params[6].free = true
# c2.vis_params[7].free = true
# c2.vis_params[8].free = true

# Put the model together
model = create_model(c1,c2)

# Note: after ceating the model, changing component parameters from free or fixed (or opposite) *will* affect the model
# We need to update the internal model parameter map using update_model(model) in such a case

display(model) # tis shows all the model parameters and does a few tests
minf, minx, cvis_model, result = fit_model_ultranest(data, model, min_num_live_points = 400, cluster_num_live_points = 200)




# Smaller numbers -> faster convergence but worse posterior characterization
# For optimization, min_num_live_points= 100-400 works fine
# For logZ estimation, >1000 would typically be OK
