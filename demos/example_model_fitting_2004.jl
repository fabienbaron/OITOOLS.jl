# Model fittng example: multiple components
# Trying to model-fit the 2004 Beauty Contest data set
using OITOOLS

data= readoifits("./data/2004-data1.oifits")[1,1];

# Setup components for the 2004 Model
c1 = create_component(type="ud", name="Star")

c1.vis_params[1].free = true; # diameter in mas
c1.vis_params[1].maxval = 2.0; # diameter in mas
c1.spectrum_params[1].free = true;

c2 = create_component(type="ring", name="pAGB-ring")
c2.spectrum_params[1].free = true;
# c2.vis_params[5].free = false
# c2.vis_params[6].free = false
# c2.vis_params[7].free = false
# c2.vis_params[8].free = false
# c2.vis_params[9].free = false

# Put the model together
model = create_model(c1,c2)
display(model) # tis shows all the model parameters and does a few tests


minf, minx, cvis_model, result = fit_model_ultranest(data, model, min_num_live_points = 400, cluster_num_live_points = 100)
# Smaller numbers -> faster convergence but worse posterior characterization
# For optimization, min_num_live_points= 100-400 works fine
# For logZ estimation, >1000 would typically be OK
