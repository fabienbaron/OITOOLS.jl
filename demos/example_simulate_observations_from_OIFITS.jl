using OITOOLS
#
#EXAMPLE 9
# Simulate an observation using an input oifits file as base
# The uv coverage and SNR of data points will be copied

# Example 1 -- using a FITS image
image_file="./data/2004true.fits"
pixsize=0.101
in_oifits  = "./data/2004-data1.oifits"
out_oifits  = "./data/2004-simulated.oifits"
simulate_from_oifits(in_oifits,out_oifits,image=image_file,pixsize=pixsize)
# One can then compare simulated data to input
data1 = readoifits(in_oifits); plot_v2(data1);
data2 = readoifits(out_oifits); plot_v2(data2);

# Example 2 - One can also use a model instead of an image
param = Dict{String,Any}(
    "star,ldlin" => 3.0,
    "star,u"     => 0.15,
    "star,f"     => 1.0,
)
model = parse_model(param, String[])
simulate_from_oifits(in_oifits, out_oifits, flat_model=model, flat_params=Float64[])
