#
# BOOTSTRAP EXAMPLE
#

include("oitools.jl");
oifitsfile = "AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1];
bootstrap_v2_fit(1000, data, visibility_ud, [8.0]);
