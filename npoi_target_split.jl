# Example to load one of NPOI files with multiple targets
include("oitools.jl");

filename = "HD48329_oidata.fits";
targetname =  "FKV0254";
master = OIFITS.load(filename);
db = OIFITS.select(master, "OI_TARGET")[1];
target_filter = db[:target_id][find(db[:target].==targetname)];
for i in target_filter
    run(`oifits-filter --clobber $filename $targetname-$i.oifits --target-id=$i --accept-flagged=0 --snr-min=0.01`);
end
allfiles = string.(targetname,"-",target_filter,".oifits")
allfiles=allfiles[isfile.(allfiles)] # remove possible empty ones
run(`rm $targetname.oifits`);
run(`oifits-merge $targetname.oifits $allfiles`);
run(`rm $allfiles`);

# Then we need to filter bad data
data = (readoifits("$targetname.oifits"))[1,1];
good = find(  (data.v2_err.>0) .& (data.v2_err.<1.0) .& (data.v2.>-0.2) .& (data.v2.<1.2) )
data.v2 = data.v2[good]
data.v2_err = data.v2_err[good]
data.v2_baseline = data.v2_baseline[good]
data.nv2 = length(data.v2)

f_chi2, params, cvis_model = fit_model_v2(data, visibility_ud, [1.0]);# diameter is the parameter
