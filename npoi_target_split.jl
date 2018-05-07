# Example to load one of NPOI files with multiple targets
using OIFITS
filename = "HD48329_oidata.fits";
targetname =  "FKV0254";
master = OIFITS.load(filename);
db = OIFITS.select(master, "OI_TARGET")[1];
target_filter = find(db[:target].==targetname);
for i in target_filter
    run(`oifits-filter --clobber $filename $targetname-$i.oifits --target-id=$i`);
end
allfiles = string.(targetname,"-",target_filter,".oifits")
run(`oifits-merge $targetname.oifits $allfiles`);
