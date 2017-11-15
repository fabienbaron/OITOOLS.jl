function check_packages()
if Pkg.installed("NFFT")==nothing
    Pkg.add("NFFT")
end

if Pkg.installed("NearestNeighbors")==nothing
    Pkg.add("NearestNeighbors")
end

if Pkg.installed("OIFITS")==nothing
    print_with_color(:red, "OIFITS package is not installed... Installing now\n");
    Pkg.add("OIFITS")
    Pkg.checkout("OIFITS","master")
    print_with_color(:red, "Adding fixes to OIFITS.jl for OIFITS v2\n");
    cp("./fixes/oifile.jl",string(Pkg.dir("OIFITS"),"/src/oifile.jl"),remove_destination=true)
end
end

check_packages()
