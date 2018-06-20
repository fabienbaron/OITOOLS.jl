function check_packages()
if Pkg.installed("NFFT")==nothing
    Pkg.add("NFFT")
    Pkg.checkout("NFFT","master")
end

if Pkg.installed("SpecialFunctions")==nothing
    Pkg.add("SpecialFunctions")
end

if Pkg.installed("NearestNeighbors")==nothing
    Pkg.add("NearestNeighbors")
end

if Pkg.installed("OIFITS")==nothing
    print_with_color(:red, "OIFITS package is not installed... Installing now\n");
    Pkg.add("OIFITS")
    #Pkg.checkout("OIFITS","master")
    print_with_color(:red, "Adding fixes to OIFITS.jl for OIFITS v2\n");
    cp("./fixes/oifile.jl",string(Pkg.dir("OIFITS"),"/src/oifile.jl"),remove_destination=true)
end

if Pkg.installed("NLopt")==nothing
Pkg.add("NLopt")
end

if Pkg.installed("OptimPack")==nothing
Pkg.add("OptimPack")
Pkg.checkout("OptimPack","nextgen")
end

end

check_packages()
