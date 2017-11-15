function check_packages()
if Pkg.installed("NFFT")==nothing
    Pkg.add("NFFT")
end

if Pkg.installed("NearestNeighbors")==nothing
    Pkg.add("NearestNeighbors")
end

if Pkg.installed("OIFITS")==nothing
    println("OIFITS package is not installed... Installing now");
    Pkg.add("OIFITS")
    Pkg.checkout("OIFITS","master")
    println("Adding fixes to OIFITS.jl for OIFITS v2")
    cp("./fixes/oifile.jl",string(Pkg.dir("OIFITS"),"/src/oifile.jl"),remove_destination=true)
end
end

check_packages()
