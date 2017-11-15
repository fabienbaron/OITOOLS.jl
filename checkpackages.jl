if Pkg.installed("NFFT")==nothing
    Pkg.add("NFFT")
end

if Pkg.installed("NearestNeighbors")==nothing
    Pkg.add("NearestNeighbors")
end

if Pkg.installed("OIFITS")==nothing
    Pkg.add("OIFITS")
    Pkg.checkout("OIFITS","master")
    println("Adding fixes to OIFITS.jl for OIFITS v2")
    cp("./fixes/oifile.jl",string(Pkg.dir("OIFITS"),"/src/"))
end
