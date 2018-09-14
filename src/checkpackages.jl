if VERSION==v"0.7.0"
using Pkg
end

function check_packages()
#SpecialFunction
if in("SpecialFunctions",keys(Pkg.installed()))
    @eval using SpecialFunctions
else
    @warn("SpecialFunctions not installed")
    Pkg.add("SpecialFunctions")
end
    
#NFFT 
if in("NFFT",keys(Pkg.installed()))
    @eval using NFFT
else
    @warn("NFFT not installed")
    Pkg.add("NFFT"); 
    Pkg.checkout("NFFT","master")
end

if in("NearestNeighbors",keys(Pkg.installed()))
    @eval using NearestNeighbors
else
    @warn("NearestNeighbors not installed")
    Pkg.add("NearestNeighbors")
end

if in("NLopt",keys(Pkg.installed()))
    @eval using NLopt
else
    @warn("NLopt not installed")
    Pkg.add("NLopt")
end

if in("OptimPackNextGen",keys(Pkg.installed()))
    @eval using OptimPackNextGen
else
    @warn("OptimPackNextGen not installed")
end

if in("OIFITS",keys(Pkg.installed()))
    @eval using OIFITS
else
    @warn("OIFITS not installed")
    print_with_color(:red, "OIFITS package is not installed... Installing now\n");
    Pkg.add("OIFITS")
    #Pkg.checkout("OIFITS","master")
    print_with_color(:red, "Adding fixes to OIFITS.jl for OIFITS v2\n");
    cp("./fixes/oifile.jl",string(Pkg.dir("OIFITS"),"/src/oifile.jl"),remove_destination=true)
end
end

check_packages()
