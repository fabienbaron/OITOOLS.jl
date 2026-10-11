# The model zoo of PLAN_ModelJacobian.md §4, shared by the benchmark and the comparisons.
using OITOOLS

const AXCIR = joinpath(pkgdir(OITOOLS), "demos", "data", "AX_Cir.oifits")

# A grey image for the image component: a 1.5 mas disc, 64 px at 0.05 mas.
const IMG = let p = joinpath(mktempdir(), "disc.fits"), nx = 64, ps = 0.05, c = nx ÷ 2 + 1
    writefits([hypot((i - c) * ps, (j - c) * ps) <= 0.75 ? 1.0 : 0.0 for i in 1:nx, j in 1:nx],
              p; pixsize = ps)
    p
end

# (name, model dict, free parameters)
const ZOO = [
    ("ud", Dict{String,Any}("s,ud" => 0.9), ["s,ud"]),
    ("gaussian", Dict{String,Any}("g,fwhm" => 0.9), ["g,fwhm"]),
    ("ldlin", Dict{String,Any}("s,ldlin" => 0.9, "s,u" => 0.3), ["s,ldlin", "s,u"]),
    ("ldquad", Dict{String,Any}("s,ldquad" => 0.9, "s,u" => 0.3, "s,w" => 0.1), ["s,ldquad", "s,u", "s,w"]),
    ("ldsqrt", Dict{String,Any}("s,ldsqrt" => 0.9, "s,u" => 0.3, "s,w" => 0.1), ["s,ldsqrt", "s,u", "s,w"]),
    ("ldpow", Dict{String,Any}("s,ldpow" => 0.9, "s,alpha" => 0.15), ["s,ldpow", "s,alpha"]),
    ("ldclaret4", Dict{String,Any}("s,ldclaret4" => 0.9, "s,c1" => 0.5, "s,c2" => -0.2, "s,c3" => 0.3, "s,c4" => -0.1),
        ["s,ldclaret4", "s,c1", "s,c2", "s,c3", "s,c4"]),
    ("binary", Dict{String,Any}("s,ud" => 0.83, "s,f" => "1 - \$c,f", "c,f" => 0.01, "c,x" => 6.2, "c,y" => -28.5),
        ["s,ud", "c,f", "c,x", "c,y"]),
    ("ring incl/pa", Dict{String,Any}("r,diamin" => 2.0, "r,diamout" => 3.0, "r,incl" => 40.0, "r,pa" => 30.0),
        ["r,diamin", "r,diamout", "r,incl", "r,pa"]),
    ("gaussian ring incl/pa", Dict{String,Any}("r,fwhmin" => 2.0, "r,fwhmout" => 3.0, "r,incl" => 40.0, "r,pa" => 30.0),
        ["r,fwhmin", "r,fwhmout", "r,incl", "r,pa"]),
    ("crescent", Dict{String,Any}("c,crin" => 2.0, "c,crout" => 3.0, "c,croff" => 0.5, "c,crprojang" => 20.0),
        ["c,crin", "c,crout", "c,croff", "c,crprojang"]),
    ("hankel profile", Dict{String,Any}("d,profile" => "exp(-(\$R/\$s)^\$p)", "d,s" => 1.0, "d,p" => 2.0, "d,udout" => 6.0),
        ["d,s", "d,p"]),
    ("hankel + az modes", Dict{String,Any}("d,profile" => "exp(-(\$R/\$s)^2)", "d,s" => 1.0, "d,udout" => 6.0,
        "d,az amp1" => 0.3, "d,az projang1" => 40.0), ["d,s", "d,az amp1", "d,az projang1"]),
    ("visfunc", Dict{String,Any}("v,visfunc" => "exp(-(\$BL*\$a)^2)", "v,a" => 0.01), ["v,a"]),
    ("image + point", Dict{String,Any}("i,image" => IMG, "i,f" => "1 - \$p,f", "p,f" => 0.2, "p,x" => 3.0, "p,y" => -2.0),
        ["p,f", "p,x", "p,y"]),
    ("chromatic (\$WL)", Dict{String,Any}("s,ud" => 0.9, "s,f" => "\$fs * (\$WL/1.6)^(-4)", "fs" => 0.6,
        "d,fwhm" => 3.0, "d,f" => "(1 - \$fs) * (\$WL/1.6)^\$slope", "slope" => 0.5), ["s,ud", "fs", "slope", "d,fwhm"]),
    ("derived chain", Dict{String,Any}("s,ud" => "2 * \$b", "b" => 0.45, "g,fwhm" => "\$b + 1", "g,f" => 0.3, "s,f" => 0.7),
        ["b"]),
    ("spatial kernel", Dict{String,Any}("s,ud" => 0.9, "spatial_kernel" => 0.3), ["s,ud", "spatial_kernel"]),
    ("binary smeared (3)", Dict{String,Any}("s,ud" => 0.83, "s,f" => "1 - \$c,f", "c,f" => 0.01, "c,x" => 6.2,
        "c,y" => -28.5, "smear" => 3), ["s,ud", "c,f", "c,x", "c,y"]),
    # Resolved stars: at 30 mas AX Cir's baselines span ζ = 8.5–23, past the kernels' short
    # series (ζ < 2) and across the switch to Hankel's expansion at ζ = 12.
    ("ldpow resolved", Dict{String,Any}("s,ldpow" => 30.0, "s,alpha" => 0.15), ["s,ldpow", "s,alpha"]),
    ("ldsqrt resolved", Dict{String,Any}("s,ldsqrt" => 30.0, "s,u" => 0.3, "s,w" => 0.1), ["s,ldsqrt", "s,u", "s,w"]),
    ("ldclaret4 resolved", Dict{String,Any}("s,ldclaret4" => 30.0, "s,c1" => 0.5, "s,c2" => -0.2, "s,c3" => 0.3,
        "s,c4" => -0.1), ["s,ldclaret4", "s,c1", "s,c2", "s,c3", "s,c4"]),
]

