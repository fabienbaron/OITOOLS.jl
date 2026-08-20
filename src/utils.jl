#
# This file includes general purpose functions
#

using LinearAlgebra,SparseArrays, PythonCall

function x_start_from_V2_dft(data, dft; λ = 1e7 , μ = 1e6 )
# estimate V and 1/sigma_V^2 from V2 and V2_err using equation 3.98a in Data Analysis (Sivia/Skilling)
nx2 = size(dft,2)
nx = Int(round(sqrt(nx2)))
V2 = 0.5*( data.v2 +  sqrt.(data.v2.^2+2*data.v2_err.^2))
V = sqrt.(V2)
W = spdiagm(0=>(1.0./V2+2*(3*V2-data.v2)./(data.v2_err.^2))) # 1/sigma^2
H = dft[data.indx_v2, :];
o = ones(nx); D_1D = spdiagm(-1=>-o[1:nx-1],0=>o); D = [kron(spdiagm(0=>ones(nx)), D_1D) ; kron(D_1D, spdiagm(0=>ones(nx)))]; DtD = D'*D;
y = real(H'*(W*H)+λ*DtD+μ*sparse(1.0I, nx2,nx2))\(real(H'*(W*V))); y=y.*(y.>=0);imdisp(reshape(y,nx,nx));
_chi2_f(y, dft, data)
end



function cdg(x) #2D array
    xvals=[i for i=1:size(x,1)]
    return [sum(xvals'*x) sum(x*xvals)]/sum(x)
end

"""
    recenter(x; mask=[], max=false)

Recenter an image by circular shifting so the centroid (or peak if `max=true`)
is at the image center. Works on 1D (vectorized square image) or 2D arrays.
If `mask` is provided, the centroid is computed from the mask instead.
"""
function recenter(x::AbstractArray{<:Real}; mask=[], max=false)
    if ndims(x)==1 # assume square array
        n=Int(sqrt(length(x)))
        xtemp=reshape(x,(n,n))
        if mask ==[]
            δ = round.(Int,(n+1)/2 .- cdg(xtemp))
        else
            δ = round.(Int,(n+1)/2 .- cdg(reshape(mask,(n,n))))
        end
        return vec(circshift(xtemp, (δ[1], δ[2])))
    else
        if mask ==[]
            δ = round.(Int,[size(x)[1] size(x)[2]]/2 .- cdg(x))
        else
            δ = round.(Int,[size(x)[1] size(x)[2]]/2 .- cdg(mask))
        end
        return circshift(x, (δ[1], δ[2]))
    end
end

function query_target_from_simbad(targetname)
    return pyimport("astroquery.simbad").Simbad.query_object(targetname)
end

"""
    sexagesimal_to_degrees(s) -> Float64

Parse a sexagesimal coordinate string into a single decimal value. Fields may be separated
by whitespace or colons, and one, two or three of them may be given.

The sign applies to the **whole** quantity, not just the leading field: `"-00 30 00"` is
`-0.5`, and `"-46 28 00.57"` is `-46.4668`. Reducing the parsed fields with a dot product
against `[1, 1/60, 1/3600]` — which is what callers used to do — gets both of those wrong,
because it negates only the degrees term.

```jldoctest
julia> sexagesimal_to_degrees("-46 28 00.5731825")
-46.46682588402778
```
"""
function sexagesimal_to_degrees(s::AbstractString)
    t = strip(s)
    isempty(t) && throw(ArgumentError("empty coordinate string"))
    neg = startswith(t, '-')
    if startswith(t, '-') || startswith(t, '+')
        t = t[nextind(t, 1):end]
    end
    parts = split(t, r"[\s:]+"; keepempty = false)
    isempty(parts) && throw(ArgumentError("cannot parse coordinate string: \"$s\""))
    value = 0.0
    scale = 1.0
    for p in parts
        value += parse(Float64, p) * scale
        scale /= 60
    end
    return neg ? -value : value
end

"""
    ra_dec_from_simbad(targetname) -> (ra_deg, dec_deg)

Resolve a target name through SIMBAD and return its J2000 right ascension and declination
in **decimal degrees**.

Both astroquery layouts are handled: releases up to ~0.4.7 return sexagesimal strings in
`RA`/`DEC` columns (RA in hours), newer ones return decimal degrees in `ra`/`dec`.

!!! note "Changed in 0.11"
    This used to return two `Vector{Float64}` of sexagesimal components, leaving the
    conversion to the caller. Every caller did it with a dot product that mishandled the
    sign of southern declinations. It now returns degrees directly.
"""
function ra_dec_from_simbad(targetname)
    res = query_target_from_simbad(targetname)
    pyis(res, pybuiltins.None) && error("SIMBAD returned no match for \"$targetname\"")
    cols = Set(pyconvert(Vector{String}, res.colnames))
    if ("ra" in cols) && ("dec" in cols)            # astroquery >= 0.4.8: decimal degrees
        return pyconvert(Float64, res["ra"][0]), pyconvert(Float64, res["dec"][0])
    elseif ("RA" in cols) && ("DEC" in cols)        # older: sexagesimal strings, RA in hours
        ra_h = sexagesimal_to_degrees(pyconvert(String, res["RA"][0]))
        dec  = sexagesimal_to_degrees(pyconvert(String, res["DEC"][0]))
        return 15 * ra_h, dec
    else
        error("SIMBAD result for \"$targetname\" has no recognised coordinate columns; " *
              "got $(sort(collect(cols)))")
    end
end

"""
    magnitudes_from_simbad(targetname)

Query SIMBAD for photometric magnitudes (V, J, H, K, L, M, N).
Returns a Dict{String,Float64} with band names as keys.
Missing magnitudes are set to NaN.
"""
function magnitudes_from_simbad(targetname)
    Simbad = pyimport("astroquery.simbad").Simbad
    # B through N: V is what the adaptive optics guides on, R appears in older catalogues,
    # and the near-infrared bands are what the observation is actually made in.
    bands = ["B", "V", "R", "I", "J", "H", "K", "L", "M", "N"]
    mags = Dict{String,Float64}(b => NaN for b in bands)
    try
        # Use TAP query to get all available fluxes for this object
        query = """
        SELECT f.filter, f.flux
        FROM basic AS b
        JOIN ident AS i ON i.oidref = b.oid
        JOIN flux AS f ON f.oidref = b.oid
        WHERE i.id = '$(replace(targetname, "'" => "''"))'
        """
        result = Simbad.query_tap(query)
        # query_tap returns Python None when there is no match; ask Python rather than
        # comparing against Julia's nothing.
        if !pyis(result, pybuiltins.None)
            filters = result["filter"]
            fluxes  = result["flux"]
            nrows   = pylen(filters)
            for row in 0:(nrows-1)          # Python row index, passed through by Py getindex
                filt = uppercase(strip(pyconvert(String, filters[row])))
                for b in bands
                    if filt == b
                        try
                            mags[b] = pyconvert(Float64, fluxes[row])
                        catch; end
                    end
                end
            end
        end
    catch err
        # NOT swallowed. A network failure and a target with no photometry are different
        # answers, and reporting them the same way is how "this star has no K magnitude"
        # comes to mean "the query never ran".
        throw(ErrorException("SIMBAD photometry query failed for \"$targetname\": " *
                             first(split(sprint(showerror, err), "\n"))))
    end
    return mags
end

"""
Photometric bands `simbad_target` asks for, in the order a panel should show them.

B through N. The near-infrared ones are what an interferometric observation is actually made
in, and V is what the adaptive optics guides on, so both matter and for different reasons.
"""
const SIMBAD_BANDS = ["B", "V", "R", "I", "J", "H", "K", "L", "M", "N"]

"""
    simbad_target(name) -> NamedTuple

Everything SIMBAD holds that an observation needs: `(; name, main_id, ra, dec, pmra, pmdec,
plx, rv, sptype, otype, mags)`.

`ra`/`dec` are DEGREES, proper motions mas/yr, parallax mas, radial velocity km/s. Missing
values are `NaN` and a missing string is `""` — SIMBAD not knowing a target's parallax is a
fact about the target, not a failure, and it has to be distinguishable from the query failing.

Beyond coordinates and magnitudes this returns what ASPRO also reads, and for the same reasons:
proper motion places the target at the epoch actually being observed, parallax turns an angular
diameter into a physical one, and the spectral type is what a surface-brightness relation needs
to predict a diameter before anything is measured.
"""
function simbad_target(targetname::AbstractString)
    Simbad = pyimport("astroquery.simbad").Simbad
    name = String(targetname)
    esc  = replace(name, "'" => "''")

    num(x) = try
        v = pyconvert(Float64, x); isfinite(v) ? v : NaN
    catch
        NaN
    end
    str(x) = try
        strip(pyconvert(String, x))
    catch
        ""
    end

    row = try
        Simbad.query_tap("""
            SELECT b.main_id, b.ra, b.dec, b.pmra, b.pmdec, b.plx_value, b.rvz_radvel,
                   b.sp_type, b.otype
            FROM basic AS b JOIN ident AS i ON i.oidref = b.oid
            WHERE i.id = '$esc'
            """)
    catch err
        throw(ErrorException("SIMBAD query failed for \"$name\": " *
                             first(split(sprint(showerror, err), "\n"))))
    end
    (pyis(row, pybuiltins.None) || pylen(row) == 0) &&
        error("SIMBAD returned no match for \"$name\"")

    mags = try
        magnitudes_from_simbad(name)
    catch
        Dict{String,Float64}()          # coordinates are still worth returning
    end

    return (; name,
              main_id = str(row["main_id"][0]),
              ra      = num(row["ra"][0]),
              dec     = num(row["dec"][0]),
              pmra    = num(row["pmra"][0]),
              pmdec   = num(row["pmdec"][0]),
              plx     = num(row["plx_value"][0]),
              rv      = num(row["rvz_radvel"][0]),
              sptype  = str(row["sp_type"][0]),
              otype   = str(row["otype"][0]),
              mags    = Dict{String,Float64}(b => get(mags, b, NaN) for b in SIMBAD_BANDS))
end

function meshgrid(xx::AbstractVector{<:Real}) #example: meshgrid([-N/2:N/2-1;]*δ);
    x = [j for i=xx, j=xx];
    return x,x'
end

function meshgrid(xx::Int64) #example: meshgrid(N);
    x = [j for i=1:xx, j=1:xx].-(div(xx,2)+1);
    return x,x'
end

function meshrad(xx::AbstractVector{<:Real})
    x = [j for i=xx, j=xx];
    return hypot.(x,x')
end

function meshrad(xx::Int64)
    x = [j for i=1:xx, j=1:xx].-(div(xx,2)+1);
    return hypot.(x,x')
end

function meshpol(xx::AbstractVector{<:Real})
    x = [j for i=1:xx, j=1:xx];
    return hypot.(x,x'), atan.(x', x)
end

function meshpol(xx::Int64)
    x = [j for i=1:xx, j=1:xx].-(div(xx,2)+1);
    return hypot.(x,x'), atan.(x', x)
end

function cart2pol(x,y)
    return hypot.(x,y), atan.(y, x)
end

function numgrad_1D(func;x=[], N=400, δ = 1e-6)
    if x==[]
        x = abs.(rand(Float64,N))
    else
        N = length(x)
    end
    numerical_g = zeros(length(x),length(func(x)))
    for i=1:N
        orig = x[i]
        x[i] = orig + 2*δ
        f2r = func(x)
        x[i] = orig + δ
        f1r = func(x)
        x[i] = orig - δ
        f1l = func(x)
        x[i] = orig - 2*δ
        f2l = func(x)
        numerical_g[i,:] = (f2l-f2r+8*(f1r-f1l))
        x[i] = orig
    end
    numerical_g ./= (12*δ)
    return numerical_g;
end

function checkgrad_1D(func;x=[], N=400, δ = 1e-6)
    if x==[]
        x = abs.(rand(N))
    else
        N = length(x)
    end
    numerical_g = similar(x)
    analytic_g = similar(x)
    for i=1:100
        orig = x[i]
        x[i] = orig + 2*δ
        f2r = func(x,analytic_g)
        x[i] = orig + δ
        f1r = func(x,analytic_g)
        x[i] = orig - δ
        f1l = func(x,analytic_g)
        x[i] = orig - 2*δ
        f2l = func(x,analytic_g)
        numerical_g[i] = (f2l-f2r+8*(f1r-f1l))
        x[i] = orig
    end
    numerical_g ./= (12*δ)
    fref = func(x,analytic_g)
    numerator = norm(analytic_g-numerical_g)
    denominator = norm(analytic_g) + norm(numerical_g)
    difference = numerator / denominator
    print(difference, "\n")
    return (numerical_g, analytic_g);
end

function checkgrad_2D(func;x=[], N=36, M= 25, δ = 1e-6)
    if x==[]
        x = abs.(rand(N,M))
    else
        N,M = size(x)
    end
    numerical_g = similar(x)
    analytic_g = vec(similar(x))
    for i=1:N
        for j=1:M
            orig = x[i,j]
            x[i,j] = orig + δ
            f1r = func(x,analytic_g)
            x[i,j] = orig - δ
            f1l = func(x,analytic_g)
            x[i,j] = orig + 2*δ
            f2r = func(x,analytic_g)
            x[i,j] = orig - 2*δ
            f2l = func(x,analytic_g)
            x[i,j] = orig + 3*δ
            f3r = func(x,analytic_g)
            x[i,j] = orig - 3*δ
            f3l = func(x,analytic_g)
            numerical_g[i,j] = (f3r-f3l+9*(f2l-f2r)+45*(f1r-f1l))
            x[i,j] = orig
        end
    end
    numerical_g = vec(numerical_g)/(60*δ)
    fref = func(x,analytic_g)
    numerator = norm(analytic_g-numerical_g)
    denominator = norm(analytic_g) + norm(numerical_g)
    difference = numerator / denominator
    print(difference, "\n")
    return (numerical_g, analytic_g);
end
