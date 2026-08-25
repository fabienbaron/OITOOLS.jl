# simbad.jl
#
# SIMBAD, in one HTTP request per target, over the TAP service CDS publishes for exactly this
# purpose (https://simbad.cds.unistra.fr/simbad/sim-tap).
#
# One request per target
# ----------------------
# SIMBAD asks clients to send a single ADQL query and let the join happen server-side, which
# is what `_TARGET_QUERY` does: `basic` for the astrometry, `flux` for the magnitudes, one
# result set. Fetching coordinates and photometry separately would double the load on a public
# service for a join the server does better.
#
# The transport is `Downloads`, a standard library, so this costs the package nothing: no
# Python import, no conda package, and no astropy table to unwrap on the way back.
#
# Why CSV rather than VOTable or JSON
# -----------------------------------
# TAP will return any of them. VOTable needs an XML parser and JSON needs a JSON parser, both
# of which would be a new dependency for a payload that is a dozen rows of a dozen columns.
# CSV needs `simbad_parse_csv` below, which is thirty lines and handles the one thing that actually
# occurs in SIMBAD output: quoted fields containing commas (`sp_type`, `main_id`).
#
# Errors arrive as XML whatever FORMAT asks for, which is how `simbad_tap` tells a rejected
# query from a result.

using Downloads

"""
    SIMBAD_TAP_URL

The TAP endpoint queries are sent to. Overridable per call through `simbad_tap`'s `url`
keyword, which is what the mirror at `simbad.u-strasbg.fr` is for.
"""
const SIMBAD_TAP_URL = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"

"""
    SIMBAD_BANDS

Photometric bands `simbad_target` reports, in the order a panel should show them.

B through N. The near-infrared ones are what an interferometric observation is actually made
in, and V is what the adaptive optics guides on, so both matter and for different reasons.
"""
const SIMBAD_BANDS = ["B", "V", "R", "I", "J", "H", "K", "L", "M", "N"]

# Everything an observation needs, in one result set. `basic` carries the astrometry and one
# row per object; `flux` carries one row per photometric band, so the join multiplies the
# astrometry across as many rows as the target has magnitudes. That redundancy is the price of
# a single request and it is a few hundred bytes.
#
# LEFT OUTER, because a target with no published photometry still has coordinates, and an
# inner join would report it as "no such target".
#
# Aliased columns throughout: the header CDS returns for a bare `b.plx_value` is not worth
# guessing at, and `ra`/`dec` are close enough to ADQL keywords to be worth renaming.
const _TARGET_QUERY = """
SELECT b.main_id AS main_id, b.ra AS ra_deg, b.dec AS dec_deg,
       b.pmra AS pmra, b.pmdec AS pmdec, b.plx_value AS plx, b.rvz_radvel AS rv,
       b.sp_type AS sptype, b.otype AS otype,
       f.filter AS band, f.flux AS mag
FROM ident AS i
JOIN basic AS b ON b.oid = i.oidref
LEFT OUTER JOIN flux AS f ON f.oidref = b.oid
WHERE i.id = '__NAME__'"""

"""
    simbad_percent_encode(s) -> String

Percent-encode `s` for a URL query string, unreserved set `A-Z a-z 0-9 - _ . ~`.

Written out rather than taken from URIs.jl: this is the only URL OITOOLS builds, and a
dependency for one function that has a fixed definition in RFC 3986 is not worth it.
"""
function simbad_percent_encode(s::AbstractString)
    out = IOBuffer()
    for b in codeunits(String(s))
        if (UInt8('A') <= b <= UInt8('Z')) || (UInt8('a') <= b <= UInt8('z')) ||
           (UInt8('0') <= b <= UInt8('9')) || b in (UInt8('-'), UInt8('_'), UInt8('.'), UInt8('~'))
            write(out, b)
        else
            print(out, '%', uppercase(string(b, base = 16, pad = 2)))
        end
    end
    return String(take!(out))
end

"""
    simbad_parse_csv(text) -> Vector{Vector{String}}

Split RFC 4180 CSV into rows of fields. Handles quoted fields, doubled quotes inside them,
embedded commas and newlines, and both line endings. Blank lines are dropped.

Every field comes back as a `String`, including the empty ones: a TAP result distinguishes
"no value" only by an empty field, and turning that into `missing` here would lose the column
alignment that [`simbad_record`](@ref) relies on.
"""
function simbad_parse_csv(text::AbstractString)
    rows = Vector{Vector{String}}()
    row  = String[]
    buf  = IOBuffer()
    inq  = false
    endrow() = (push!(row, String(take!(buf))); push!(rows, row); row = String[])

    i, n = firstindex(text), lastindex(text)
    while i <= n
        c = text[i]
        i = nextind(text, i)
        if inq
            if c != '"'
                write(buf, c)
            elseif i <= n && text[i] == '"'      # "" inside a quoted field is one quote
                write(buf, '"'); i = nextind(text, i)
            else
                inq = false
            end
        elseif c == '"'
            inq = true
        elseif c == ','
            push!(row, String(take!(buf)))
        elseif c == '\r'
            i <= n && text[i] == '\n' && (i = nextind(text, i))
            endrow()
        elseif c == '\n'
            endrow()
        else
            write(buf, c)
        end
    end
    (position(buf) > 0 || !isempty(row)) && endrow()

    return filter(r -> !(length(r) == 1 && isempty(r[1])), rows)
end

"""
    simbad_tap(adql; url = SIMBAD_TAP_URL, timeout = 30) -> (colnames, rows)

Run one ADQL query against SIMBAD's TAP service and return its column names and rows, all as
`String`s.

```julia
cols, rows = simbad_tap("SELECT TOP 5 main_id, ra, dec FROM basic WHERE plx_value > 500")
```

This is the general escape hatch; [`simbad_target`](@ref) is the query OITOOLS actually needs.
The full schema is at <https://simbad.cds.unistra.fr/simbad/tap/tapsearch.html>.

Three failures are reported differently on purpose, because they call for different actions:
the service being unreachable (network), the service refusing the query (ADQL), and the query
succeeding with no rows (the target does not exist). Only the last is the caller's data
problem.
"""
function simbad_tap(adql::AbstractString; url::AbstractString = SIMBAD_TAP_URL,
                    timeout::Real = 30)
    full = string(url, "?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=",
                  simbad_percent_encode(strip(adql)))
    io  = IOBuffer()
    res = Downloads.request(full; output = io, timeout = Float64(timeout), throw = false)
    body = String(take!(io))

    if res isa Downloads.RequestError
        throw(ErrorException("SIMBAD is unreachable ($url): $(res.message)"))
    end

    # A rejected query comes back as a VOTable carrying the parser's complaint, whatever FORMAT
    # asked for -- and sometimes with status 200. Leading '<' is the reliable tell.
    if startswith(lstrip(body), '<')
        m = match(r"<INFO[^>]*value=\"ERROR\"[^>]*>(.*?)</INFO>"s, body)
        detail = m === nothing ? first(strip(body), 300) : strip(m.captures[1])
        throw(ErrorException("SIMBAD rejected the query: $detail"))
    end
    if res.status >= 400
        throw(ErrorException("SIMBAD returned HTTP $(res.status): $(first(strip(body), 300))"))
    end

    return simbad_tap_parse(body)
end

"""
    simbad_tap_parse(body) -> (colnames, rows)

Split a TAP CSV response into its header and its data rows.

Separate from [`simbad_tap`](@ref) so the parsing can be exercised against a saved response
with no network involved; `test/test_simbad.jl` reads one from `test/references`.
"""
function simbad_tap_parse(body::AbstractString)
    rows = simbad_parse_csv(body)
    isempty(rows) && return (String[], Vector{String}[])
    return (String[strip(c) for c in rows[1]], rows[2:end])
end

"""
    simbad_record(name, colnames, rows) -> NamedTuple

Turn the result of the [`simbad_target`](@ref) query into the record that function returns.

Split out from the query so the parsing can be tested against a saved response instead of the
network; `test/test_simbad.jl` does exactly that.

Columns are looked up by name, not position, so adding one to `_TARGET_QUERY` cannot silently
shift the others.
"""
function simbad_record(name::AbstractString, colnames::AbstractVector,
                       rows::AbstractVector)
    isempty(rows) && error("SIMBAD returned no match for \"$name\"")
    col = Dict(lowercase(String(c)) => i for (i, c) in enumerate(colnames))

    cell(r, k) = (i = get(col, k, 0); (i == 0 || i > length(r)) ? "" : strip(r[i]))
    str(r, k)  = String(cell(r, k))
    num(r, k)  = (s = cell(r, k); isempty(s) ? NaN : something(tryparse(Float64, s), NaN))

    # SIMBAD not knowing a magnitude and SIMBAD not being asked are the same NaN here, but a
    # band that is absent from the answer must never read as 0.0 -- zero is Vega-bright.
    mags = Dict{String,Float64}(b => NaN for b in SIMBAD_BANDS)
    for r in rows
        b = uppercase(str(r, "band"))
        haskey(mags, b) || continue
        v = num(r, "mag")
        isnan(v) || (mags[b] = v)
    end

    first_row = rows[1]
    return (; name    = String(name),
              main_id = str(first_row, "main_id"),
              ra      = num(first_row, "ra_deg"),
              dec     = num(first_row, "dec_deg"),
              pmra    = num(first_row, "pmra"),
              pmdec   = num(first_row, "pmdec"),
              plx     = num(first_row, "plx"),
              rv      = num(first_row, "rv"),
              sptype  = str(first_row, "sptype"),
              otype   = str(first_row, "otype"),
              mags)
end

"""
    simbad_target(name; timeout = 30) -> NamedTuple

Everything SIMBAD holds that an observation needs, in **one** request:
`(; name, main_id, ra, dec, pmra, pmdec, plx, rv, sptype, otype, mags)`.

`ra`/`dec` are DEGREES, proper motions mas/yr, parallax mas, radial velocity km/s. Missing
values are `NaN` and a missing string is `""` — SIMBAD not knowing a target's parallax is a
fact about the target, not a failure, and it has to be distinguishable from the query failing.
`mags` always has all ten [`SIMBAD_BANDS`](@ref) as keys, `NaN` where there is no measurement.

Beyond coordinates and magnitudes this returns what ASPRO also reads, and for the same reasons:
proper motion places the target at the epoch actually being observed, parallax turns an angular
diameter into a physical one, and the spectral type is what a surface-brightness relation needs
to predict a diameter before anything is measured.

`name` is matched against SIMBAD's identifier table, so anything SIMBAD lists works —
`"Vega"`, `"alf Lyr"`, `"HD 172167"`.

```julia
t = simbad_target("Vega")
t.ra, t.dec        # 279.2347, 38.7837 degrees
t.mags["K"]        # 0.129
```

This is the one to call when more than one field is wanted: [`ra_dec_from_simbad`](@ref) and
[`magnitudes_from_simbad`](@ref) are thin wrappers around it, so asking for both costs two
requests where this costs one.
"""
function simbad_target(name::AbstractString; timeout::Real = 30)
    n = strip(String(name))
    isempty(n) && throw(ArgumentError("simbad_target needs a target name"))
    # ADQL string literals escape a quote by doubling it. Nothing else needs escaping, and the
    # value never reaches a shell.
    adql = replace(_TARGET_QUERY, "__NAME__" => replace(n, "'" => "''"))
    cols, rows = simbad_tap(adql; timeout)
    return simbad_record(n, cols, rows)
end

"""
    ra_dec_from_simbad(name) -> (ra_deg, dec_deg)

Resolve a target name through SIMBAD and return its J2000 right ascension and declination
in **decimal degrees**.

Call [`simbad_target`](@ref) instead if the magnitudes or the parallax are wanted too; this
issues its own request.
"""
function ra_dec_from_simbad(name::AbstractString)
    t = simbad_target(name)
    return t.ra, t.dec
end

"""
    magnitudes_from_simbad(name) -> Dict{String,Float64}

Photometric magnitudes for a target, keyed by [`SIMBAD_BANDS`](@ref). Bands SIMBAD has no
measurement for are `NaN`, never `0.0`.

A network or query failure throws. That distinction is the point: a target with no K magnitude
and a query that never ran are different answers, and reporting them the same way is how
"this star has no K magnitude" comes to mean "SIMBAD was down".
"""
magnitudes_from_simbad(name::AbstractString) = simbad_target(name).mags
