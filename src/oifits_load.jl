# Opening an OIFITS file: `load_oifits`, the one entry point every reader goes through
# (`readoifits`, `list_oifits_targets`, `write_oifits`' template, the validator, the GUI).
#
# OIFITS.jl reads a file whole and stops at the first table it cannot read. Files in the wild
# break the standard in ways that leave most of their data readable (PLAN_C14.md): a table
# declaring a revision whose columns it lacks, or a revision that does not exist; GRAVITY's
# OI_FLUX naming its data `FLUX`; an OI_VIS with amplitudes only; rows whose TARGET_ID is not in
# OI_TARGET; a FITS file inside a ZIP archive. So the file is read here table by table, through
# OIFITS.jl's own per-table reader, and a table that fails is read the closest way that works --
# that table only -- and the way is reported.

import AstroFITS, OIFITS, p7zip_jll

const _OI_DEPS  = ("OI_TARGET", "OI_ARRAY", "OI_WAVELENGTH", "OI_CORR")
const _OI_BLOCK = Dict("OI_TARGET" => OIFITS.OI_TARGET, "OI_ARRAY" => OIFITS.OI_ARRAY,
                       "OI_WAVELENGTH" => OIFITS.OI_WAVELENGTH, "OI_CORR" => OIFITS.OI_CORR,
                       "OI_VIS" => OIFITS.OI_VIS, "OI_VIS2" => OIFITS.OI_VIS2, "OI_T3" => OIFITS.OI_T3,
                       "OI_FLUX" => OIFITS.OI_FLUX, "OI_INSPOL" => OIFITS.OI_INSPOL)
# The highest revision OIFITS.jl defines for each table.
const _OI_MAXREV = Dict("OI_TARGET" => 2, "OI_ARRAY" => 2, "OI_WAVELENGTH" => 2, "OI_CORR" => 1,
                        "OI_VIS" => 2, "OI_VIS2" => 2, "OI_T3" => 2, "OI_FLUX" => 1, "OI_INSPOL" => 1)
# A column under another name in some pipelines' files (GRAVITY's OI_FLUX).
const _OI_ALIAS = Dict("FLUXDATA" => "FLUX")
# The two halves of an amplitude/phase pair: a table carrying one half only is read with the other
# NaN, which `readoifits` holds as "not measured".
const _OI_PAIR  = Dict("VISAMP" => "VISPHI", "VISAMPERR" => "VISPHIERR",
                       "VISPHI" => "VISAMP", "VISPHIERR" => "VISAMPERR",
                       "T3AMP" => "T3PHI", "T3AMPERR" => "T3PHIERR",
                       "T3PHI" => "T3AMP", "T3PHIERR" => "T3AMPERR")

"""
    load_oifits(path; notes = String[]) -> OIFITS.OIDataSet

Read an OIFITS file as an `OIDataSet`, table by table, accepting what the standard's readers
refuse when the data are still readable. Each table is read at the revision it declares; if that
fails, at the other revisions OIFITS.jl defines (highest first); if none reads it, with a missing
column tolerated -- `FLUX` read as OI_FLUX's `FLUXDATA`, one half of an amplitude/phase pair left
NaN. Rows whose TARGET_ID is not in the file's OI_TARGET are dropped (a table left empty is
skipped), as PMOIRED does. A ZIP archive holding one FITS file is read through that file.

Each such step is described in `notes`, one line per table. A table that cannot be read at all
throws, naming the first reason.
"""
function load_oifits(filename::AbstractString; notes::Vector{String} = String[])
    path, tmp = filename, nothing
    if _is_zip(filename)
        tmp, member = _unzip_fits(filename)
        path = tmp
        push!(notes, "a ZIP archive; read its member $(member)")
    end
    try
        AstroFITS.FitsFile(path) do f
            ds = OIFITS.OIDataSet()
            ids = Set{Int}()            # the file's own TARGET_IDs; `push!` renumbers them
            for pass in 1:2, i in 2:length(f)
                hdu = f[i]
                hdu isa AstroFITS.FitsTableHDU || continue
                ext = something(hdu.extname, "")
                haskey(_OI_BLOCK, ext) || continue
                (ext in _OI_DEPS) == (pass == 1) || continue
                lab = _table_label(ext, hdu)
                db = _read_oiblock(ext, hdu, lab, notes)
                if pass == 2
                    db = _known_targets(db, ids, lab, notes)
                    db === nothing && continue
                end
                ext == "OI_TARGET" && union!(ids, (e.target_id for e in db.list))
                push!(ds, db)
            end
            ds
        end
    finally
        tmp === nothing || rm(tmp; force = true)
    end
end

_table_label(ext, hdu) = (ins = get(String, hdu, "INSNAME", ""); isempty(ins) ? ext : "$(ext) $(ins)")

_why(e) = e isa OIFITS.MissingColumn ? "no $(e.col) column" :
          e isa MethodError ? "a revision OIFITS does not define" : sprint(showerror, e)

# One table: its declared revision, then the others OIFITS.jl defines, then field by field.
function _read_oiblock(ext, hdu, lab, notes)
    T = _OI_BLOCK[ext]
    first_err = try
        return read(T, hdu)
    catch e
        e
    end
    declared = get(Int, hdu, "OI_REVN", 0)
    for r in _OI_MAXREV[ext]:-1:1
        r == declared && continue
        try
            db = read(T, hdu; hack_revn = r)
            push!(notes, "$(lab): read as revision $(r) (declared $(declared): $(_why(first_err)))")
            return db
        catch
        end
    end
    T === OIFITS.OI_TARGET && throw(first_err)      # read by its own reader, not field by field
    for r in _OI_MAXREV[ext]:-1:1
        try
            return _read_tolerant(T, hdu, r, lab, notes)
        catch
        end
    end
    throw(first_err)
end

# A table read field by field from its revision-`r` format, as OIFITS.jl's own reader does, with a
# column under its alias, and a missing half of an amplitude/phase pair set to NaN.
function _read_tolerant(T, hdu, r, lab, notes)
    db = T(undef); db.revn = r
    cols = hdu.column_names
    spec_of = Dict{String,Any}(); halves = String[]; said = String[]
    for spec in OIFITS.get_format(db)
        S = spec.symb; spec_of[spec.name] = spec
        if spec.rank == 0
            S === :revn && continue
            try
                setfield!(db, S, get(fieldtype(T, S), hdu, spec.name))
            catch e
                (spec.optional && e isa KeyError) || rethrow()
            end
            continue
        end
        name = spec.name
        if !(name in cols) && haskey(_OI_ALIAS, name) && _OI_ALIAS[name] in cols
            push!(said, "$(_OI_ALIAS[name]) read as $(name)")
            name = _OI_ALIAS[name]
        end
        if name in cols
            setfield!(db, S, OIFITS.read_column(fieldtype(T, S), hdu, name))
        elseif spec.optional
        elseif haskey(_OI_PAIR, name) && _OI_PAIR[name] in cols
            push!(halves, name)
        else
            throw(OIFITS.MissingColumn(name, something(hdu.extname, "")))
        end
    end
    for name in halves
        setfield!(db, spec_of[name].symb, fill(NaN, size(getfield(db, spec_of[_OI_PAIR[name]].symb))))
    end
    isempty(halves) || push!(said, "no " * join(halves, "/") * " column (left missing)")
    push!(notes, "$(lab): read as revision $(r), " * join(said, ", "))
    return db
end

# The rows of a data table whose TARGET_ID is in OI_TARGET; `nothing` when there are none. Every
# per-row field is subset along its last axis, from the table's own format.
function _known_targets(db, ids, lab, notes)
    keep = findall(in(ids), db.target_id)
    n = length(db.target_id)
    length(keep) == n && return db
    if isempty(keep)
        push!(notes, "$(lab): skipped, no row has a TARGET_ID listed in OI_TARGET")
        return nothing
    end
    for spec in OIFITS.get_format(db)
        (spec.rank == 0 || !isdefined(db, spec.symb)) && continue
        x = getfield(db, spec.symb)
        setfield!(db, spec.symb, x isa AbstractVector ? x[keep] :
                                 x[ntuple(_ -> Colon(), ndims(x) - 1)..., keep])
    end
    push!(notes, "$(lab): $(n - length(keep)) of $(n) rows dropped, their TARGET_ID is not in OI_TARGET")
    return db
end

# A ZIP archive starts with "PK\x03\x04"; some pipelines deliver their .fits files zipped.
_is_zip(path) = open(io -> (b = read(io, 4); length(b) == 4 && b == UInt8['P', 'K', 3, 4]), path)

"""
    _unzip_fits(path) -> (tmpfile, member)

Extract the one FITS file a ZIP archive holds to a temporary file (the caller removes it), with
the `7z` that Julia ships (`p7zip_jll`). An archive with no FITS file, or several, is refused,
naming them.
"""
function _unzip_fits(path::AbstractString)
    z = p7zip_jll.p7zip()
    lst = try
        read(pipeline(`$z l -slt -ba $path`; stderr = devnull), String)
    catch
        throw(ArgumentError("$(path) starts like a ZIP archive but cannot be read as one"))
    end
    members = [String(strip(l[8:end])) for l in split(lst, '\n') if startswith(l, "Path = ")]
    fits = filter(m -> occursin(r"\.(oifits|fits|fit)$"i, m), members)
    length(fits) == 1 || throw(ArgumentError("$(path) is a ZIP archive holding " *
        (isempty(fits) ? "no FITS file" : "$(length(fits)) FITS files ($(join(fits, ", ")))") *
        "; unzip it and open the file wanted"))
    tmp = tempname() * ".fits"
    open(io -> run(pipeline(`$z x -so $path $(fits[1])`; stdout = io)), tmp, "w")
    return tmp, fits[1]
end
