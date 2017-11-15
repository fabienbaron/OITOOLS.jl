#
# oifile.jl --
#
# Implement reading/writing of OI-FITS data from/to FITS files.
#
#------------------------------------------------------------------------------
#
# This file is part of OIFITS.jl which is licensed under the MIT "Expat"
# License:
#
# Copyright (C) 2015-2017: Éric Thiébaut.
#
#------------------------------------------------------------------------------

using FITSIO
using FITSIO.Libcfitsio

# Read a column from a table (this is low-level API, it is expected that the
# FITS file is open, this must be done by the caller with fits_assert_open).
function read_column(ff::FITSFile, colnum::Integer)
    # Get the type and the dimensions of the data stored in the column.
    (typecode, repcnt, width) = fits_get_eqcoltype(ff, colnum)
    dims = fits_read_tdim(ff, colnum)
    nrows = fits_get_num_rows(ff)

    # Allocate the array and read the column contents.
    T = cfitsio_datatype(typecode)
    if T <: AbstractString
        # Column contains an array of strings.  Strip the leading dimension
        # which is the maximum length of each strings.  On return trailing
        # spaces are removed (they are insignificant according to the FITS
        # norm).
        if length(dims) == 1
            dims = nrows
        else
            dims[1:end-1] = dims[2:end]
            dims[end] = nrows
        end
        data = Array{T}(dims...)
        fits_read_col(ff, colnum, 1, 1, data)
        return map(rstrip, data)
    elseif T == Void
        error("unsupported column data")
    else
        # Column contains numerical data.
        if length(dims) == 1 && dims[1] == 1
            # Result will be a simple vector.
            dims = nrows
        else
            # Result will be a multi-dimensional array.
            push!(dims, nrows)
        end
        data = Array{T}(dims...)
        fits_read_col(ff, colnum, 1, 1, data)
        return data
    end
end

# Get the type of the data-block.
function get_dbtype(hdr::FITSHeader)
    if get_hdutype(hdr) == :binary_table
        extname = fixname(get_string(hdr, "EXTNAME", ""))
        if startswith(extname, "OI_")
            return Symbol(replace(extname, r"[^A-Z0-9_]", '_'))
        end
    end
    :unknown
end

const _COMMENT = Set(["HISTORY", "COMMENT"])

function get_value(hdr::FITSHeader, key::AbstractString)
    haskey(hdr, key) || error("missing FITS keyword \"$key\"")
    hdr[key]
end

function get_value(hdr::FITSHeader, key::AbstractString, def)
    haskey(hdr, key) ? hdr[key] : def
end

function get_comment(hdr::FITSHeader, key::AbstractString)
    haskey(hdr, key) || error("missing FITS keyword \"$key\"")
    FITSIO.get_comment(hdr, key)
end

function get_comment(hdr::FITSHeader, key::AbstractString, def::AbstractString)
    haskey(hdr, key) ? get_comment(hdr, key) : def
end

for (fn, T, S) in ((:get_integer, Integer,        Int),
                   (:get_real,    Real,           Float64),
                   (:get_logical, Bool,           Bool),
                   (:get_string,  AbstractString, AbstractString))
    if S == T
        @eval begin
            function $fn(hdr::FITSHeader, key::AbstractString, def::$T)
                val = haskey(hdr, key) ? hdr[key] : def
                isa(val, $T) || error("bad type for FITS keyword \"$key\"")
                return val
            end
            function $fn(hdr::FITSHeader, key::AbstractString)
                haskey(hdr, key) || error("missing FITS keyword \"$key\"")
                val = hdr[key]
                isa(val, $T) || error("bad type for FITS keyword \"$key\"")
                return val
            end
        end
    else
        @eval begin
            function $fn(hdr::FITSHeader, key::AbstractString, def::$T)
                val = haskey(hdr, key) ? hdr[key] : def
                isa(val, $T) || error("bad type for FITS keyword \"$key\"")
                return convert($S, val)
            end
            function $fn(hdr::FITSHeader, key::AbstractString)
                haskey(hdr, key) || error("missing FITS keyword \"$key\"")
                val = hdr[key]
                isa(val, $T) || error("bad type for FITS keyword \"$key\"")
                return convert($S, val)
            end
        end
    end
end

# Returns invalid result if not a valid OI-FITS data-block.
# Unless quiet is true, print warn message.
function check_datablock(hdr::FITSHeader; quiet::Bool=false)
    # Values returned in case of error.
    dbname = ""
    dbrevn = -1
    dbdefn = nothing

    # Use a while loop to break out whenever an error occurs.
    while get_hdutype(hdr) == :binary_table
        # Get extension name.
        extname = get_value(hdr, "EXTNAME", nothing)
        if ! isa(extname, AbstractString)
            quiet || warn(extname == nothing ? "missing keyword \"EXTNAME\""
                                             : "EXTNAME value is not a string")
            break
        end
        extname = fixname(extname)
        startswith(extname, "OI_") || break
        dbname = extname
        if ! haskey(_DATABLOCKS, dbname)
            quiet || warn("unknown OI-FITS data-block \"$extname\"")
            break
        end

        # Get revision number.
        revn = 1#get_value(hdr, "OI_REVN", nothing)
        if ! isa(revn, Integer)
            quiet || warn(revn == nothing ? "missing keyword \"OI_REVN\""
                                          : "\"OI_REVN\" value is not an integer")
            break
        end
        dbrevn = revn
        if dbrevn <= 0
            quiet || warn("invalid \"OI_REVN\" value ($dbrevn)")
            break
        end
        if dbrevn > length(_FORMATS)
            quiet || warn("unsupported \"OI_REVN\" value ($dbrevn)")
            break
        end
        if ! haskey(_FORMATS[dbrevn], dbname)
            quiet || warn("unknown OI-FITS data-block \"$extname\"")
        end
        dbdefn = _FORMATS[dbrevn][dbname]
        break
    end
    return (dbname, dbrevn, dbdefn)
end

function hash_column_names(hdr::FITSHeader)
    columns = Dict{Name,Int}()
    hdutype = get_hdutype(hdr)
    if hdutype == :binary_table || hdutype == :ascii_table
        ncols = get_integer(hdr, "TFIELDS", 0)
        for k in 1:ncols
            ttype = get_string(hdr, "TTYPE$k")
            columns[fixname(ttype)] = k
        end
    end
    columns
end

function get_file_handle(hdu::HDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    return hdu.fitsfile
end

function read_datablock(hdu::HDU; quiet::Bool=false)
    read_datablock(hdu, read_header(hdu), quiet=quiet)
end

function read_datablock(hdu::HDU, hdr::FITSHeader; quiet::Bool=false)
    ff = get_file_handle(hdu)
    (dbtype, revn, defn) = check_datablock(hdr, quiet=quiet)
    defn == nothing && return nothing
    columns = hash_column_names(hdr)
    nerrs = 0
    data = Dict{Symbol,Any}(:revn => revn)
    for field in defn.fields
        spec = defn.spec[field]
        name = spec.name
        if spec.keyword
            value = get_value(hdr, name, nothing)
            if value == nothing
                warn("missing keyword \"$name\" in OI-FITS $dbtype data-block")
                nerrs += 1
            else
                data[field] = value
            end
        else
            colnum = get(columns, name, 0)
            if colnum < 1
                warn("missing column \"$name\" in OI-FITS $dbtype data-block")
                nerrs += 1
            else
                data[field] = read_column(ff, colnum)
            end
        end
    end
    nerrs > 0 && error("bad OI-FITS $dbtype data-block")
    return build_datablock(dbtype, revn, data)
end

function load(filename::AbstractString; quiet::Bool=true)
    return load(FITS(filename, "r"), quiet=quiet)
end

function load(f::FITS; quiet::Bool=true)
    master = new_master()

    # Read all contents, skipping first HDU.
    for hdu in 2:length(f)
        db = read_datablock(f[hdu], quiet=quiet)
        if db == nothing
            quiet || println("skipping HDU $hdu (no OI-FITS data)")
            continue
        end
        dbname = _EXTNAMES[typeof(db)]
        quiet || println("reading OI-FITS $dbname in HDU $hdu")
        attach!(master, db)
    end
    update!(master)
    return master
end
