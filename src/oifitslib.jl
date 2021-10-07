# Julia's command line wrapper for John Young's oifitslib library
# Library available here https://github.com/jsy1001/oifitslib

using Glob
function oifits_check(;infile::String="")
    if infile==""
        println("Syntax: ")
        println(raw"oifits_check(infile=\"./dir/input.oifits\")")
    else
        command = string("oifits-check ", infile)
        println("$command")
        run(`sh -c $command`)
    end
end

function oifits_merge(;outfile::String="", infiles::Union{String, Vector{String}}="")
    if typeof(infiles)==Vector{String}
        infiles_str = join(infiles, " ")
    else
        infiles_str = join(glob(infiles)," ")
    end
    if (infiles=="")||(outfile=="")
        println("Syntax for oifits_merge: ")
        println(raw"oifits_merge(outfile=\"./pathto/out.oifits\", infiles=\"./pathto/input1.oifits /pathto/intput2.oifits\"])")
        println(raw"oifits_merge(outfile=\"./pathto/out.oifits\", infiles=\"./pathto/input*.oifits\"")
    else
        # automatic removal if merged file already exists
        if isfile(outfile)
            println("outfile already exists -- removing");
            rm(outfile)
        end
        command = string("oifits-merge ", outfile, " ", infiles)
        println("$command")
        run(`sh -c $command`)
    end
end


function oifits_filter(;outfile::String="", infile::String="", arrname::String="", insname::String="", corrname::String="",target_id::Int=-1,
    mjd_min::Float64=0.0, mjd_max::Float64=1e99, wave_min::Float64=0.0, wave_max::Float64=1e99, bas_min::Float64=0.0, bas_max::Float64=1e99,
    uvrad_min::Float64=0.0, uvrad_max::Float64=1e99, snr_min::Float64=0.0, snr_max::Float64=1e99,
    accept_vis::Bool=true, accept_vis2::Bool=true, accept_t3amp::Bool=true, accept_t3phi::Bool=true, accept_flux::Bool=true, accept_flagged::Bool=true)
    if (infile=="")||(outfile=="")
        println("Syntax for oifits_filter: ")
        println(raw"oifits_filter(outfile=\"./pathto/out.oifits\", infile=\"./pathto/input1.oifits\"])")
    else
        if infile == outfile
            error("oifits_filter error: do not use the same infile and outfile")
        end
        if !isfile(infile)
            error("oifits_filter error: infile does not exist")
        end
        options = " -o " #in OITOOLS we use clobber by default
        options = string(options, "--accept-vis=$(Int(accept_vis)) --accept-vis2=$(Int(accept_vis2)) --accept-t3amp=$(Int(accept_t3amp)) --accept-t3phi=$(Int(accept_t3phi)) --accept-flux=$(Int(accept_flux)) --accept-flagged=$(Int(accept_flagged))")
        if arrname !=""
            options = string(options, " --arrname=$arrname")
        end
        if insname !=""
            options = string(options, " --insname=$insname")
        end
        if corrname !=""
            options = string(options, " --corrname=$corrname")
        end
        if target_id !=-1
            options = string(options, " --target-id=$(target_id)")
        end

        if (mjd_min !=0) || (mjd_max !=199)
            options = string(options, " --mjd-min=$mjd_min --mjd-max=$mjd_max")
        end

        if (wave_min !=0) || (wave_max !=199)
            options = string(options, " --wave-min=$wave_min --wave-max=$wave_max")
        end

        if (bas_min !=0) || (bas_max !=199)
            options = string(options, " --bas-min=$bas_min --bas-max=$bas_max")
        end

        if (uvrad_min !=0) || (uvrad_min !=199)
            options = string(options, " --uvrad-min=$uvrad_min --uvrad-max=$uvrad_max")
        end

        if (snr_min !=0) || (snr_max !=199)
            options = string(options, " --snr-min=$snr_min --snr-max=$snr_max")
        end

        command = string("oifits-filter ", options, " ", infile, " ", outfile)
        println("$command")
        run(`sh -c $command`)
    end
end
