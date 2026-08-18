#
# CHARA observation planning GUI (TclTk + embedded matplotlib)
#
include("../src/OITOOLS.jl")
using .OITOOLS, Dates
using PythonPlot
using TclTk

# ── CHARA combiner choices ────────────────────────────────────────────────────
const CHARA_COMBINERS = ["MIRCX", "MYSTIC", "SPICA"]
const CHARA_TEL_NAMES = ["S1", "S2", "E1", "E2", "W1", "W2"]
const CHARA_REF_INDEX = 6  # W2 is the default reference

# ── SIMBAD cache ──────────────────────────────────────────────────────────────
# Stores (ra, dec, mags_dict) per target
const simbad_cache = Dict{String, Tuple{Float64, Float64, Dict{String,Float64}}}()
const MAG_BANDS = ["V", "J", "H", "K", "L", "M", "N"]

# ── Temp file for plot PNG ────────────────────────────────────────────────────
const PLOT_FILE = joinpath(tempdir(), "chara_plan_gui.png")
const PLOT_FIGSIZE = (10, 4)

# ── Global state ──────────────────────────────────────────────────────────────
# We keep references to the Tk image and label so we can update them
const gui_state = Dict{Symbol, Any}()

function set_status!(interp, msg)
    interp["status_text"] = msg
end

function get_obsdate(interp)
    return DateTime(Date(convert(String, interp["date_var"]), dateformat"yyyy-mm-dd"))
end

function set_obsdate!(interp, d::Date)
    interp["date_var"] = Dates.format(d, dateformat"yyyy-mm-dd")
end

function step_date!(interp, delta::Int)
    d = Date(convert(String, interp["date_var"]), dateformat"yyyy-mm-dd")
    set_obsdate!(interp, d + Dates.Day(delta))
end

function get_facility()
    read_facility_file("CHARA")
end

function get_elevation(interp, varname, default)
    try
        v = parse(Float64, convert(String, interp[varname]))
        return (0.0 <= v <= 90.0) ? v : default
    catch
        return default
    end
end

# ── Render empty night ────────────────────────────────────────────────────────
function render_empty_night(interp)
    try
        obsdate = get_obsdate(interp)
        facility = get_facility()
        ioff()
        empty_night(facility, obsdate; figsize=PLOT_FIGSIZE, savefile=PLOT_FILE)
        ion()
        update_plot_image!(interp)
        interp["pop_text"] = ""
        interp["mag_text"] = ""
        set_status!(interp, "Night of $(Date(obsdate)) — enter a target and click Plan")
    catch e
        set_status!(interp, "Error: $(sprint(showerror, e))")
    end
end

# ── Update the TkPhoto from the PNG file ──────────────────────────────────────
function update_plot_image!(interp)
    imgname = gui_state[:plot_image_name]
    interp(:image, :create, :photo, imgname, "-file", PLOT_FILE)
end

# ── Date changed callback ────────────────────────────────────────────────────
function on_date_changed(interp::TclInterp, args::TclObj)
    # If we have a cached target, recompute full plan; otherwise just redraw empty night
    targetname = String(strip(convert(String, interp["target_var"])))
    key = lowercase(targetname)
    if !isempty(targetname) && haskey(simbad_cache, key)
        run_plan(interp, args)
    else
        render_empty_night(interp)
    end
end

# ── Main plan callback ───────────────────────────────────────────────────────
function run_plan(interp::TclInterp, args::TclObj)
    try
        targetname = String(strip(convert(String, interp["target_var"])))
        if isempty(targetname)
            render_empty_night(interp)
            return
        end

        obsdate = get_obsdate(interp)
        combiner = convert(String, interp["combiner_var"])
        show_alt = convert(Int, interp["details_var"]) == 1

        # Telescope config from checkboxes
        config = Int[]
        for i in 1:6
            val = convert(Int, interp["tel$(i)_var"])
            push!(config, val)
        end
        active = findall(config .> 0)
        if length(active) < 2
            set_status!(interp, "Error: select at least 2 telescopes")
            return
        end
        if config[CHARA_REF_INDEX] > 0
            config[CHARA_REF_INDEX] = 2
        else
            config[active[end]] = 2
        end

        # Resolve target (cached)
        key = lowercase(targetname)
        if haskey(simbad_cache, key)
            ra, dec, mags = simbad_cache[key]
        else
            set_status!(interp, "Querying SIMBAD for $(targetname)...")
            ra, dec = ra_dec_from_simbad(targetname)   # decimal degrees
            mags = magnitudes_from_simbad(targetname)
            simbad_cache[key] = (ra, dec, mags)
        end
        ra_h = ra / 15.0
        interp["radec_text"] = string("RA ", round(ra_h, digits=4), "h  Dec ", round(dec, digits=4), "°")
        # Display magnitudes
        mag_parts = [isnan(mags[b]) ? "" : "$(b)=$(round(mags[b], digits=1))"
                     for b in MAG_BANDS]
        interp["mag_text"] = join(filter(!isempty, mag_parts), "  ")

        set_status!(interp, "Computing observability...")

        alt_limit = get_elevation(interp, "altlimit_var", 45.0)
        alt_max   = get_elevation(interp, "altmax_var", 80.0)

        facility = get_facility()
        obs = night_observability(facility, ra, dec, obsdate;
                                  alt_limit=alt_limit, alt_max=alt_max, moon_min_sep=30.0)

        if isempty(obs.good_alt)
            interp["pop_text"] = ""
            # Still show the empty night
            ioff()
            empty_night(facility, obsdate; figsize=PLOT_FIGSIZE, savefile=PLOT_FILE)
            ion()
            update_plot_image!(interp)
            set_status!(interp, "Target never rises above $(round(Int,alt_limit))° this night")
            return
        end

        set_status!(interp, "Searching best POP...")
        results = best_pop(facility, dec, obs.ha, config; n_best=1, min_minutes=1)

        if isempty(results)
            interp["pop_text"] = "No valid POP"
            ioff()
            # Show night with altitude but no delay
            gantt_onenight(targetname, obsdate, obs.lst, obs.lst_midnight,
                           obs.az, obs.alt, obs.good_alt, Int[];
                           good_twilight=obs.good_twilight,
                           figsize=PLOT_FIGSIZE, savefile=PLOT_FILE, show_alt=show_alt)
            ion()
            update_plot_image!(interp)
            set_status!(interp, "No valid POP configuration found")
            return
        end

        pop = results[1].pop
        active_now = findall(config .> 0)
        pop_str = join(["$(CHARA_TEL_NAMES[i]):POP$(pop[i])" for i in active_now], "  ")
        moon_str = round(minimum(obs.moon_sep), digits=1)
        fli_str  = round(obs.moon_fli * 100, digits=1)

        interp["pop_text"] = pop_str

        # Render plot to file
        ioff()
        obs_plan(targetname, facility, ra, dec, obsdate, pop, config;
                 alt_limit=alt_limit, alt_max=alt_max,
                 figsize=PLOT_FIGSIZE, savefile=PLOT_FILE, show_alt=show_alt)
        ion()
        update_plot_image!(interp)

        set_status!(interp, "$(results[1].score) min in delay — Moon sep $(moon_str)°, FLI $(fli_str)%")

    catch e
        set_status!(interp, "Error: $(sprint(showerror, e))")
    end
end

# ── Build GUI ─────────────────────────────────────────────────────────────────
function chara_plan_gui()
    interp = tk_start()
    top = TkToplevel(interp)
    interp(:wm, :title, top, "AsprOI")
    interp(:wm, :resizable, top, 0, 0)

    frame = TtkFrame(top, padding=(12, 12, 12, 12))
    frame.grid(column=0, row=0, sticky="nwes")

    row = 0

    # ── Target ────────────────────────────────────────────────────────────
    TtkLabel(frame, text="Target:").grid(column=0, row=row, sticky="e", padx=5, pady=3)
    interp["target_var"] = "AZ Cyg"
    target_entry = TtkEntry(frame, width=12, textvariable="target_var")
    target_entry.grid(column=1, row=row, sticky="w", padx=5, pady=3)
    interp["radec_text"] = ""
    TtkLabel(frame, textvariable="radec_text", foreground="gray").grid(
        column=2, row=row, columnspan=3, sticky="w", padx=5, pady=3)

    row += 1

    # ── Magnitudes ─────────────────────────────────────────────────────────
    TtkLabel(frame, text="Mags:").grid(column=0, row=row, sticky="e", padx=5, pady=1)
    interp["mag_text"] = ""
    TtkLabel(frame, textvariable="mag_text", foreground="gray",
             font="TkSmallCaptionFont").grid(
        column=1, row=row, columnspan=4, sticky="w", padx=5, pady=1)

    row += 1

    # ── Date (default: 2026-06-03) ────────────────────────────────────────
    TtkLabel(frame, text="Date:").grid(column=0, row=row, sticky="e", padx=5, pady=3)

    interp["date_var"] = "2026-06-03"

    date_frame = TtkFrame(frame)
    date_frame.grid(column=1, row=row, columnspan=4, sticky="w", padx=5, pady=3)

    prev_day_cb = TclTk.Callback((i, a) -> (step_date!(i, -1); on_date_changed(i, a)))
    next_day_cb = TclTk.Callback((i, a) -> (step_date!(i, +1); on_date_changed(i, a)))

    TtkButton(date_frame, text="<", width=2, command=prev_day_cb.name).grid(
        column=0, row=0, padx=2)
    TtkEntry(date_frame, textvariable="date_var", width=12).grid(
        column=1, row=0, padx=2)
    TtkButton(date_frame, text=">", width=2, command=next_day_cb.name).grid(
        column=2, row=0, padx=2)

    row += 1

    # ── Combiner ──────────────────────────────────────────────────────────
    TtkLabel(frame, text="Combiner:").grid(column=0, row=row, sticky="e", padx=5, pady=3)
    interp["combiner_var"] = "MIRCX"
    combo = TtkCombobox(frame, textvariable="combiner_var", width=12,
                        values=join(CHARA_COMBINERS, " "), state="readonly")
    combo.grid(column=1, row=row, columnspan=2, sticky="w", padx=5, pady=3)

    row += 1

    # ── Telescopes ────────────────────────────────────────────────────────
    TtkLabel(frame, text="Telescopes:").grid(column=0, row=row, sticky="e", padx=5, pady=3)
    tel_frame = TtkFrame(frame)
    tel_frame.grid(column=1, row=row, columnspan=4, sticky="w", padx=5, pady=3)

    for i in 1:6
        interp["tel$(i)_var"] = "1"
        TtkCheckbutton(tel_frame, text=CHARA_TEL_NAMES[i],
                       variable="tel$(i)_var").grid(column=i-1, row=0, padx=3)
    end

    row += 1

    # ── Elevation limits ─────────────────────────────────────────────────
    TtkLabel(frame, text="Elevation:").grid(column=0, row=row, sticky="e", padx=5, pady=3)
    elev_frame = TtkFrame(frame)
    elev_frame.grid(column=1, row=row, columnspan=4, sticky="w", padx=5, pady=3)

    interp["altlimit_var"] = "45"
    interp["altmax_var"] = "80"
    TtkLabel(elev_frame, text="Min:").grid(column=0, row=0, padx=(0,2))
    altlimit_entry = TtkEntry(elev_frame, textvariable="altlimit_var", width=4)
    altlimit_entry.grid(column=1, row=0, padx=2)
    TtkLabel(elev_frame, text="Max:").grid(column=2, row=0, padx=(8,2))
    altmax_entry = TtkEntry(elev_frame, textvariable="altmax_var", width=4)
    altmax_entry.grid(column=3, row=0, padx=2)

    row += 1

    # ── Plan button ───────────────────────────────────────────────────────
    plan_cb = TclTk.Callback(run_plan)
    TtkButton(frame, text="Plan!", command=plan_cb.name).grid(
        column=0, row=row, columnspan=5, pady=10)

    row += 1

    # ── Best POP display ──────────────────────────────────────────────────
    TtkLabel(frame, text="Best POP:").grid(column=0, row=row, sticky="e", padx=5, pady=3)
    interp["pop_text"] = ""
    TtkLabel(frame, textvariable="pop_text", font="TkFixedFont",
             foreground="blue").grid(
        column=1, row=row, columnspan=4, sticky="w", padx=5, pady=3)

    row += 1

    # ── Plot image ────────────────────────────────────────────────────────
    imgname = "plotimg"
    interp(:image, :create, :photo, imgname, "-width", 1, "-height", 1)
    gui_state[:plot_image_name] = imgname
    plot_label = TkLabel(frame, image=imgname, borderwidth=0)
    plot_label.grid(column=0, row=row, columnspan=5, padx=5, pady=5)

    row += 1

    # ── Details checkbox ──────────────────────────────────────────────────
    interp["details_var"] = "0"
    details_cb = TclTk.Callback(on_date_changed)
    TtkCheckbutton(frame, text="Details", variable="details_var",
                   command=details_cb.name).grid(
        column=0, row=row, columnspan=5, sticky="w", padx=5, pady=3)

    row += 1

    # ── Status bar ────────────────────────────────────────────────────────
    interp["status_text"] = "Ready"
    TtkLabel(frame, textvariable="status_text", foreground="gray").grid(
        column=0, row=row, columnspan=5, sticky="we", padx=5, pady=3)

    # Bind Enter key
    bind(top, "<Return>", plan_cb.name)

    # Elevation entries trigger recalculation on Enter or focus-out
    elev_cb = TclTk.Callback(on_date_changed)
    bind(altlimit_entry, "<Return>", elev_cb.name)
    bind(altlimit_entry, "<FocusOut>", elev_cb.name)
    bind(altmax_entry, "<Return>", elev_cb.name)
    bind(altmax_entry, "<FocusOut>", elev_cb.name)

    interp(:focus, target_entry)

    # Render initial empty night
    render_empty_night(interp)

    return interp
end

chara_plan_gui()
