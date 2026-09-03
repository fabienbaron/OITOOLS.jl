# The window: Qt callback registration and `gui()`.
#
# Split from the rest of src/gui/ because everything else describes what to draw, while this
# builds the thing that draws it and hands control to Qt's event loop.

function __init__()
    # Callbacks must be registered before any QML file that calls them is loaded.
    QML.@qmlfunction(shell_open, shell_set_view, shell_select_dataset,
                     shell_dataset_names, shell_script, shell_export, shell_pick,
                     shell_console, shell_ready, shell_reset_zoom, shell_pick_text,
                     shell_observables, shell_plot_kinds, shell_grouping_noun, shell_gantt_hover, shell_reconstruct, shell_imaging_summary,
                     shell_image_defaults, shell_ft_setup, shell_chi2_breakdown,
                     shell_facilities, shell_telescopes, shell_gantt, shell_best_pops,
                     shell_shift_date, shell_sim_source, shell_model_image, shell_simbad, shell_ui_scale,
                     shell_fit_model, shell_fit_values, shell_fit_rows,
                     shell_image_colormaps, shell_image_colormap, shell_recenter_image,
                     shell_show_start_image, shell_save_image, shell_engine_output,
                     shell_job_poll, shell_job_stop, shell_job_running,
                     shell_reset_image_zoom,
                     shell_set_plot_scale, shell_set_marker_size,
                     shell_save_settings, shell_load_settings, shell_plot_scale,
                     shell_set_zoom_step,
                     shell_model_rows, shell_model_components, shell_model_inspection,
                     shell_model_chi2, shell_set_param, shell_model_warnings, shell_fit_output,
                     shell_component_kinds, shell_add_component, shell_remove_component,
                     shell_rename_component, shell_save_figure,
                     shell_expression_keywords, shell_check_expression,
                     shell_profile_params, shell_add_profile_params, shell_profile_grid,
                     shell_profile_curves,
                     shell_flux_constraint,
                     shell_default_bounds, shell_check_constraints,
                     shell_profile_templates, shell_profile_template,
                     shell_component_geometry, shell_set_component_geometry,
                     shell_ld_laws, shell_ld_law, shell_set_ld_law,
                     shell_az_modes, shell_add_az_mode, shell_remove_az_mode,
                     shell_model_residuals, shell_model_sed, shell_model_depends,
                     shell_sparco_converged,
                     shell_result_ensemble, shell_show_result,
                     shell_bin_info, shell_rebin,
                     shell_chi2_map_info, shell_free_names, shell_nested_backend, shell_sim_band,
                     shell_open_model, shell_save_model,
                     shell_import_pmoired, shell_export_pmoired,
                     # The in-window file picker: QtQuick.Dialogs leaves its window mapped on
                     # some systems, so the picker is drawn inside our own window and needs
                     # Julia for everything the toolkit would otherwise supply.
                     picker_list, picker_places, picker_join, picker_parent,
                     picker_start, picker_kind, picker_would_overwrite, picker_examples)
    return nothing
end

"""
Folder the file dialog opens in.

Next to the last dataset loaded if there is one, then `\$OITOOLSGUI_DATA_DIR`, then the
package's `demos/data`. That last is the user-facing data directory — it holds Alpha Cen A and
B, the 2004 beauty-contest set, and `BC2026/`, whose files are the only ones shipped with
differential phase and OI_FLUX and therefore the only ones on which the `diffphi` and `flux`
views show anything.

It used to open in `test/gui/data`, which is the automated tests' fixture directory: two files,
neither with an OI_FLUX table. The env var exists so those tests can still pin it.

Opening in the filesystem root or the home directory is a small tax paid on every single open,
and it makes automated clicking guesswork.
"""
function _initial_folder(session::Session)
    if !isempty(session.datasets)
        d = dirname(abspath(session.datasets[end].path))
        isdir(d) && return "file://" * d
    end
    forced = get(ENV, "OITOOLSGUI_DATA_DIR", "")
    isempty(forced) || (isdir(forced) && return "file://" * abspath(forced))
    for sub in (joinpath("demos", "data"), joinpath("test", "gui", "data"))
        p = joinpath(pkgdir(OITOOLS), sub)
        isdir(p) && return "file://" * p
    end
    return "file://" * pwd()
end

"Tab names, in the order Main.qml lists them."
const TAB_NAMES = ("exploring", "observing", "modeling", "imaging")

"The bare noun for each tab, accepted as well: nobody types the gerund from memory."
const TAB_ALIASES = ("explore", "observe", "model", "image")

"""
Which perspective the window opens on, from `\$OITOOLSGUI_TAB` (default Explore).

Its reason for existing is testing: a tab that is never current is constructed but never laid
out, so any layout warning it would emit never appears. Naming a tab is what makes an
automated run render it. An unknown name is reported rather than silently ignored, since a
typo would otherwise look like the setting having no effect.
"""
function _initial_tab()
    want = lowercase(strip(get(ENV, "OITOOLSGUI_TAB", "")))
    isempty(want) && return 0
    i = something(findfirst(==(want), TAB_NAMES), findfirst(==(want), TAB_ALIASES), nothing)
    if i === nothing
        @warn "OITOOLSGUI_TAB is not a tab name; opening on Exploring" got = want options = TAB_NAMES
        return 0
    end
    return i - 1        # QML indexes from zero
end

function OITOOLS.gui(session::Session = Session();
                        qmlfile::AbstractString = joinpath(pkgdir(OITOOLS), "src", "gui", "qml", "Main.qml"),
                        autoquit_ms::Integer = 0,
                        on_ready = nothing)
    check_qt_conflict()
    # Backstop for sessions started by hand rather than through bin/oitoolsgui.jl. GLMakie is
    # imported by now but has not built a context yet, and Mesa reads the variables at context
    # creation, so this is still early enough. A no-op when the launcher already ran it.
    READY_HOOK[] = on_ready
    gfx = configure_graphics!()

    # Backstop for a session started by hand rather than through bin/oitoolsgui.jl, where
    # `prefer_native_wayland!` may never have been called. GLMakie is loaded by now, so GLFW's
    # platform is settled and cannot be changed -- but Qt's can, and the two must not differ.
    # Asking GLFW what it actually got is the whole check: on a Wayland session that answers
    # X11, Qt has to follow it to xcb or the process ends up with two windowing systems.
    #
    # A no-op through the launcher, which has already matched them, and whenever
    # QT_QPA_PLATFORM is set.
    on_x11 = GLMakie.GLFW.GetPlatform() == GLMakie.GLFW.PLATFORM_X11
    qt = configure_qt_platform!(; match_x11 = on_x11, verbose = false)

    uiscale = ui_scale_override()
    @info "UI scale: $(uiscale.reason)"

    # Fill the glyph atlas before ANY figure exists. Text first rasterised after the Qt window
    # is up renders corrupted; the same text set beforehand is clean. See `PLOT_GLYPHS` for the
    # measurements, and `test/gui/glyphtest.jl` for the harness that isolated it.
    #
    # A plain statement, deliberately: inside a `@debug` string this would never run, because
    # Julia's logging macros skip their interpolations when the level is off.
    nglyphs = prewarm_glyphs!()

    fig = Makie.Figure(fonts = PLOT_FONTS)
    ax  = Makie.Axis(fig[1, 1]; xlabel = "", ylabel = "", title = "OITOOLS")
    style_axis!(ax)
    # Build every plot BEFORE the window exists. Once QML is up, GLMakie allocates GPU buffers
    # on insertion and the context belongs to Qt's render thread -- see src/livecanvas.jl.
    canvas = build_canvas(fig, ax)

    # A second figure for the Image perspective. Sharing one canvas between Exploring and
    # Imaging would mean a reconstruction wiping the plot a user is reading, and a tab whose
    # own result appears somewhere else. Both figures are built here, before the window: the
    # GL context belongs to Qt's render thread afterwards, and inserting a plot then is what
    # allocates buffers with no context bound.
    imfig = Makie.Figure(fonts = PLOT_FONTS)
    imax  = Makie.Axis(imfig[1, 1]; xlabel = "α (mas)", ylabel = "δ (mas)",
                       title = "reconstruction")
    style_axis!(imax)
    imcanvas = build_canvas(imfig, imax)

    # And a third for the Gantt. Three figures rather than one shared: a reconstruction must
    # not wipe the plot being read on another tab, and a night plan is not a scatter — it needs
    # its own axis limits, its own y ticks and its own set of plots.
    gfig = Makie.Figure(fonts = PLOT_FONTS)
    gax  = Makie.Axis(gfig[1, 1]; xlabel = "LST (h)", title = "")
    style_axis!(gax)
    gantt = build_gantt(gfig, gax)

    dfig = Makie.Figure(fonts = PLOT_FONTS)
    dax  = Makie.Axis(dfig[1, 1])
    style_axis!(dax)
    delayplot = build_delay_plot(dfig, dax)

    # And one for the model rendering. `build_canvas` gives an axis with a heatmap already on
    # it, which is all this needs; `show_image!` then draws into it.
    mfig = Makie.Figure(fonts = PLOT_FONTS)
    max_ = Makie.Axis(mfig[1, 1]; title = "model")
    style_axis!(max_)
    modelcanvas = build_canvas(mfig, max_)

    # And a sixth for the χ² map. Its own figure for the same reason as the Gantt: it is not a
    # scatter and not an image -- it needs a colorbar of its own, contour levels of its own, and
    # axes labelled with parameter names rather than coordinates.
    ofig = Makie.Figure(fonts = PLOT_FONTS)
    oax  = Makie.Axis(ofig[1, 1])
    style_axis!(oax)
    chi2map = build_chi2_map(ofig, oax)

    # And a seventh for the radial-profile preview. Two axes side by side rather than two
    # figures: I(r) and V(B) are read together -- the question is always what a change to the
    # profile does to the signature the data actually measures.
    # Stacked, not side by side: this figure lives in the ~300dp preview column beside the
    # profile editor, and two axes across that width leaves the r ticks an unreadable smear.
    # And an eighth for the model residuals. It shares the χ² map's rectangle -- only one of
    # the two is ever visible -- but not its figure: that one is a heatmap with a colorbar in
    # `fig[1, 2]`, and this is a scatter with a legend there.
    rfig = Makie.Figure(fonts = PLOT_FONTS)
    rax  = Makie.Axis(rfig[1, 1])
    style_axis!(rax)
    residplot = build_residuals(rfig, rax)

    # And a ninth for the SED, in the same rectangle again. Its own figure rather than the
    # residuals': that one is a scatter against baseline, this is a set of curves against
    # wavelength, and an axis cannot be both.
    sfig = Makie.Figure(fonts = PLOT_FONTS)
    sax  = Makie.Axis(sfig[1, 1])
    style_axis!(sax)
    sedplot = build_sed(sfig, sax)

    prfig = Makie.Figure(fonts = PLOT_FONTS)
    priax = Makie.Axis(prfig[1, 1])
    prvax = Makie.Axis(prfig[2, 1])
    style_axis!(priax); style_axis!(prvax)
    profileplot = build_profile_plot(prfig, priax, prvax)

    sh  = ShellState(session, fig, ax, nothing, String[], Any[], 0, :uv, :baseline, false, false,
                     "no dataset loaded", String[], canvas, "", nothing, imcanvas, nothing, gantt,
                     delayplot, modelcanvas, nothing, Any[], chi2map, profileplot, residplot, sedplot,
                     "", "", nothing)
    SHELL[] = sh
    install_interactions!(sh)

    if !isempty(session.datasets)
        sh.current = 1
        sh.status  = refresh_plot!(sh)
    end
    # Seed the console before the window exists, so it opens with a transcript rather than an
    # empty box -- including the graphics and scale decisions, which are exactly what you want
    # to see when the window looks wrong.
    console!(sh, "OITOOLSGUI ready")
    console!(sh, "graphics: $(gfx.reason)")
    console!(sh, "qt platform: $(qt.reason)")
    console!(sh, "glyph atlas: $(nglyphs) glyphs pre-warmed")
    console!(sh, "ui scale: $(uiscale.reason)")
    for (i, d) in enumerate(session.datasets)
        console!(sh, "load_dataset!(session, \"$(d.path)\")"; kind = :cmd)
        i == sh.current && console!(sh, sh.status)
    end
    isempty(session.datasets) && console!(sh, "no dataset loaded — use Open OIFITS…")

    # `fig`, not `fig.scene`: QMLMakie's own examples pass the Figure, and MakieArea resolves
    # the scene itself (Makie.get_scene). Handing it the root Scene directly skips whatever
    # the figure-level path sets up.
    QML.loadqml(qmlfile;
                plot            = fig,
                imagePlot       = imfig,
                ganttPlot       = gfig,
                delayPlot       = dfig,
                modelPlot       = mfig,
                chi2MapPlot     = ofig,
                profilePlot     = prfig,
                residPlot       = rfig,
                sedPlot         = sfig,
                autoQuitMs      = Int(autoquit_ms),
                initialTab      = _initial_tab(),
                uiScaleOverride = uiscale.scale,
                initialFolder   = _initial_folder(session),
                fullscreenOnStart = get(ENV, "OITOOLSGUI_FULLSCREEN", "") == "1",
                initialStatus   = sh.status)
    QML.exec()

    # An automated run has no way to read the console pane once the window is gone, so hand it
    # over on the way out. Only active when asked for; production runs never write a file.
    dump = get(ENV, "OITOOLSGUI_CONSOLE_DUMP", "")
    if !isempty(dump)
        try
            open(dump, "w") do io
                for l in sh.console; println(io, l); end
                # Not console output: the last identified point, so an automated run can
                # check that picking works without reading a screenshot.
                println(io, "PICKTEXT: ", sh.picktext)
            end
        catch err
            @warn "could not write console dump" path = dump exception = err
        end
    end
    return session
end
