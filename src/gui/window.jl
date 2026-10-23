# The window: Qt callback registration and `gui()`.
#
# Split from the rest of src/gui/ because everything else describes what to draw, while this
# builds the thing that draws it and hands control to Qt's event loop.

function __init__()
    # Callbacks must be registered before any QML file that calls them is loaded. Each one goes
    # through `_guarded`: QML.jl calls Julia from a C function pointer with no `try`, so an error
    # escaping a callback would unwind through Qt's event loop and freeze the window.
    for f in (shell_open, shell_set_view, shell_select_dataset,
                     shell_dataset_names, shell_current_dataset, shell_dataset_summary,
                     shell_filter_choices, shell_filter_preview, shell_apply_filter,
                     shell_filter_spec, shell_filter_preview_clear,
                     shell_set_select_tool, shell_selection_text, shell_clear_selection,
                     shell_selection_spec, shell_export_oifits, shell_script, shell_export, shell_pick,
                     shell_console, shell_ready, shell_about, shell_reset_zoom, shell_pick_text,
                     shell_observables, shell_plot_kinds, shell_grouping_noun, shell_gantt_hover, shell_reconstruct, shell_imaging_summary,
                     shell_image_defaults, shell_ft_setup, shell_chi2_breakdown,
                     shell_facilities, shell_telescopes, shell_telescope_order, shell_gantt, shell_best_pops, shell_pop_plans,
                     shell_shift_date, shell_sim_source, shell_simulate,
                     shell_model_image, shell_simbad, shell_ui_scale,
                     shell_fit_model, shell_fit_values, shell_fit_rows, shell_fit_correlation,
                     shell_corner_plot, shell_save_corner,
                     shell_image_colormaps, shell_image_colormap, shell_recenter_image,
                     shell_show_start_image, shell_vi_prior_sample,
                     shell_save_image, shell_engine_output,
                     shell_channels, shell_show_channel,
                     shell_tempering_diagnostics,
                     shell_job_poll, shell_job_stop, shell_job_running,
                     shell_reset_image_zoom,
                     shell_set_plot_scale, shell_set_marker_size,
                     shell_save_settings, shell_load_settings, shell_reset_settings,
                     shell_settings_path,
                     shell_controls_styles, shell_default_controls_style,
                     shell_ui_font_files, shell_plot_scale,
                     shell_version,
                     shell_set_zoom_step,
                     shell_model_rows, shell_model_components, shell_model_inspection,
                     shell_model_chi2, shell_set_param, shell_model_warnings, shell_fit_output,
                     shell_component_kinds, shell_add_component, shell_remove_component, shell_sparco_fluxes,
                     shell_rv_open, shell_rv_summary, shell_rv_close, shell_orbit_component, shell_orbit_view,
                     shell_add_image_component, shell_image_info, shell_set_image,
                     shell_remove_param, shell_add_param, shell_link_params, shell_param_name_problem,
                     shell_param_removable, shell_remove_unused,
                     shell_rename_component, shell_save_figure,
                     shell_expression_keywords, shell_expression_params, shell_check_expression,
                     shell_profile_params, shell_add_profile_params, shell_profile_grid,
                     shell_profile_curves,
                     shell_flux_constraint,
                     shell_default_bounds, shell_check_constraints,
                     shell_profile_templates, shell_profile_template,
                     shell_component_geometry, shell_set_component_geometry,
                     shell_ld_laws, shell_ld_law, shell_set_ld_law,
                     shell_kind_labels, shell_model_description,
                     shell_az_modes, shell_add_az_mode, shell_remove_az_mode,
                     shell_model_residuals, shell_model_sed, shell_model_depends,
                     shell_set_residual_mode, shell_set_overlay_mode,
                     shell_optional_engines, shell_compare_available,
                     shell_sparco_converged,
                     shell_result_ensemble, shell_show_result,
                     shell_bin_info, shell_rebin,
                     shell_chi2_map_info, shell_free_names,
                     shell_companion_candidates, shell_companion_grid, shell_companion_search,
                     shell_companion_map_view, shell_companion_progress, shell_grid_progress,
                     shell_multistart, shell_multistart_minima, shell_multistart_select,
                     shell_multistart_view, shell_multistart_names, shell_multistart_pmoired,
                     shell_lines, shell_add_line, shell_remove_line, shell_line_by_velocity,
                     shell_companion_defaults,
                     shell_smearing, shell_set_smearing, shell_nested_backend, shell_sim_band,
                     shell_open_model, shell_save_model,
                     shell_tools_source, shell_tools_headers, shell_tools_hdu, shell_tools_validate,
                     shell_tools_bin, shell_tools_targets, shell_tools_export,
                     shell_tools_split_preview, shell_tools_split, shell_tools_merge,
                     shell_cube_open, shell_cube_show, shell_cube_axes, shell_cube_show_at, shell_cube_tiles, shell_cube_power, shell_table_show, shell_table_sort, shell_table_filter, shell_table_channel, shell_table_comment, shell_table_widths, shell_tools_group,
                     shell_fits_keywords, shell_fits_apply, shell_fits_keyword_source, shell_fits_save,
                     shell_fit_each, shell_epochs_table, shell_save_epochs,
                     shell_import_pmoired, shell_import_pmoired_text, shell_export_pmoired,
                     # The in-window file picker: QtQuick.Dialogs leaves its window mapped on
                     # some systems, so the picker is drawn inside our own window and needs
                     # Julia for everything the toolkit would otherwise supply.
                     picker_list, picker_places, picker_join, picker_parent, picker_basename, picker_suggest,
                     shell_initial_folder, shell_set_work_directory,
                     picker_start, picker_kind, picker_would_overwrite, picker_examples,
                     shell_walk_request, shell_walk_reply)
        QML.qmlfunction(string(nameof(f)), (args...) -> _guarded(f, args...))
    end
    return nothing
end

"""
    _guarded(f, args...)

Call the Qt callback `f`, and turn an error into the panels' error string instead of letting it
reach Qt. The error goes to the console pane as `! name: cause` and, with its stack trace, to
stderr; the return value is that same `"! name: cause"`, which the panels already read as a
failure. A callback whose caller expects a number or a Bool gets the string instead, which
JavaScript coerces: one control shows a wrong value, the window keeps running.
"""
function _guarded(f, args...)
    try
        return f(args...)
    catch err
        msg = string("! ", nameof(f), ": ", _cause(err))
        # stderr first: the console pane needs a running shell, and the callback may have failed
        # because there is none.
        Base.showerror(stderr, err, catch_backtrace())
        println(stderr)
        try
            SHELL[] === nothing || console!(SHELL[], msg; kind = :err)
        catch
        end
        return msg
    end
end

"""
Folder the file dialog opens in.

Next to the last dataset loaded if there is one, then `\$OITOOLSGUI_DATA_DIR`, then the
package's `demos/data`. That last is the user-facing data directory — it holds Alpha Cen A and
B, the 2004 beauty-contest set, and `Beauty_Contest_2026/`, whose files are the only ones shipped with
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
        isdir(d) && return file_url(d)
    end
    forced = get(ENV, "OITOOLSGUI_DATA_DIR", "")
    isempty(forced) || (isdir(forced) && return file_url(abspath(forced)))
    # Same order as `picker_start`, and for the same reason: a saved work directory outranks the
    # shipped examples, which are only where a first session has nothing better to offer.
    w = gui_work_directory()
    isempty(w) || return file_url(w)
    for sub in (joinpath("demos", "data"), joinpath("test", "gui", "data"))
        p = OITOOLS.resource(sub)
        p === nothing || return file_url(p)
    end
    return file_url(pwd())
end

"""
    shell_initial_folder() -> String

Where a file dialog should open, worked out NOW rather than when the window was built.

QML had the startup value as a context property and passed it to every Open: a work directory
chosen mid-session could not reach it, so Open kept returning to the demo data until the
application was restarted while Save -- which asks for the folder afresh -- followed along.
"""
function shell_initial_folder()
    sh = SHELL[]
    sh === nothing && return picker_start("")
    return _initial_folder(sh.session)
end

"Tab names, in the order Main.qml lists them."
const TAB_NAMES = ("exploring", "observing", "modeling", "imaging", "tools")

"The bare noun for each tab, accepted as well: nobody types the gerund from memory."
const TAB_ALIASES = ("explore", "observe", "model", "image", "tool")

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

"""
    _main_qml() -> String

`Main.qml`, from wherever the resources are.

The one resource with no fallback: without it there is no window, so a missing file is reported
here, naming the roots that were searched. `loadqml`'s own error names only the last path tried,
which in a relocated application is a directory on the machine that built it.
"""
function _main_qml()
    p = OITOOLS.resource("src", "gui", "qml", "Main.qml")
    p === nothing && error("""
        Main.qml was not found. The GUI's QML files are shipped resources, and none of the
        resource roots holds them:

        $(join("    " .* OITOOLS._resource_roots(), "\n"))

        Set \$$(OITOOLS.RESOURCE_DIR_VAR) to the directory that contains src/gui/qml.
        """)
    return p
end

"Qt Quick Controls styles the bundled Qt carries, in the order the settings panel lists them."
const CONTROLS_STYLES = ("Basic", "Fusion", "Universal", "Material", "Imagine", "FluentWinUI3")

# What the window ships with. Fusion rather than Basic -- Basic is Qt's own default on Linux
# and draws its checkboxes and spin boxes noticeably larger than the rest of the window is
# scaled for. Named once here, since the settings panel offers a reset to it.
const DEFAULT_CONTROLS_STYLE = "Fusion"

"""
    shell_controls_styles() -> String

The style in force, then every style that can be chosen, one per line.

QML asks rather than carrying its own copy of the list: a second list would be free to drift
from `CONTROLS_STYLES`, and offering a style the bundled Qt does not have would fail silently
at the next launch rather than at the click.
"""
shell_controls_styles() =
    join((get(ENV, "QT_QUICK_CONTROLS_STYLE", DEFAULT_CONTROLS_STYLE), CONTROLS_STYLES...), '\n')

"""
    shell_default_controls_style() -> String

The style the window ships with, for the settings panel's reset.

Asked for rather than repeated in QML: the reset has to restore what an unconfigured install
would run, and a literal in the panel would go on claiming the old answer after this changed.
"""
shell_default_controls_style() = DEFAULT_CONTROLS_STYLE

"""
    shell_ui_font_files() -> String

The font files the window hands to Qt itself, one `file://` URL per line, regular then bold.

**Qt reads families from the SYSTEM font database, which is not where our fonts live.** A
family named in the settings panel resolves only if the machine happens to have it installed --
stock Windows has no Noto at all -- so an unconfigured window looks different on every
platform. These two faces ride in the MakieAssets artifact, which is already bundled for
`PLOT_FONT`, so a `FontLoader` on them costs nothing to ship and gives the same UI font
everywhere.

Empty when the assets are not where Makie says they are, in which case QML falls back to the
platform's own font rather than to a family name that would not resolve.
"""
function shell_ui_font_files()
    urls = String[]
    for f in ("NotoSans-Regular.ttf", "NotoSans-Bold.ttf")
        path = try
            Makie.assetpath("fonts", f)
        catch err
            @debug "Makie assets are not where they were expected" err file = f
            ""
        end
        isempty(path) || !isfile(path) || push!(urls, file_url(path))
    end
    return join(urls, '\n')
end

"""
    apply_controls_style!() -> String

Choose the Qt Quick Controls style, from the settings file if one was saved.

**The style is geometry as well as colour.** Basic is Qt's own default on Linux, and draws its
checkboxes and spin boxes larger than the rest of the window is scaled for; Fusion draws the
same controls appreciably smaller and more like a desktop toolkit, which is why it is
`DEFAULT_CONTROLS_STYLE`; Material and Universal are touch-sized and follow their own theming
rather than the palette; FluentWinUI3 is Windows 11's, which tracks the SYSTEM light/dark
setting and is what made the window come out half dark on a dark desktop.

Read from `gui_settings_file()` rather than from a constant so the settings panel can offer it.

`\$QT_QUICK_CONTROLS_STYLE` still wins: an environment that has already chosen is never
overridden, which keeps the escape hatch for a machine where a style misbehaves.

It can be set here at all — unlike the WINDOWING platform, which cannot — because the style is
read when the first Controls component is instantiated, at `loadqml`, not when
`QGuiApplication` is constructed. That is still ahead of us even in a compiled application.
Changing it therefore takes effect on the next launch, not the current one.
"""
# Qt on Windows reads the environment through msvcrt's copy (`_wgetenv_s`, `getenv` in
# Qt6Core.dll), which `ENV` (SetEnvironmentVariableW) does not update; msvcrt's `_wputenv`
# updates that copy and the process environment both.
function _setenv_for_qt!(name::AbstractString, value::AbstractString)
    ENV[name] = value
    Sys.iswindows() && ccall((:_wputenv, "msvcrt"), Cint, (Cwstring,), string(name, '=', value))
    return value
end

function apply_controls_style!()
    haskey(ENV, "QT_QUICK_CONTROLS_STYLE") && return ENV["QT_QUICK_CONTROLS_STYLE"]
    want = DEFAULT_CONTROLS_STYLE
    try
        path = gui_settings_file()
        if isfile(path)
            saved = get(TOML.parsefile(path), "controls_style", "")
            saved isa AbstractString && saved in CONTROLS_STYLES && (want = saved)
        end
    catch err
        @debug "could not read the saved controls style; using the default" err DEFAULT_CONTROLS_STYLE
    end
    return _setenv_for_qt!("QT_QUICK_CONTROLS_STYLE", want)
end

# When `exec()` returns, Qt may already have destroyed the window's GL context (it has on
# Windows), and GLMakie's atexit `closeall()` would then free GPU buffers through it and print a
# `ContextNotAvailable` trace. The window's screens are dropped from GLMakie's registries
# instead, with their texture-atlas entries; their GL objects went with the context.
function _forget_qml_screens!()
    isqml(s) = s isa GLMakie.Screen{QMLMakie.QMLWindow}
    try
        for s in filter(isqml, collect(GLMakie.ALL_SCREENS))
            s.glscreen.context.valid = false
            delete!(GLMakie.ALL_SCREENS, s)
            delete!(GLMakie.SCREEN_REUSE_POOL, s)
        end
        filter!(!isqml, GLMakie.SINGLETON_SCREEN)
        for key in collect(keys(GLMakie.atlas_texture_cache))
            key[2] isa QMLMakie.QMLWindow || continue
            tex, callback = pop!(GLMakie.atlas_texture_cache, key)
            tex.id = 0
            Makie.remove_font_render_callback!(key[1], callback)
        end
    catch err
        @debug "could not detach the window's screens from GLMakie" err
    end
    return nothing
end

function OITOOLS.gui(session::Session = Session();
                        qmlfile::AbstractString = _main_qml(),
                        autoquit_ms::Integer = 0,
                        on_ready = nothing)
    check_qt_conflict()
    # Backstop for sessions started by hand rather than through bin/oitoolsgui.jl. GLMakie is
    # imported by now but has not built a context yet, and Mesa reads the variables at context
    # creation, so this is still early enough. A no-op when the launcher already ran it.
    READY_HOOK[] = on_ready
    apply_controls_style!()
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

    # @debug, not @info: the usual answer is "unset, QML scales from the screen", which is a
    # line of terminal output every launch saying nothing happened. It is worth reading when a
    # window comes up the wrong size, so it stays reachable with
    # `JULIA_DEBUG=OITOOLSGUIExt`, and the settings panel shows the value regardless.
    uiscale = ui_scale_override()
    @debug "UI scale: $(uiscale.reason)"

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

    # And a seventh for the radial-profile preview. Two axes in one figure rather than two
    # figures: I(r) and V(B) are read together -- the question is always what a change to the
    # profile does to the signature the data actually measures. Side by side, because the
    # preview takes whatever width the pane has beside a fixed-width editor, and a wide frame
    # stacked into two strips leaves each plot a sliver of its own height.
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
    prvax = Makie.Axis(prfig[1, 2])
    style_axis!(priax); style_axis!(prvax)
    profileplot = build_profile_plot(prfig, priax, prvax)

    # And one for the Orbit view: the sky track beside the radial velocities.
    orbfig = Makie.Figure(fonts = PLOT_FONTS)
    orbsky = Makie.Axis(orbfig[1, 1]); orbrv = Makie.Axis(orbfig[1, 2])
    style_axis!(orbsky); style_axis!(orbrv)
    ORBIT_CANVAS[] = build_orbit_view(orbfig, orbsky, orbrv)

    # And one for the Tools cube viewer, drawn by `show_image!` as the reconstruction is.
    cfig = Makie.Figure(fonts = PLOT_FONTS)
    cax  = Makie.Axis(cfig[1, 1]; xlabel = "α (mas)", ylabel = "δ (mas)", title = "cube")
    style_axis!(cax)
    CUBE_CANVAS[] = build_canvas(cfig, cax)
    # And the cube viewer's channels side by side.
    tfig = Makie.Figure(fonts = PLOT_FONTS)
    CUBE_TILES[] = build_cube_tiles(tfig)
    # And the Tools tab's Data view of a table.
    oitable = build_oitable_model()

    sh  = ShellState(session, fig, ax, nothing, String[], Any[], 0, :uv, :baseline, false,
                     :none, :none, false,
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
        console!(sh, "load_dataset!(session, $(_literal(sources(d))))"; kind = :cmd)
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
                cubePlot        = cfig,
                cubeTilesPlot   = tfig,
                oitableModel    = oitable,
                orbitPlot       = orbfig,
                autoQuitMs      = Int(autoquit_ms),
                initialTab      = _initial_tab(),
                uiScaleOverride = uiscale.scale,
                initialFolder   = _initial_folder(session),
                fullscreenOnStart = get(ENV, "OITOOLSGUI_FULLSCREEN", "") == "1",
                walkDir         = get(ENV, "OITOOLSGUI_WALK", ""),
                # The startup scale and settings lines are diagnostics, not news: they print on
                # every launch and say nothing a working window has not already shown. Kept
                # behind a switch because the scale line is how the DPI exponent gets refitted
                # -- it reports the screen's real physical dpi, which is the measurement.
                verboseStartup  = get(ENV, "OITOOLSGUI_VERBOSE", "") == "1",
                initialStatus   = sh.status)
    QML.exec()
    _forget_qml_screens!()

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
