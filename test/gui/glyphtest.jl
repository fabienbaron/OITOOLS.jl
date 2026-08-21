# Glyph-corruption diagnostic.
#
# Reported on three machines with different graphics stacks: characters in titles, axis labels
# and legends come out with strokes of a neighbouring atlas cell drawn through them, while tick
# labels stay clean. On one report the legend's scatter MARKERS are misplaced too, which puts
# the fault in Makie's texture atlas rather than in text shaping.
#
# Ruled out already, each by measurement:
#
#   * the font, the string and Makie      -- the identical figure saved to PNG by plain GLMakie
#                                            is clean (run this with BACKEND=png to confirm)
#   * Mesa's d3d12 driver                 -- GALLIUM_DRIVER=llvmpipe corrupts as well
#   * WSL                                 -- a plain Linux machine shows it too
#   * inserting glyphs at run time        -- pre-filling the atlas before the window did not fix it
#
# So this script exists to vary ONE thing at a time and see which one matters.
#
# Usage
# -----
#     julia --project=bin test/gui/glyphtest.jl
#
# Knobs, all environment variables, all with a default that reproduces what the GUI does:
#
# | variable | values | what it tests |
# |---|---|---|
# | `BACKEND` | `qml` (default), `png` | Qt's shared GL context, or plain GLMakie to a file. The control. |
# | `FONT`    | `dejavu` (default), `makie` | Makie pre-renders `a-z A-Z 0-9 . -` into its cached atlas ONLY for its own default fonts. Pinning DejaVu Sans means every glyph is inserted at run time. |
# | `FIGS`    | `1` (default) .. `6` | How many Makie figures share one Qt GL context. The GUI builds six. |
# | `WHEN`    | `before` (default), `after` | Whether the labels are set before the window exists or changed after it, which is when the GUI sets most of them. |
# | `PREWARM` | `0` (default), `1` | Fill the atlas with the whole character set before any figure is built. |
# | `SETUP`   | `1` (default), `0` | Run the launcher's `configure_graphics!` + `prefer_native_wayland!` before GLMakie loads. On Wayland, `0` leaves GLFW on XWayland while Qt is native. |
# | `QUITMS`  | milliseconds | Close automatically, so a harness can screenshot the window. |
#
# Read the window like this:
#
#   * the TICK LABELS are the control — they are digits, and they have been clean everywhere
#   * the title, the axis labels and the legend text are where corruption has been seen
#   * the marker row is there because markers come from the same atlas; if they are misplaced
#     or wrong, the fault is the atlas and not the text path
#   * the ASCII block prints every printable character, so the corrupted ones can be named
#     rather than described

using Printf

const BACKEND = lowercase(get(ENV, "BACKEND", "qml"))
const FONTSEL = lowercase(get(ENV, "FONT",    "dejavu"))
const NFIGS   = clamp(something(tryparse(Int, get(ENV, "FIGS", "1")), 1), 1, 6)
const WHEN    = lowercase(get(ENV, "WHEN",    "before"))
const PREWARM = get(ENV, "PREWARM", "0") in ("1", "true", "yes")
const QUITMS  = something(tryparse(Int, get(ENV, "QUITMS", "0")), 0)
const SETUP   = get(ENV, "SETUP", "1") in ("1", "true", "yes")

# The launcher's load order, reproduced exactly, because it is part of what is under test.
#
# `bin/oitoolsgui.jl` calls `configure_graphics!()` and then `prefer_native_wayland!()` BEFORE
# `using GLMakie`, since GLFW reads its windowing choice when the first GL context is created.
# Skip it and, on a Wayland session, GLFW takes XWayland while Qt takes Wayland natively --
# two windowing systems in one process sharing one GL context. `SETUP=0` is how to find out
# whether that pairing is what corrupts the atlas.
if SETUP
    using OITOOLS
    configure_graphics!()
    using GLFW_jll
    prefer_native_wayland!()
end

using GLMakie
# QML at top level, NOT inside the backend branch below: `QML.@qmlfunction` is expanded when
# that branch is lowered, which happens before a `using` inside the same branch has run.
using QML, QMLMakie
const Makie = GLMakie.Makie

# The GUI pins one family for every slot. That is what invalidates Makie's pre-rendered atlas,
# so it is the first thing worth toggling.
const PINNED = (regular = "DejaVu Sans", bold = "DejaVu Sans",
                italic = "DejaVu Sans", bold_italic = "DejaVu Sans")

# ── the strings under test ───────────────────────────────────────────────────
#
# Taken from the GUI verbatim where possible, so a corruption seen here is the same corruption.
const TITLE  = "v2 — AZCYG2011_FINAL2018.oifits"      # em dash, underscore, mixed case
const XLABEL = "Baseline (Mλ)"                         # the label reported as "(⁄⁄⁄)"
const YLABEL = "V²"                                    # superscript
const ASCII1 = "!\"#\$%&'()*+,-./0123456789:;<=>?"
const ASCII2 = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_"
const ASCII3 = "`abcdefghijklmnopqrstuvwxyz{|}~"
const GREEK  = "λ μ α δ χ Δ θ σ π  ° ² ³ ± × · — … ≤ ≥ ≈"
const LEGEND = ["E1-E2", "S1-E1", "S2-E1", "S2-W2", "W1-S1"]

const LEGCOLORS = [Makie.RGBAf(0, 0, 0, 1), Makie.RGBAf(0, 0, 1, 1),
                   Makie.RGBAf(0.8, 0.8, 0.8, 1), Makie.RGBAf(0, 1, 1, 1),
                   Makie.RGBAf(0.4, 0.4, 0.4, 1)]

"""
    prewarm!(fonts)

Insert every character this script draws into Makie's texture atlas up front.

Makie's `render_default_glyphs!` covers only `a-z`, `A-Z`, `0-9`, `.`, `-` and the Unicode
minus, and only for the fonts in its own default theme — so with a pinned family nothing is
pre-rendered and every glyph is added on first draw. This makes that difference testable.
"""
function prewarm!(fonts)
    atlas = Makie.get_texture_atlas()
    chars = Set{Char}()
    for s in (TITLE, XLABEL, YLABEL, ASCII1, ASCII2, ASCII3, GREEK, LEGEND...)
        union!(chars, collect(s))
    end
    n = 0
    for name in unique(values(fonts))
        font = try Makie.to_font(name) catch; continue end
        for c in chars
            try; Makie.insert_glyph!(atlas, c, font); n += 1; catch; end
        end
    end
    return n
end

"""
    build_figure(fonts; labelled) -> (figure, apply!)

The figure under test, and a closure that puts the labels on it.

Splitting the two is the point: `apply!` can be called before the window exists or after it,
which is the `WHEN` knob. The GUI sets almost every label after, from `draw!`.
"""
function build_figure(fonts; labelled::Bool)
    fig = Makie.Figure(; fonts, size = (1100, 760))
    ax  = Makie.Axis(fig[1, 1])

    # Tick labels are the control: digits only, and clean on every machine so far. Fixed
    # values so the same numbers appear in every run and two screenshots can be compared.
    ax.xticks = (collect(0.0:25.0:200.0), string.(0:25:200))
    ax.yticks = (collect(0.0:0.2:1.0), [@sprintf("%.1f", v) for v in 0.0:0.2:1.0])
    Makie.xlims!(ax, 0, 200); Makie.ylims!(ax, -0.05, 1.05)

    # Markers, from the same atlas as the glyphs. If these are misplaced or wrong while the
    # text is merely ugly, the atlas is implicated and the text path is not.
    for (i, c) in enumerate(LEGCOLORS)
        Makie.scatter!(ax, fill(20.0 * i, 5), collect(0.75:0.05:0.95);
                       color = c, markersize = 12)
    end
    Makie.scatter!(ax, collect(10.0:10.0:190.0), fill(0.62, 19);
                   color = :black, marker = :cross, markersize = 14)
    Makie.scatter!(ax, collect(10.0:10.0:190.0), fill(0.56, 19);
                   color = :black, marker = :rect, markersize = 10)

    # A hand-drawn legend, exactly as src/gui/livecanvas.jl draws one: marker plus text in the
    # same axis. This is the construct the corruption was reported on.
    legx = collect(20.0:35.0:160.0)
    Makie.scatter!(ax, legx, fill(0.40, length(legx));
                   color = LEGCOLORS, marker = :circle, markersize = 12)
    legtext = Makie.Observable(fill("", length(legx)))
    Makie.text!(ax, [Makie.Point2f(x, 0.40) for x in legx]; text = legtext,
                align = (:left, :center), offset = (12, 0), fontsize = 14)

    # The character inventory, so a corrupted glyph can be NAMED.
    rows = Makie.Observable(fill("", 4))
    Makie.text!(ax, [Makie.Point2f(5, y) for y in (0.28, 0.20, 0.12, 0.04)];
                text = rows, align = (:left, :center), fontsize = 15)

    apply! = function ()
        ax.title  = TITLE
        ax.xlabel = XLABEL
        ax.ylabel = YLABEL
        legtext[] = LEGEND
        rows[]    = [ASCII1, ASCII2, ASCII3, GREEK]
        return nothing
    end
    labelled && apply!()
    return fig, apply!
end

# ── what was actually run, printed so a screenshot can be matched to a config ──
function report()
    println("─"^72)
    @printf("BACKEND=%s  FONT=%s  FIGS=%d  WHEN=%s  PREWARM=%s  SETUP=%s\n",
            BACKEND, FONTSEL, NFIGS, WHEN, PREWARM, SETUP)
    for v in ("GALLIUM_DRIVER", "LIBGL_ALWAYS_SOFTWARE", "MESA_LOADER_DRIVER_OVERRIDE",
              "__GLX_VENDOR_LIBRARY_NAME", "DISPLAY", "WAYLAND_DISPLAY")
        haskey(ENV, v) && @printf("  %-28s %s\n", v, ENV[v])
    end
    # Which glyphs Makie pre-rendered, and for which fonts — the asymmetry this script exists
    # to test. `MAKIE_DEFAULT_THEME.fonts` is what `render_default_glyphs!` uses.
    try
        deffonts = unique(string.(Makie.to_value.(values(Makie.MAKIE_DEFAULT_THEME.fonts))))
        println("  Makie default fonts          ", join(deffonts, ", "))
    catch
    end
    println("  pinned fonts                 ", FONTSEL == "dejavu" ? PINNED.regular : "(Makie default)")
    println("─"^72)
end

report()
fonts = FONTSEL == "makie" ? Makie.MAKIE_DEFAULT_THEME.fonts : PINNED
PREWARM && println("pre-warmed ", prewarm!(fonts), " glyph insertions")

if BACKEND == "png"
    # The control. No Qt, no shared context: plain GLMakie straight to a file.
    fig, apply! = build_figure(fonts; labelled = WHEN == "before")
    WHEN == "before" || apply!()
    out = get(ENV, "OUT", joinpath(@__DIR__, "glyphtest.png"))
    Makie.save(out, fig)
    println("wrote ", out)
else
    # Figure 1 carries the test content; 2..N are decoys that merely EXIST, as the GUI's other
    # perspectives do. If corruption appears only once several figures share the context, that
    # is the answer.
    figs = Any[]
    applies = Any[]
    for i in 1:NFIGS
        f, a = build_figure(fonts; labelled = i == 1 && WHEN == "before")
        push!(figs, f); push!(applies, a)
    end
    while length(figs) < 6                      # QML always binds six properties
        push!(figs, Makie.Figure(; fonts))
    end

    # `after` is the case that matters: the GUI sets its labels once the window is up, from
    # `draw!`, and by then the GL context belongs to Qt's render thread.
    #
    # Driven by a QML Timer, NOT by `@async`. `QML.exec()` blocks inside C++ without yielding
    # to Julia's scheduler, so a task spawned before it never runs — the first version of this
    # script did exactly that and quietly tested nothing.
    function apply_labels()
        applies[1]()
        println("labels applied AFTER the window was shown")
        return ""
    end
    QML.@qmlfunction apply_labels

    qml = joinpath(@__DIR__, "glyphtest.qml")
    QML.loadqml(qml; plot1 = figs[1], plot2 = figs[2], plot3 = figs[3],
                     plot4 = figs[4], plot5 = figs[5], plot6 = figs[6],
                     nfigs = NFIGS, autoQuitMs = QUITMS,
                     applyAfter = (WHEN == "after"),
                     caption = "BACKEND=qml FONT=$FONTSEL FIGS=$NFIGS " *
                               "WHEN=$WHEN PREWARM=$PREWARM SETUP=$SETUP")
    QML.exec()
end
