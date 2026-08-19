# Manual checks

Everything here needs a human. Everything that does **not** is a script, and lives elsewhere:

| | what it covers |
|---|---|
| `test/gui/runtests.jl` | 308 assertions, run under Xvfb |
| `test/gui/gui_xvfb.jl` | 11 scripted actions against a live window |
| `test/gui/gui_click.sh` | real clicks: file dialog, ComboBox popups, redraw after both |
| `test/gui/gui_click_wayland.sh` | the same, on GPU hardware and native Wayland |

Deliberately short: a checklist nobody finishes is worse than no checklist. Markdown rather
than `.jl` because it contains no code — `bin/oitoolsgui.jl` is the launcher and *is* code.

```
julia --project=bin bin/oitoolsgui.jl
```

The GUI lives in `OITOOLSGUIExt`, an extension of OITOOLS activated by GLMakie + QMLMakie +
QML. `using OITOOLS` alone loads neither Makie nor Qt.

## 1. Mouse interaction with the plot

Nothing automated touches these. They are handled inside Makie, below the level our callbacks
see, so a script cannot confirm them.

- [ ] **Scroll** over the plot zooms.
- [ ] **Drag** pans.
- [ ] **Left click** on a point fills the line under the plot with its identification —
      station pair, baseline in Mλ, wavelength, MJD, value and error. (Automated in
      `gui_click.sh`, which asserts a point is identified; check it reads sensibly.)
- [ ] **Right click** resets the view to the whole dataset.

`MakieArea` gates its input on having keyboard focus, and a `MouseArea` of ours over the
scene would eat the press that grants it. If any of these is dead, say *which* — "click does
nothing" and "zoom does nothing" point at different causes.

## 2. How it looks

- [ ] **Window size** on first open. Default 1280×1040 logical, scaled by `OITOOLSGUI_SCALE`
      (default 1.25) and clamped to 95% of the screen.
- [ ] **Plot fonts and markers.** Default 1.7× the publication sizes. Tune with
      `OITOOLSGUI_PLOT_SCALE=2.0`, or `1.0` to see exactly what a paper figure looks like.
- [ ] **Baseline legend** is compact and centred under the plot, not spanning the full width.
- [ ] **Colorbar** appears for `wav` and `mjd`, and is *completely absent* for `baseline` —
      no ghost 0..1 scale beside the plot.
- [ ] **uv coverage is isotropic.** Isotropic means equal *scale* on both axes — what
      `axis("equal")` gives in oiplot. The box is not square unless the data is.

## 3. Only on your machine

- [ ] **Which renderer you actually got:**
      ```
      QSG_INFO=1 julia --project=bin bin/oitoolsgui.jl 2>&1 | grep -i -A3 OpenGL
      ```
      Want `D3D12 (NVIDIA ...)` under WSL. `llvmpipe` means the GPU is not being used. The
      `libEGL` warnings are noise either way.

- [ ] **Wayland vs XWayland.** On a Wayland session (WSLg included) GLMakie is forced onto
      native Wayland, because GLFW.jl otherwise hard-codes X11 and you end up with Qt on
      Wayland and GLMakie on XWayland — two windowing systems in one process. To compare:
      ```
      OITOOLSGUI_GLFW_X11=1 julia --project=bin bin/oitoolsgui.jl
      ```
      If the GUI behaves better with that set, say so and it becomes the default on WSL.

- [ ] **Startup hangs with no window and no error?** Look for a stale CondaPkg lock. The
      launcher warns about it now, and it is the one failure that looks exactly like a
      graphics hang:
      ```
      rm -f .CondaPkg/lock
      ```
