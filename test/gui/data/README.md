# Test data

Two OIFITS files, copied from OITOOLS, so the GUI has something to open without depending on
a sibling clone being present or on a particular machine's paths. They are what
`test/gui_xvfb.jl` points the file dialog at, and what the manual checklist uses.

| file | what it exercises |
|------|-------------------|
| `2004-data1.oifits` | monochromatic, single epoch — the "wav" and "mjd" colour modes have nothing to spread over, which is the degenerate case that used to throw |
| `v1295Aql.oifits`   | polychromatic, several epochs — real colour-by-wavelength and colour-by-MJD, and enough baselines that the legend wraps |

Kept small deliberately: the whole point is that a GUI smoke test finishes in seconds.
