using Pkg
Pkg.activate(@__DIR__)
CI = get(ENV, "CI", nothing) == "true"
# Plotting lives in OITOOLSPythonPlotExt, not in OITOOLS itself. Loading PythonPlot is what
# brings that module into existence; without it every `@docs` block naming a plotting
# function fails with "no docs found", because the docstrings are attached to methods that
# have not been defined yet. Agg first: makedocs must not open a window, and the default
# backend probe is also what used to map a conflicting Qt into the process.
ENV["MPLBACKEND"] = get(ENV, "MPLBACKEND", "Agg")
using Documenter, OITOOLS, PythonPlot
# The Makie plotting API and the corner plot live in their own extensions, for the same reason
# the matplotlib one does: `@docs` blocks naming their functions find nothing unless the module
# that defines them exists.
using CairoMakie, PairPlots

const PLOTEXT = Base.get_extension(OITOOLS, :OITOOLSPythonPlotExt)
PLOTEXT === nothing && error("OITOOLSPythonPlotExt did not load; plotting docs would be lost")
const MAKIEEXT = Base.get_extension(OITOOLS, :OITOOLSMakieExt)
MAKIEEXT === nothing && error("OITOOLSMakieExt did not load; the Makie plotting docs would be lost")
const PAIREXT = Base.get_extension(OITOOLS, :OITOOLSPairPlotsExt)
PAIREXT === nothing && error("OITOOLSPairPlotsExt did not load; corner_figure would be undocumented")

makedocs(;
    modules=[OITOOLS, PLOTEXT, MAKIEEXT, PAIREXT],
    sitename = "OITOOLS",
    checkdocs = :exports,
    doctest = false,
    format = Documenter.HTML(;
        prettyurls = CI,
        collapselevel = 2,
        # api/modeling.md is a single page carrying the whole model-fitting API, docstrings
        # included, and passes Documenter's 100 KiB advisory. Raised rather than split: the
        # value of that page is that everything about fitting a model is on it. The hard
        # threshold, the one that fails the build, stays at its 200 KiB default.
        size_threshold_warn = 150 * 1024,
    ),
    authors = "Fabien Baron and contributors",
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
        "Guides" => [
            "Demo Scripts"           => "examples/intro.md",
            "Reading OIFITS"         => "examples/reading.md",
            "Plotting"               => "examples/plotting.md",
            "Model Fitting"          => "examples/modeling.md",
            "Uncertainties"          => "examples/uncertainties.md",
            "Simulating"             => "examples/simulating.md",
            "Imaging"                => "examples/imaging.md",
            "SPARCO"                 => "examples/sparco.md",
            "Hybrid (Model + Image)" => "examples/hybrid.md",
            "Maximum Entropy"        => "examples/maxent.md",
            "BSDMM Reconstruction"   => "examples/bsdmm.md",
            "SQUEEZE (MCMC)"         => "examples/squeeze.md",
        ],
        "API Reference" => [
            "OIFITS Handling"            => "api/oifits.md",
            "Plotting"                   => "api/plotting.md",
            "Model Fitting"              => "api/modeling.md",
            "Imaging (Setup)"            => "api/imaging_setup.md",
            "Imaging (Gradient Descent)" => "api/imaging_gd.md",
            "Imaging (BSDMM)"           => "api/imaging_bsdmm.md",
            "Imaging (Maximum Entropy)"  => "api/imaging_maxent.md",
            "Imaging (SPARCO)"           => "api/imaging_sparco.md",
            "Imaging (+ Modeling)"       => "api/imaging_hybrid.md",
            "Imaging (Variational Inference)" => "api/imaging_vi.md",
            "Imaging (MCMC)"             => "api/imaging_mcmc.md",
            "Observation Planning"       => "api/planning.md",
        ],
    ],
)

if CI
    deploydocs(;
        repo   = "github.com/fabienbaron/OITOOLS.jl",
        target = "build",
    )
end
