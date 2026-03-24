using Pkg
Pkg.activate(@__DIR__)
CI = get(ENV, "CI", nothing) == "true"
using Documenter, OITOOLS

makedocs(;
    modules=[OITOOLS],
    sitename = "OITOOLS",
    checkdocs = :exports,
    doctest = false,
    format = Documenter.HTML(;
        prettyurls = CI,
        collapselevel = 2,
    ),
    authors = "Fabien Baron and contributors",
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
        "Demo Scripts" => "examples/intro.md",
        "Reading OIFITS" => "examples/reading.md",
        "Plotting" => "examples/plotting.md",
        "Model Fitting" => "examples/modeling.md",
        "Simulating" => "examples/simulating.md",
        "Imaging" => "examples/imaging.md",
        "BSDMM Reconstruction" => "examples/bsdmm.md",
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
