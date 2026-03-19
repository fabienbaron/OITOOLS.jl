using Pkg
Pkg.activate(@__DIR__)
CI = get(ENV, "CI", nothing) == "true"
ENV["PYTHON"] = ""
Pkg.build("PyCall")
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
        "Examples" => [
            "examples/intro.md",
            "examples/reading.md",
            "examples/plotting.md",
            "examples/modeling.md",
            "examples/simulating.md",
            "examples/imaging.md",
        ],
        "API Reference" => [
            "OIFITS Handling" => "api/oifits.md",
            "Plotting"        => "api/plotting.md",
            "Model Fitting"   => "api/modeling.md",
            "Imaging"         => "api/imaging.md",
            "Observation Planning" => "api/planning.md",
        ],
    ],
)

if CI
    deploydocs(;
        repo   = "github.com/fabienbaron/OITOOLS.jl",
        target = "build",
    )
end
