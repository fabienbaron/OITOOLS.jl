using Pkg
Pkg.activate(@__DIR__)
CI = get(ENV, "CI", nothing) == "true"
ENV["PYTHON"] = ""
Pkg.build("PyCall")
using Documenter, OITOOLS

makedocs(;
    modules=[OITOOLS],
    sitename = "OITOOLS",
    doctest = false,
    format = Documenter.HTML(;
        prettyurls = CI,
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
    ],
)

if CI
    deploydocs(;
        repo   = "github.com/fabienbaron/OITOOLS.jl",
        target = "build",
    )
end
