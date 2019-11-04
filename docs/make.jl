using Documenter, OITOOLS

DEPLOYDOCS = (get(ENV, "CI", nothing) == "true")

makedocs(
    modules=[OITOOLS],
    sitename = "OITOOLS",
    format = Documenter.HTML(
        prettyurls = DEPLOYDOCS,
    ),
    authors = "Fabien Baron and contributors",
    pages = [ "Home" => "index.md",
    "Installation" => "install.md",
    "Examples" => Any[
                "examples/intro.md",
                "examples/reading.md",
                "examples/plotting.md",
                "examples/modeling.md",
                "examples/simulating.md",
                "examples/imaging.md"]
                ]
)

if DEPLOYDOCS
    deploydocs(
        repo = "github.com/fabienbaron/OITOOLS.jl",
    )
end
