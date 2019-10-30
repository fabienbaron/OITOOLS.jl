using Documenter, OITOOLS

DEPLOYDOCS = (get(ENV, "CI", nothing) == "true")

makedocs(
    modules=[OITOOLS],
    sitename = "OITOOLS",
    format = Documenter.HTML(
        prettyurls = DEPLOYDOCS,
    ),
    authors = "Fabien Baron and contributors",
    pages = ["index.md", "install.md", "introduction.md", "examples.md"]
)

if DEPLOYDOCS
    deploydocs(
        repo = "github.com/fabienbaron/OITOOLS.jl",
    )
end
