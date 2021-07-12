using Pkg

if lowercase(get(ENV, "CI", "false")) == "true"
    let basepython = get(ENV, "PYTHON", "python3")
        envpath = joinpath(@__DIR__, "env")
        run(`pip3 install --user virtualenv`)
        run(`virtualenv --python=$basepython $envpath`)
        if Sys.iswindows()
            python = joinpath(envpath, "Scripts", "python.exe")
        else
            python = joinpath(envpath, "bin", "python3")
        end
        run(`$python -m pip install numpy`)
        run(`$python -m pip install scipy`)
        run(`$python -m pip install matplotlib`)
        run(`$python -m pip install ultranest`)
        run(`$python -m pip install astroquery`)
        ENV["PYTHON"] = python
        Pkg.build("PyCall")
    end
else
    ENV["PYTHON"]=""
    Pkg.add("PyCall")
    Pkg.build("PyCall")
    Pkg.add("Conda")
    using Conda
    Conda.add("astroquery", channel="astropy")
    Conda.add("ultranest", channel="conda-forge")
end
