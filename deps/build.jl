using Pkg

if lowercase(get(ENV, "CI", "false")) == "true"
    #ENV["PYTHON"]=""; using Pkg; Pkg.add("PyCall"); Pkg.build("PyCall"); Pkg.add("PyPlot"); Pkg.add("Conda"); using Conda; Conda.add_channel("anaconda"); Conda.update(); Conda.add("numpy"); Conda.add("scipy"); Conda.add("matplotlib");
    let basepython = get(ENV, "PYTHON", "python2")
        envpath = joinpath(@__DIR__, "env")
        run(`pip install --user virtualenv`)
        run(`virtualenv --python=$basepython $envpath`)
        if Sys.iswindows()
            python = joinpath(envpath, "Scripts", "python.exe")
        else
            python = joinpath(envpath, "bin", "python2")
        end
        run(`$python -m pip install numpy`)
        run(`$python -m pip install scipy`)
        run(`$python -m pip install matplotlib`)

        ENV["PYTHON"] = python
        Pkg.build("PyCall")
    end
end
