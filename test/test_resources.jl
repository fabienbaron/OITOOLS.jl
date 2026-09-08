# Shipped resources: where QML, the instrument configs and the demo data are found.
#
# The failure this guards against is invisible in a checkout, which is why it is tested rather
# than trusted: `pkgdir` resolves on the machine that built an application and nowhere else, so
# a relocated OITOOLS finds none of its data files. A missing config directory raises nothing
# at all — the facility list is simply empty — so only an assertion catches it.

@testset "shipped resources" begin
    root = pkgdir(OITOOLS)

    @testset "a checkout resolves everything from the package directory" begin
        @test resource_dir() == root
        for sub in (("src", "configs"), ("src", "gui", "qml", "Main.qml"),
                    ("demos", "data"), ("demos", "models"))
            @test OITOOLS.resource(sub...) == joinpath(root, sub...)
        end
        @test OITOOLS.resource("src", "no-such-resource") === nothing
    end

    @testset "the config directory is a function, not a baked constant" begin
        # A `const` would be evaluated at precompile time and carry the building machine's
        # path into the sysimage, which is exactly what an application cannot correct.
        @test isdir(OITOOLS._configs_dir())
        @test OITOOLS._configs_dir() == joinpath(root, "src", "configs")
        @test !isempty(list_configs())
    end

    @testset "a forced root wins, and is fallen through when it lacks the file" begin
        mktempdir() do dir
            mkpath(joinpath(dir, "src", "configs"))
            write(joinpath(dir, "src", "configs", "MADEUP.toml"), "name = \"MADEUP\"\n")
            withenv(OITOOLS.RESOURCE_DIR_VAR => dir) do
                @test resource_dir() == dir
                # taken from the forced root
                @test OITOOLS._configs_dir() == joinpath(dir, "src", "configs")
                @test OITOOLS.resource("src", "configs", "MADEUP.toml") ==
                      joinpath(dir, "src", "configs", "MADEUP.toml")
                # absent there, so the checkout still supplies it rather than being hidden
                @test OITOOLS.resource("demos", "models") == joinpath(root, "demos", "models")
                @test OITOOLS.resource("src", "gui", "qml", "Main.qml") ==
                      joinpath(root, "src", "gui", "qml", "Main.qml")
            end
        end
        # and the variable does not leak
        @test resource_dir() == root
    end

    @testset "the roots are ordered bundle before checkout" begin
        roots = OITOOLS._resource_roots()
        @test roots[end] == root
        @test any(r -> endswith(r, joinpath("share", "oitools")), roots)
        mktempdir() do dir
            withenv(OITOOLS.RESOURCE_DIR_VAR => dir) do
                @test first(OITOOLS._resource_roots()) == dir
            end
        end
    end
end
