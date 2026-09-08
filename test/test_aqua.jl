# Package hygiene. Cheap, and it catches a class of problem the other suites cannot see:
# an export that names nothing, a dependency with no version bound, method ambiguities.
#
# `test_deps_compat` is the one that earns its keep. Without a bound, resolution is free to
# pick a BREAKING release of any dependency into a user's install, and the first sign of it is
# their broken environment rather than anything CI would show.
#
# `test_piracies` matters for a second reason beyond correctness: type piracy invalidates
# precompiled code in other packages, which is recompilation paid at every load.

using Aqua

@testset "package hygiene (Aqua)" begin
    # Two of these are declared for reasons no static check can see, and both are load-bearing:
    #
    #   CondaPkg  Preferences resolves only for DIRECT dependencies, so [preferences.CondaPkg]
    #             in Project.toml does nothing unless CondaPkg is one. Reached indirectly
    #             through PythonCall it returns nothing, with no error.
    #   Pkg       reached lazily through `Base.require` in gui_launcher.jl's `_offer_to_install`,
    #             which keeps Pkg off the load path until someone answers the prompt.
    Aqua.test_all(OITOOLS; stale_deps = (ignore = [:CondaPkg, :Pkg],),
                           # The GUI stack is a weak dependency and is not loaded here, so the
                           # ambiguity walk would report on packages this test never brought in.
                           ambiguities = (broken = false,))
end
