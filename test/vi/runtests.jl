# Variational inference: the suites that came across with OIVI.
#
# Gated on VarInf the way the tempering tests are gated on Pigeons, and skipped on CI
# (`CI=true`, which GitHub Actions sets): the limb-darkening profile's bagged runs take most of
# an hour on a 2-core runner. They are the real check on the shared adjoint (the FD / JVP / VJP
# suite fails on a wrong sign long before any image looks odd) and on VarInf's correlated field,
# so run them locally after changing either:
#
#     julia --project=bin test/vi/runtests.jl

using Test

const _HAVE_VARINF = try
    @eval using VarInf
    true
catch
    false
end

@testset "variational inference" begin

if !_HAVE_VARINF
    @info "VarInf not available — skipping the VI suites (unregistered weak dependency)"
    @test true
elseif get(ENV, "CI", "") == "true"
    @info "CI — skipping the VI suites; run test/vi/runtests.jl locally"
    @test true
else
    for f in ("test_core.jl", "test_precision.jl", "test_diffphase.jl", "test_pointsource.jl",
              "test_protocol_smoke.jl", "test_protocol_reconstruct.jl", "test_protocol_e2e.jl",
              "test_ldprofile.jl")
        # Each suite in its own module. They came from separate `julia test/x.jl` runs and
        # each defines its own `const oifitsfile`, `const data` and so on at top level;
        # included into one namespace the second definition is an error, not a shadow.
        @testset "$f" begin
            m = Module(Symbol("VI_", replace(f, "." => "_")))
            Core.eval(m, :(using Test))
            Base.include(m, joinpath(@__DIR__, f))
        end
    end
end

end
