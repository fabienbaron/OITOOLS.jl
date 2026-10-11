# Levenberg-Marquardt (`_lm_bounded`, src/fit_lm.jl) against LsqFit's, which it transcribes.
# Where no bound is involved and the χ² test is off, the two must reach the same point; where a
# parameter is held on a bound they must not, and ours must reach the bounded minimum.
using OITOOLS, Test, LinearAlgebra

@testset "Levenberg-Marquardt with bounds that hold" begin
    O = OITOOLS
    LF = O.LsqFit
    # Three Gaussians on a background, 10 parameters, residuals in units of the noise.
    t64 = collect(range(-5, 5; length = 400))
    gm(t, p) = p[10] .+ p[1] .* exp.(-(t .- p[2]).^2 ./ (2p[3]^2)) .+
                        p[4] .* exp.(-(t .- p[5]).^2 ./ (2p[6]^2)) .+
                        p[7] .* exp.(-(t .- p[8]).^2 ./ (2p[9]^2))
    function gj(t, p)
        J = zeros(eltype(p), length(t), 10)
        for (a, c, s) in ((1, 2, 3), (4, 5, 6), (7, 8, 9))
            e = exp.(-(t .- p[c]).^2 ./ (2p[s]^2))
            J[:, a] = e
            J[:, c] = p[a] .* e .* (t .- p[c]) ./ p[s]^2
            J[:, s] = p[a] .* e .* (t .- p[c]).^2 ./ p[s]^3
        end
        J[:, 10] .= 1
        return J
    end
    ptrue = [1.0, -2.0, 0.6, 0.7, 0.5, 1.0, 0.4, 2.5, 0.3, 0.05]
    y64 = gm(t64, ptrue) .+ 0.01 .* sin.(37 .* t64)
    p0 = [0.8, -1.7, 0.8, 0.9, 0.2, 0.8, 0.5, 2.2, 0.4, 0.1]

    @testset "$T" for T in (Float64, Float32)
        t, y, σ = T.(t64), T.(y64), T(0.01)
        f(p) = (gm(t, p) .- y) ./ σ
        j(p) = gj(t, p) ./ σ
        x0 = T.(p0)
        free_lo, free_hi = fill(T(-Inf), 10), fill(T(Inf), 10)
        lsq(lo, hi) = LF.curve_fit((tt, p) -> gm(tt, p) ./ σ, (tt, p) -> gj(tt, p) ./ σ, t, y ./ σ,
                                   x0; lower = lo, upper = hi, maxIter = 1000, store_trace = true)

        # No bound, χ² test off: LsqFit's steps. Its trace has an iteration-0 entry. In Float32
        # its step test (1e-8, below Float32's resolution) cannot fire and it iterates on at the
        # same point; ours stops at eps(Float32).
        a = O._lm_bounded(f, j, x0, free_lo, free_hi; maxIter = 1000, chi2_tol = 0)
        b = lsq(free_lo, free_hi)
        @test a.converged && isempty(a.held)
        @test T == Float64 ? a.iterations == length(b.trace) - 1 :
                             a.iterations < length(b.trace) - 1
        @test a.x ≈ b.param rtol = (T == Float64 ? 1e-12 : 1e-5)
        @test eltype(a.x) == T

        # The χ² test stops earlier, within its tolerance of the same minimum.
        c = O._lm_bounded(f, j, x0, free_lo, free_hi)
        @test c.reason === :chi2 && c.iterations < a.iterations
        @test sum(abs2, c.resid) - sum(abs2, a.resid) < 1e-3 * sum(abs2, a.resid) / length(t)

        # The background wants 0.05; held at 0.08 it is on its bound, and the minimum is the one
        # with it fixed there. LsqFit, which clips instead of holding, does not reach it.
        lo = copy(free_lo); lo[10] = T(0.08)
        h = O._lm_bounded(f, j, x0, lo, free_hi)
        @test h.converged && h.held == [10] && h.x[10] == T(0.08)
        fixed = O._lm_bounded(p -> f(vcat(p, T(0.08))), p -> j(vcat(p, T(0.08)))[:, 1:9],
                              x0[1:9], free_lo[1:9], free_hi[1:9])
        @test sum(abs2, h.resid) ≈ sum(abs2, fixed.resid) rtol = 1e-4
        @test sum(abs2, h.resid) <= sum(abs2, lsq(lo, free_hi).resid) * (1 + 1e-6)
        @test h.iterations < 20

        # The cap still stops a fit, and says so.
        m = O._lm_bounded(f, j, x0, free_lo, free_hi; maxIter = 2)
        @test m.reason === :maxiter && !m.converged && m.iterations == 2
        @test_throws ArgumentError O._lm_bounded(f, j, x0, lo .+ T(1), free_hi)
    end
end
