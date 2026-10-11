# hankel_nodes.jl
#
# The quadrature a radial-profile component is integrated with: where its nodes go, how many
# there are, and what the profile itself says about where they must go.
#
# A profile is integrated in s = r / r_max ∈ [s_min, 1], μ = √(1 − s²), with nodes and weights
# such that ∫ f(s) s ds ≈ Σ gₖ f(sₖ). The r_max² of r dr appears in H and in N and cancels in
# V = H/N, so the transform takes r = r_max·s and the weights g as they are.
#
# The rule is panelled Gauss–Legendre:
#   - inside s_s = 0.9, GL in s, on panels between the profile's breakpoints;
#   - on the rim [s_s, 1], μ = t², s = √(1 − t⁴), s ds = 2t³ dt, GL in t: a limb-darkened edge
#     μ^α becomes t^{2α}, and the standard laws become polynomials in t;
#   - a centre where the profile is singular (I ∝ s^-p), s = s₁·u^k on [0, s₁], k = 3/(2 − p),
#     which makes I s ds and its derivative in p smooth in u.
# GL in s alone converges slowly on the rim's infinite slope; GL in t alone puts a √ at the
# centre, which odd azimuthal orders feel. Split this way the same count reaches 1e-8 in V,
# ∂V/∂ζ and ∂V/∂α for every order, the five LD laws, a Gaussian and a uniform disc.
#
# The breakpoints come from the profile (`_hankel_probe`): jumps, kinks, and steep but smooth
# transitions, found on a fine grid at the current parameters and located by bisection. Each
# panel is then smooth, and the count converges as for a smooth profile. How many nodes is
# decided by `_hankel_converge`, doubling from the lobe rule `_hankel_n` until V and ∂V/∂p
# agree. PLAN_HankelQuadratureStudy.md, stages 1c, 1d and 1e, has the measurements.

const _HANKEL_SPLIT = 0.9      # s_s: where the inner panels end and the rim starts
const _HANKEL_RIM   = 0.3      # fraction of the nodes given to the rim
const _HANKEL_CAP   = 4096     # the most nodes the convergence check goes to

"""
    _hankel_n(lobes) -> Int

The number of nodes a smooth profile needs for |ΔV| ≤ 1e-8 on baselines up to
`lobes = B_max · D` (D the outer diameter, B in cycles per unit of D). Measured from 1 to 100
lobes, values and derivatives, and short nowhere: 32 nodes up to 7 lobes, then ~1.1 per lobe.
It is where `_hankel_converge` starts.
"""
_hankel_n(lobes::Real) = 4 * cld(ceil(Int, max(32, 1.1 * lobes + 24)), 4)

# `_HankelLayout` (what the probe found) is defined in model_kernels.jl, beside the cache that
# keeps it.
_hankel_moving(l::_HankelLayout) = any(x -> x != 0, l.dbdp)

# The profile at normalised radii s, with μ given (the rim grid knows μ = t² exactly).
function _hankel_profile(profile_fn, pp, r_max, s::AbstractVector, μ::AbstractVector)
    v = profile_fn(r_max .* s, μ, pp...)
    return v isa Number ? fill(v, length(s)) : v
end
_hankel_profile(profile_fn, pp, r_max, s::AbstractVector) =
    _hankel_profile(profile_fn, pp, r_max, s, sqrt.(max.(zero(eltype(s)), 1 .- s .^ 2)))
_hankel_profile1(profile_fn, pp, r_max, s::T) where {T} =
    Float64(first(_hankel_profile(profile_fn, pp, r_max, T[s])))

# The median of a short vector.
_hankel_med(v) = (w = sort(v); w[cld(length(w), 2)])

"""
    _hankel_probe(profile_fn, pp, r_max, s_min, T) -> _HankelLayout

Where the profile is not smooth, at the parameters `pp`, from the profile alone. Two grids,
each tested on its own: 4097 points even in s on [s_min, 0.99], 512 even in t (μ = t²) over
the rim, where s crowds. A **jump** is a |ΔI| above 8× the median of |ΔI| on each side (5
values, the immediate neighbours left out) and above 64·eps(T)·max|I|, located by bisection;
if it does not survive bisection it is a steep but smooth **transition**, bracketed out to
where I is within 1e-10 of the jump of each plateau. A **kink** is a local maximum of |Δ²I|
against the same two-sided medians, located by zooming. A **singular centre** (s_min = 0) is
I(0) not finite, or I ∝ s^-p between s = 1e-4 and 1e-3 with p > 0.05. Then the derivative of
each breakpoint with respect to each parameter, by locating it again at pp ± h.
"""
function _hankel_probe(profile_fn, pp::AbstractVector, r_max::Real, s_min::Real,
                       ::Type{T}) where {T<:AbstractFloat}
    rm = T(r_max); ppT = T[ForwardDiff.value(x) for x in pp]
    f1(s, q = ppT) = _hankel_profile1(profile_fn, q, rm, T(s))
    smin = Float64(s_min)

    sing = false; pexp = 0.0
    if smin == 0
        a, b = f1(1e-4), f1(1e-3)
        if !isfinite(f1(0.0)) || (a > 0 && b > 0 && log10(a / b) > 0.05)
            sing = true
            pexp = (a > 0 && b > 0) ? clamp(log10(a / b), 0.05, 1.9) : 1.0
        end
    end

    grids = Tuple{Vector{Float64},Vector{Float64}}[]          # (s, μ)
    if smin < 0.99
        s1 = collect(range(smin, 0.99, length = 4097))
        push!(grids, (s1, sqrt.(1 .- s1 .^ 2)))
        t = collect(range((1 - 0.99^2)^0.25, 0.0, length = 513)[1:end-1])
        push!(grids, (sqrt.(1 .- t .^ 4), t .^ 2))
    else
        t = collect(range((1 - smin^2)^0.25, 0.0, length = 4097)[1:end-1])
        push!(grids, (sqrt.(1 .- t .^ 4), t .^ 2))
    end

    breaks = Float64[]; kinds = Symbol[]
    for (s, μ) in grids
        # The profile, evaluated in T, held in Float64: exact, and the differences below stay
        # exact. Typed, since the profile is called through `Function` and is not inferred.
        I::Vector{Float64} = Float64.(_hankel_profile(profile_fn, ppT, rm, T.(s), T.(μ)))
        lo = (sing && s[1] == smin) ? 3 : 1
        fin = filter(isfinite, view(I, lo:length(I)))
        isempty(fin) && continue
        scale = maximum(abs, fin)
        scale == 0 && continue
        floor_ = 64 * eps(T) * scale                           # the profile's rounding noise
        d1 = diff(I); d2 = diff(d1)
        function spike(d, k)
            L = [abs(d[j]) for j in max(1, k - 6):k-2 if isfinite(d[j])]
            R = [abs(d[j]) for j in k+2:min(length(d), k + 6) if isfinite(d[j])]
            (isempty(L) || isempty(R)) && return false
            return abs(d[k]) > max(floor_, 8 * max(_hankel_med(L), _hankel_med(R)))
        end
        for k in lo+1:length(d1)
            isfinite(d1[k]) && abs(d1[k]) > 1e-9 * scale || continue
            # cheap first: larger than both values two points away, before the medians
            k > 2 && k < length(d1) - 1 &&
                abs(d1[k]) <= 8 * min(abs(d1[k-2]), abs(d1[k+2])) && continue
            spike(d1, k) || continue
            b, kind, lo_, hi_ = _hankel_locate_jump(f1, s[k], s[k+1],
                                                    s[max(1, k - 1)], s[min(length(s), k + 2)])
            if kind === :jump
                push!(breaks, b); push!(kinds, :jump)
            else
                append!(breaks, (lo_, b, hi_)); append!(kinds, (:transition, :transition, :transition))
            end
        end
        for k in lo+2:length(d2)-1
            abs(d2[k]) >= abs(d2[k-1]) && abs(d2[k]) > abs(d2[k+1]) || continue
            k + 2 <= length(d2) && abs(d2[k]) <= 8 * min(abs(d2[k-2]), abs(d2[k+2])) && continue
            hk = abs(s[k+1] - s[k])
            abs(d2[k]) > 1e-9 * scale && abs(d2[k]) > 1e3 * hk^2 * scale && spike(d2, k) || continue
            sk = s[k+1]
            any(abs(b - sk) < 3 * abs(s[k+2] - s[k]) for b in breaks) && continue
            push!(breaks, _hankel_locate_kink(profile_fn, ppT, rm, minmax(s[k], s[k+2])..., T))
            push!(kinds, :kink)
        end
    end
    p = sortperm(breaks)
    keep = [i == 1 || breaks[p[i]] - breaks[p[i-1]] > 1e-14 for i in eachindex(p)]
    bs = breaks[p][keep]; ks = kinds[p][keep]

    # How each jump and kink moves with each parameter: located again at pp ± h, in a bracket
    # a few probe spacings wide. A transition is left where it is: it only refines panels.
    dbdp = zeros(length(bs), length(ppT))
    w = 4 * (0.99 - min(smin, 0.98)) / 4096
    for (i, b) in enumerate(bs), j in eachindex(ppT)
        ks[i] === :transition && continue
        h = cbrt(eps(T)) * max(1, abs(ppT[j]))     # above T's resolution of the bisection
        bp = _hankel_relocate(profile_fn, setindex!(copy(ppT), ppT[j] + h, j), rm, b, ks[i], w, T)
        bm = _hankel_relocate(profile_fn, setindex!(copy(ppT), ppT[j] - h, j), rm, b, ks[i], w, T)
        (isfinite(bp) && isfinite(bm)) || continue
        d = (bp - bm) / 2h
        dbdp[i, j] = abs(d) > 1e-9 ? d : 0.0
    end
    return _HankelLayout(bs, ks, sing, pexp, dbdp)
end

# A jump between a and b (the probe interval), by bisection. If it does not survive
# bisection it is a steep smooth transition: bracket it, out to where I is within 1e-10 of
# the jump of the plateau on each side (sl, sr: the probe points beyond).
function _hankel_locate_jump(f1, a, b, sl, sr)
    a, b = minmax(a, b); sl, sr = minmax(sl, sr)
    fa, fb = f1(a), f1(b); J = fb - fa
    for _ in 1:60
        c = (a + b) / 2; fc = f1(c)
        abs(fc - fa) > abs(fb - fc) ? (b = c; fb = fc) : (a = c; fa = fc)
    end
    c = (a + b) / 2
    abs(fb - fa) > abs(J) / 2 && return c, :jump, c, c
    Il, Ir = f1(sl), f1(sr)
    lo, hi = sl, c
    for _ in 1:60
        m = (lo + hi) / 2; abs(f1(m) - Il) <= 1e-10 * abs(J) ? (lo = m) : (hi = m)
    end
    left = lo
    lo, hi = c, sr
    for _ in 1:60
        m = (lo + hi) / 2; abs(f1(m) - Ir) <= 1e-10 * abs(J) ? (hi = m) : (lo = m)
    end
    return c, :transition, left, hi
end

# A kink in [a, b]: zoom on the largest second difference, 10 rounds of 65 points.
function _hankel_locate_kink(profile_fn, pp, rm, a, b, ::Type{T}) where {T}
    for _ in 1:10
        z = collect(range(a, b, length = 65))
        fz::Vector{Float64} = Float64.(_hankel_profile(profile_fn, pp, rm, T.(z)))
        j = argmax(abs.(diff(diff(fz)))); w = z[2] - z[1]
        a, b = z[j] - w, z[j+2] + w
    end
    return (a + b) / 2
end

# A jump or kink found at b, located again with the parameters `pp`, in [b − w, b + w], the
# bracket widened ×4 up to 64w if it is not there. NaN if it is still not found (it moved
# further, or went away): the caller probes again.
function _hankel_relocate(profile_fn, pp, rm, b, kind::Symbol, w, ::Type{T}) where {T}
    kind === :transition && return b
    f1(s) = _hankel_profile1(profile_fn, pp, rm, T(s))
    for k in 0:3
        ww = w * 4^k
        a, c = max(0.0, b - ww), min(1.0, b + ww)
        if kind === :jump
            x, kd, _, _ = _hankel_locate_jump(f1, a, c, a, c)
            kd === :jump && return x
        else
            # a kink: the largest second difference inside, away from the bracket's ends
            x = _hankel_locate_kink(profile_fn, pp, rm, a, c, T)
            a + ww / 8 < x < c - ww / 8 && return x
        end
    end
    return NaN
end

"""
    _hankel_breaks_now(layout, profile_fn, pp, r_max, s_min, T) -> breaks or nothing

The layout's breakpoints at the parameters `pp`, those that move located again next to
where they were; `nothing` if one cannot be found there, and the profile must be probed again.
"""
function _hankel_breaks_now(l::_HankelLayout, profile_fn, pp, r_max, s_min, ::Type{T}) where {T}
    _hankel_moving(l) || return l.breaks
    ppT = T[ForwardDiff.value(x) for x in pp]
    w = 4 * (0.99 - min(Float64(s_min), 0.98)) / 4096
    out = copy(l.breaks)
    for i in eachindex(out)
        any(x -> x != 0, view(l.dbdp, i, :)) || continue
        b = _hankel_relocate(profile_fn, ppT, T(r_max), l.breaks[i], l.kinds[i], w, T)
        isfinite(b) || return nothing
        out[i] = b
    end
    return out
end

"""
    _hankel_nodes(T, s_min, n; layout, breaks, n0) -> (s, μ, g)

The panelled rule's nodes on [s_min, 1]: radii s, μ = √(1 − s²) and weights g with
∫ f(s) s ds ≈ Σ gₖ f(sₖ), as `Vector{T}` (or of the breakpoints' type, when they carry Duals:
the nodes then move with the parameters, and the quadrature differentiates the breakpoints'
motion too). Each panel's count is set at `n0` (by length, at least 6) and scaled by n/n0, so
that doubling n refines every panel. Without a layout: no breakpoints, no singular centre.
"""
_hankel_nodes(::Type{T}, s_min::Real, n::Int; layout::Union{Nothing,_HankelLayout} = nothing,
              breaks::AbstractVector = layout === nothing ? Float64[] : layout.breaks,
              n0::Int = n) where {T<:AbstractFloat} =
    _hankel_nodes_impl(T, s_min, n, layout, breaks, n0)

# The body, positional: the keyword call with `breaks` given was not inferred.
function _hankel_nodes_impl(::Type{T}, s_min::Real, n::Int, layout::Union{Nothing,_HankelLayout},
                            breaks::AbstractVector, n0::Int) where {T<:AbstractFloat}
    n >= 2 || throw(ArgumentError("a profile component needs at least 2 nodes, got $n"))
    S, R = _hankel_node_types(T, eltype(breaks))
    sing = layout !== nothing && layout.sing && s_min == 0
    pexp = layout === nothing ? 0.0 : layout.pexp
    smin = Float64(s_min)
    f = n / n0
    ss = max(_HANKEL_SPLIT, smin)
    nrim0 = smin >= _HANKEL_SPLIT ? n0 : max(8, round(Int, _HANKEL_RIM * n0))
    ninner0 = max(8, n0 - nrim0)
    s = S[]; μ = S[]; g = S[]
    # the rim, in t; t runs down as s runs up
    redges = vcat(S(ss), sort([b for b in breaks if ss < b < 1]), S(1))
    ttot = (1 - ss^2)^0.25
    for i in 1:length(redges)-1
        ta, tb = (1 - redges[i+1]^2)^0.25, (1 - redges[i]^2)^0.25
        m = max(2, round(Int, f * max(6, round(Int, nrim0 * ForwardDiff.value(tb - ta) / ttot))))
        x, w = gausslegendre(m)
        t = (tb - ta) / 2 .* (x .+ 1) .+ ta
        append!(s, sqrt.(1 .- t .^ 4)); append!(μ, t .^ 2); append!(g, (tb - ta) .* w .* t .^ 3)
    end
    if smin < _HANKEL_SPLIT
        edges = vcat(S(smin), sort([b for b in breaks if smin < b < ss]), S(ss))
        if sing                       # a singular innermost panel, [0, s₁]
            s1 = min(S(0.1), edges[2])
            edges = vcat(S(0), s1, [e for e in edges[2:end] if e > s1])
        end
        L = ss - smin
        for i in 1:length(edges)-1
            a, b = edges[i], edges[i+1]
            b > a || continue
            m = max(2, round(Int, f * max(6, round(Int, ninner0 * ForwardDiff.value(b - a) / L))))
            x, w = gausslegendre(m)
            if i == 1 && sing
                k = 3 / (2 - pexp)    # I s ds ∝ u², and its p-derivative ∝ u² ln u
                u = (x .+ 1) ./ 2
                sv = b .* u .^ k
                append!(s, sv); append!(g, b^2 * k .* u .^ (2k - 1) .* w ./ 2)
            else
                sv = (b - a) / 2 .* (x .+ 1) .+ a
                append!(s, sv); append!(g, (b - a) / 2 .* w .* sv)
            end
            append!(μ, sqrt.(1 .- sv .^ 2))
        end
    end
    # A node below ∛floatmin(T), or whose weight underflows, in T (next to a singular centre,
    # s = s₁·u^k) contributes nothing T can hold, and would evaluate the profile where it
    # overflows.
    smin_T = cbrt(floatmin(ForwardDiff.valtype(R)))
    keep = [ForwardDiff.value(R(x)) > smin_T && ForwardDiff.value(R(w)) != 0 for (x, w) in zip(s, g)]
    return Vector{R}(s[keep]), Vector{R}(μ[keep]), Vector{R}(g[keep])
end

# The type `_hankel_nodes` computes in and the one it returns: Float64 and T for plain
# breakpoints; with Dual breakpoints (moving with the parameters), their type throughout.
_hankel_node_types(::Type{T}, ::Type{<:AbstractFloat}) where {T} = (Float64, T)
_hankel_node_types(::Type{T}, ::Type{B}) where {T,B} = (promote_type(T, B), promote_type(T, B))

"""
    _hankel_converge(nodes, profile_fn, pp, B, n₀; tol, cap=4096, nsub=64)
        -> (m, agreement, converged)

How many times `n₀` a profile needs on the baselines `B`: V and ∂V/∂p (order 0, the parameters
`pp` at their current values, ∂I/∂p by ForwardDiff) at n and 2n agree to `tol` (∂V/∂p relative
to max(1, |∂V/∂p|)), starting from n = n₀ and doubling. It is judged on `nsub` of the
baselines, spread evenly through them sorted, the longest included: on 11 profiles (steps,
kinks, transitions, singular centres), 622 and 2000 baselines, Float64 and Float32, 64 chose
the count all of them did, in a tenth to a thirtieth of the time. Returns the multiple m of
the smaller count that agreed with the larger one, the last agreement reached, and whether it
reached `tol` before `cap` (at least 2n₀) nodes; when it did not, m·n₀ is the cap.
`nodes(n) -> (r, μ, g)` gives the component's nodes.
"""
function _hankel_converge(nodes, profile_fn, pp::AbstractVector, B::AbstractVector, n0::Int;
                          tol::Real, cap::Int = _HANKEL_CAP, nsub::Int = 64)
    isempty(B) && return (1, 0.0, true)
    Bv = _hankel_subset([abs(ForwardDiff.value(b)) for b in B], nsub)
    ppv = [ForwardDiff.value(x) for x in pp]
    # its own names: an assignment to `V` here would be one to the `V` below
    function vals(n)
        rn, μn, gn = nodes(n)
        Ifun(q) = (v = profile_fn(rn, μn, q...); v isa Number ? fill(v, length(rn)) : v)
        In = Ifun(ppv)
        Kn = _build_hkernel(rn, gn, Bv, 0, false)
        Nn = dot(gn, In)
        Vn = (Kn * In) ./ Nn
        isempty(ppv) && return (Vn, zeros(eltype(Vn), length(Bv), 0))
        dI = ForwardDiff.jacobian(Ifun, ppv)
        return Vn, (Kn * dI .- Vn .* (dI' * gn)') ./ Nn
    end
    capn = max(cap, 2n0)
    m = 1
    V, Vp = vals(n0)
    agree = Inf
    while 2m * n0 <= capn
        V2, Vp2 = vals(2m * n0)
        dp = isempty(Vp2) ? 0.0 : maximum(abs.(Vp2 .- Vp)) / max(1, maximum(abs, Vp2))
        agree = max(Float64(maximum(abs.(V2 .- V))), Float64(dp))
        agree <= tol && return (m, agree, true)
        m *= 2
        V, Vp = V2, Vp2
    end
    return (m, agree, false)
end

# k of the baselines B, evenly through them sorted, the longest included (all, if not more).
function _hankel_subset(B::AbstractVector, k::Int)
    s = sort(unique(B))
    length(s) <= k && return s
    return s[round.(Int, range(1, length(s), length = k))]
end

_hankel_tol(::Type{Float32}) = 1e-6
_hankel_tol(::Type{<:AbstractFloat}) = 1e-8

"""
    _hankel_describe(layout) -> String

What the probe found, for messages: "a jump at r = 3 mas", "a singular centre (I ∝ r^-1.5)".
"""
function _hankel_describe(l::_HankelLayout, r_max)
    parts = String[]
    for (b, k) in zip(l.breaks, l.kinds)
        k === :transition && continue
        push!(parts, "a $(k) at r = $(round(b * r_max; sigdigits = 4)) mas")
    end
    any(==(:transition), l.kinds) && push!(parts, "a steep transition")
    l.sing && push!(parts, "a singular centre (I ∝ r^-$(round(l.pexp; sigdigits = 3)))")
    return isempty(parts) ? "nothing unusual" : join(parts, ", ")
end

"""
    hankel_profile_vis(profile_fn, pp, r_min, r_max, B; nr=0, T=Float64)
        -> (V, n, layout, converged, agreement)

The visibility of a profile (order 0) on the baselines `B` (cycles/mas), with the quadrature a
model component uses: the probe, the panelled rule, and the count from the convergence check
(or `nr` nodes if `nr > 0`). For a preview, where there is a profile and no component.
"""
function hankel_profile_vis(profile_fn, pp::AbstractVector, r_min::Real, r_max::Real,
                            B::AbstractVector; nr::Int = 0, T::Type{<:AbstractFloat} = Float64)
    smin = r_min / r_max
    lay = _hankel_probe(profile_fn, pp, r_max, smin, T)
    n0 = nr > 0 ? nr : _hankel_n(isempty(B) ? 0 : 2 * r_max * maximum(abs, B))
    nodes(n) = ((sn, μn, gn) = _hankel_nodes(T, smin, n; layout = lay, n0); (r_max .* sn, μn, gn))
    m, agree, ok = nr > 0 ? (1, 0.0, true) :
                   _hankel_converge(nodes, profile_fn, pp, B, n0; tol = _hankel_tol(T))
    r, μ, g = nodes(m * n0)
    I = _hankel_profile(profile_fn, T.(pp), T(r_max), r ./ T(r_max), μ)
    return _hankel_wvis(I, r, g, T.(B)), length(r), lay, ok, agree
end
