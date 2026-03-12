"""
MaximENT.jl

Pure-Julia port of MaximENT — Quantified Maximum Entropy Image Reconstruction
(Gull, Charter & Skilling, 1990–1991, Maximum Entropy Data Consultants Ltd).

Original Fortran architecture:
  • A flat workspace array `ST` indexed by 40 numbered "areas" (hidden /
    visible / data spaces), managed through MFETCH/MWHILE/VSTORE blocking
    machinery needed when RAM was scarce.
  • Two Fortran COMMON blocks (`MEMCOM`, `MEINFO`) held all control scalars.
  • OPUS/TROPUS/ICF/MEMEX were global subroutine callbacks.

Julia design:
  • `MaximENTState`  — named Vector fields replacing the numbered area system.
  • `MaximENTParams` — mutable struct replacing the COMMON blocks.
  • User supplies `linfwd!`, `linadj!`, `fwd!` as ordinary Julia functions.
  • The entire MFETCH/MWHILE/VFETCH/VSTORE block-management layer is dropped;
    Julia operates on whole arrays directly.
  • VRAND replaced by Julia's built-in `randn!` / `rand`.
  • MEMDET (bidiagonal eigensolver) replaced by `svd(Bidiagonal(…))`, which
    calls LAPACK's dedicated bidiagonal SVD (`dbdsqr`) directly — more accurate
    and simpler than explicitly forming L·Lᵀ as a `SymTridiagonal`.
  • ICF = identity (as in the original BSMEM application), so hidden-space
    and visible-space have the same dimension.

Entropy modes (methd1):
  1 = Standard          S = Σ(m − h + h·log(h/m))
  2 = Positive/negative (bipolar, for signed maps)
  3 = L1L2w  S = Σ(m·log(1 + h/m) − h)
  4 = Quadratic         S = −Σ h²/(2m)

Likelihood modes (methd2):
  1 = Gaussian (chi-squared)
  2 = Poisson

Stopping criterion (methd0):
  1 = Classic, known noise (σ fixed)
  2 = Classic, unknown noise (σ estimated)
  3 = Alpha given externally
  4 = Historic χ² = N_data
"""

using LinearAlgebra, Random, Printf

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────
const NSIZE = 8      # Maximum entries in the ω(α) interpolation table
const LMAX  = 30     # Maximum conjugate-gradient iterations
const TINY  = 1e-30  # CG overflow guard
const EPS   = 1e-20  # General numerical floor
const FLOOR = 1e-12  # Positivity floor for image update (VUPDT1)
const PILOG = 0.9189385332046727  # log(√(2π))

# ─────────────────────────────────────────────────────────────────────────────
# Control parameters  (replaces ICOM / RCOM COMMON blocks)
# ─────────────────────────────────────────────────────────────────────────────
"""
    MaximENTParams(; kwargs...)

All control parameters for the MaximENT engine.

# Keyword arguments
- `methd`  = [1,1,1]  : [methd0, methd1, methd2] method switches
- `nrand`  = 3        : random vectors used in evidence estimation
- `iseed`  = 0        : RNG seed
- `aim`    = 1.0      : target stopping criterion (AIM ≥ 0)
- `rate`   = 0.5      : dimensionless step-size limit (RATE > 0)
- `utol`   = 0.1      : user tolerance (0 < UTOL < 1)
- `alpha`  = 1.0      : initial regularisation parameter

Per-pixel model and per-datum accuracy vectors are always read from
`s.model` and `s.acc_vec` in the associated `MaximENTState`.
"""
mutable struct MaximENTParams
    # Method switches
    methd0 ::Int    # Stopping criterion
    methd1 ::Int    # Entropy type
    methd2 ::Int    # Likelihood type
    mackay_alpha ::Bool    # Use MacKay fixed-point alpha update when true
    mackay_damp  ::Float64  # Damping exponent for MacKay update (0 < damp ≤ 1)
    ritz_alpha   ::Bool     # Use Ritz-value-based alpha update when true
    # Counters
    nrand  ::Int    # Random vectors for evidence
    iseed  ::Int    # RNG seed
    ntrans ::Int    # Accumulated transform count
    istat  ::Int    # Composite status code
    ntable ::Int    # Current ω(α) table length
    # Scalars
    aim    ::Float64
    rate   ::Float64
    utol   ::Float64   # User tolerance
    tol    ::Float64   # Arithmetic tolerance ≈ machine epsilon
    alpha  ::Float64   # Regularisation parameter
    scale  ::Float64   # Estimated data scale
    hhigh  ::Float64   # Upper bubble bound
    # ω(α) interpolation table
    xtable ::Vector{Float64}   # log(α) values
    ytable ::Vector{Float64}   # log(ω) values
    vtable ::Vector{Float64}   # intrinsic variance values
    # Random number generator
    rng       ::Xoshiro
    saved_rng ::Xoshiro
end

function MaximENTParams(;
        methd ::Vector{Int}  = [1, 1, 1],
        nrand ::Int          = 3,
        iseed ::Int          = 0,
        aim   ::Float64      = 1.0,
        rate  ::Float64      = 0.5,
        utol  ::Float64      = 0.1,
        alpha        ::Float64      = 1.0,
        mackay_alpha ::Bool         = false,
        mackay_damp  ::Float64      = 0.5,
        ritz_alpha   ::Bool         = false)
    # Machine epsilon via the Skilling formula: |(4/3 - 1) - 1/3|
    tol = abs((4.0/3.0 - 1.0) - 1.0/3.0)
    rng = Xoshiro(iseed)
    MaximENTParams(
        methd[1], methd[2], methd[3],
        mackay_alpha, mackay_damp, ritz_alpha,
        nrand, iseed, 0, 0, 0,
        aim, rate, utol, tol, alpha, 1.0, 1.0e20,
        zeros(NSIZE), zeros(NSIZE), fill(EPS, NSIZE),
        rng, copy(rng)
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# State workspace  (replaces the numbered ST area system)
# ─────────────────────────────────────────────────────────────────────────────
"""
    MaximENTState(nhid, ndat)

Working storage for one MaxEnt reconstruction problem.

Areas correspond to the original Fortran numbering:
  Hidden space (size nhid):   <1> sqrt_metric, <2> gradS, <3> model,
                              <4> gradL, <5> h
  Visible space (size nhid):  <11> f  (= h when ICF = identity)
  Data space (size ndat):     <21> data, <22> acc_vec, <23>..<28> workspaces
"""
mutable struct MaximENTState
    sqrt_metric ::Vector{Float64}  # <1>  √[metric]
    gradS       ::Vector{Float64}  # <2>  −∇S  (workspace)
    model       ::Vector{Float64}  # <3>  Default model m
    gradL       ::Vector{Float64}  # <4>  −∇L  (workspace)
    h           ::Vector{Float64}  # <5>  Hidden distribution
    f           ::Vector{Float64}  # <11> Visible distribution f = ICF·h
    data        ::Vector{Float64}  # <21> Observed data D
    acc_vec     ::Vector{Float64}  # <22> Accuracies 1/σ
    d_w1        ::Vector{Float64}  # <23>
    residuals   ::Vector{Float64}  # <24> Normalised residuals
    d_w2        ::Vector{Float64}  # <25>
    d_w3        ::Vector{Float64}  # <26> (random vectors)
    d_w4        ::Vector{Float64}  # <27>
    d_w5        ::Vector{Float64}  # <28>
    ritz_λ      ::Vector{Float64}  # Ritz eigenvalues (lower bidiagonal SVD²)
    ritz_w      ::Vector{Float64}  # Ritz weights (projection × x0sq × trial wt)
end

function MaximENTState(nhid::Int, ndat::Int)
    MaximENTState(
        zeros(nhid), zeros(nhid), zeros(nhid),
        zeros(nhid), zeros(nhid), zeros(nhid),
        zeros(ndat), zeros(ndat), zeros(ndat),
        zeros(ndat), zeros(ndat), zeros(ndat),
        zeros(ndat), zeros(ndat),
        Float64[], Float64[]
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Low-level transform wrappers
# These correspond to Fortran MEMOP / MEMTR / MEMEX, with the ICF = identity.
# ─────────────────────────────────────────────────────────────────────────────

# MEMOP(J,K,FLAG=true): out_d = [acc] · Opus · √[metric] · in_h
function _apply_fwd!(s, p, in_h, out_d, linfwd!)
    tmp = s.sqrt_metric .* in_h
    linfwd!(tmp, out_d)
    out_d .*= s.acc_vec
    p.ntrans += 1
end

# MEMTR(K,J,FLAG=true): out_h = √[metric] · Tropus · [acc] · in_d
function _apply_adj!(s, p, in_d, out_h, linadj!)
    linadj!(in_d .* s.acc_vec, out_h)
    out_h .*= s.sqrt_metric
    p.ntrans += 1
end

# full forward operator: out_d = fwd!(in_h)  (nonlinear experiment, or linfwd! for linear case)
function _apply_full_fwd!(s, p, in_h, out_d, fwd!)
    fwd!(in_h, out_d)
    p.ntrans += 1
end

# ─────────────────────────────────────────────────────────────────────────────
# Entropy functions  (MEMENT / MENTx)
# Fills s.sqrt_metric (<1>) and s.gradS (<2>), returns (S, summet)
# ─────────────────────────────────────────────────────────────────────────────

function entropy_gradient!(s::MaximENTState, p::MaximENTParams)
    if     p.methd1 == 1; return _ment1!(s, p)
    elseif p.methd1 == 2; return _ment2!(s, p)
    elseif p.methd1 == 3; return _ment3!(s, p)
    else;                  return _ment4!(s, p)
    end
end

# Standard entropy: S = Σ(m − h + h·log(h/m))
# √[metric] = √h,   −∇S = log(h/m)
function _ment1!(s, p)
    @. s.gradS       = log(s.h / s.model)
    @. s.sqrt_metric = sqrt(s.h)
    summet = sum(s.h)
    S      = summet - sum(s.model) - dot(s.h, s.gradS)
    return S, summet
end

# Positive/negative (bipolar) entropy (VENT2)
# Works on pairs (h⁺, h⁻) encoded as the single hidden variable h.
@inline function _ment2_element(H, M)
    NTERM = 6
    A = 2M
    X = H / A
    Y = sqrt(1 + X^2)
    Z = if X > -10
        Y + X
    else
        rx2 = 1 / X^2
        Z   = 1.0
        for iterm in NTERM:-1:2
            Z = 1 + (1.5 - iterm) * rx2 * Z / iterm
        end
        -Z / (2X)
    end
    logZ = log(Z)
    return A*Y, logZ, A*(Y - 1 - X*logZ), sqrt(A*Y)   # (AY, gradS, Si, sqrt_metric)
end

function _ment2!(s, p)
    S = 0.0; summet = 0.0
    @inbounds for i in eachindex(s.h)
        AY, logZ, Si, meti = _ment2_element(s.h[i], s.model[i])
        summet           += AY
        s.gradS[i]        = logZ
        S                += Si
        s.sqrt_metric[i]  = meti
    end
    return S, summet
end

# L1L2w:  S = Σ(m·log(1 + h/m) − h)
# √[metric] = (m+h)/√m,   −∇S = h/(m+h)
function _ment3!(s, p)
    mh = s.model .+ s.h
    @. s.gradS       = s.h / mh
    @. s.sqrt_metric = mh / sqrt(s.model)
    summet = dot(s.sqrt_metric, s.sqrt_metric)
    S      = mapreduce((mi, mhi, hi) -> mi * log(mhi/mi) - hi, +, s.model, mh, s.h)
    return S, summet
end

# Quadratic regularisation:  S = −Σ h²/(2m)
# √[metric] = √m,   −∇S = h/m
function _ment4!(s, p)
    @. s.gradS       = s.h / s.model
    @. s.sqrt_metric = sqrt(s.model)
    summet = sum(s.model)
    S      = -dot(s.h, s.gradS) / 2
    return S, summet
end

# ─────────────────────────────────────────────────────────────────────────────
# Likelihood functions  (MEMCHI / MEMPOI)
# Returns (alhood, data_count, units)
# Fills s.residuals (<24>) and optionally s.acc_vec (<22>)
# ─────────────────────────────────────────────────────────────────────────────

# Gaussian: residuals = (D − F)·acc,  alhood = ½ Σ residuals²
function likelihood_gradient!(s::MaximENTState, p::MaximENTParams)
    D = s.data; F = s.d_w2; r = s.residuals
    a = s.acc_vec .+ EPS
    @. r = (D - F) * a
    return dot(r, r) / 2, Float64(count(>(EPS), a)), sum(log, a)
end

# Poisson: updates acc_vec, residuals, alhood
function poisson_gradient!(s::MaximENTState, p::MaximENTParams)
    D   = s.data; F = s.d_w2
    PI2 = 0.15915494309189534   # 1/(2π)
    @. s.acc_vec   = sqrt(D + 1) / (F + EPS)
    @. s.residuals = (D - F) * s.acc_vec
    alhood = mapreduce((d, f) -> f - d + d * log(d / (f + EPS) + EPS), +, D, F)
    units  = -0.5 * sum(d -> log(d + PI2), D)
    return alhood, Float64(length(D)), units
end

# METEST: 1 − cos(∠(−∇S, −∇L))
function gradient_angle(s::MaximENTState)
    ss = dot(s.gradS, s.gradS)
    sc = dot(s.gradS, s.gradL)
    cc = dot(s.gradL, s.gradL)
    denom = sqrt(max(ss * cc, EPS^2))
    return 1.0 - sc / denom
end

# ─────────────────────────────────────────────────────────────────────────────
# ω(α) interpolation table  (MEMLAW / MEMLAR)
# ─────────────────────────────────────────────────────────────────────────────

# Insert (or update) the entry (alpha, omega, var) into the table.
function alpha_table_update!(p::MaximENTParams, alpha::Float64, omega::Float64,
                 var::Float64, init::Bool)
    xnew = log(alpha)
    ynew = log(max(omega, EPS))
    vnew = max(var, EPS)
    init && (p.ntable = 0)
    # Update existing entry if alpha matches
    for i in 1:p.ntable
        if xnew == p.xtable[i]
            s = vnew / (vnew + p.vtable[i])
            x = p.vtable[i] / (vnew + p.vtable[i])
            p.ytable[i] = ynew * x + p.ytable[i] * s
            p.vtable[i] = vnew * x
            return
        end
    end
    # Delete worst entry if table is full
    if p.ntable >= NSIZE
        S = max(var, EPS)
        j = -1
        for i in 1:p.ntable
            x  = log(alpha) - p.xtable[i]
            d  = p.vtable[i] + x^4
            if d > S; S = d; j = i; end
        end
        j == -1 && return
        p.ntable -= 1
        for i in j:p.ntable
            p.xtable[i] = p.xtable[i+1]
            p.ytable[i] = p.ytable[i+1]
            p.vtable[i] = p.vtable[i+1]
        end
    end
    # Insert in sorted order
    j = 1
    for i in 1:p.ntable
        xnew > p.xtable[i] && (j = i + 1)
    end
    for i in p.ntable:-1:j
        p.xtable[i+1] = p.xtable[i]
        p.ytable[i+1] = p.ytable[i]
        p.vtable[i+1] = p.vtable[i]
    end
    p.xtable[j] = xnew; p.ytable[j] = ynew; p.vtable[j] = vnew
    p.ntable += 1
end

# Weighted linear regression on the table → estimate of log(ω(alpha))
function alpha_table_interpolate(p::MaximENTParams, alpha::Float64)
    p.ntable == 1 && return p.ytable[1], 0.0
    # First pass: weighted mean
    W = 0.0; WX = 0.0; WY = 0.0
    for i in 1:p.ntable
        x  = log(alpha) - p.xtable[i]
        wt = 1.0 / max(p.vtable[i] + x^4, p.tol)
        W  += wt; WX += wt * x; WY += wt * p.ytable[i]
    end
    xbar = WX / W; ybar = WY / W
    # Second pass: slope
    WXX = 0.0; WXY = 0.0
    for i in 1:p.ntable
        x  = log(alpha) - p.xtable[i]
        wt = 1.0 / max(p.vtable[i] + x^4, p.tol)
        dx = x - xbar; dy = p.ytable[i] - ybar
        WXX += wt * dx^2; WXY += wt * dx * dy
    end
    yval  = ybar - xbar * (WXX > p.tol ? WXY / WXX : 0.0)
    sigma = 1.0 / sqrt(max(W, EPS))
    return yval, sigma
end

# ─────────────────────────────────────────────────────────────────────────────
# Alpha control  (MEMLA0 / MEMLA)
# ─────────────────────────────────────────────────────────────────────────────

function alpha_init!(p::MaximENTParams, summet, gradl, omega, var)
    alpha_table_update!(p, 1.0e20, omega, var, true)
    if omega < 1.0
        p.alpha = (gradl / sqrt(summet)) / p.rate
        return false   # lcode = false → alpha changed
    else
        p.alpha = 1.0e20
        return true
    end
end

function alpha_update!(p::MaximENTParams, summet, omega, var, gradl, init, lcodeb, lcodet)
    alpha_table_update!(p, p.alpha, omega, var, init)
    lcodeb && lcodet || return true   # lcodea = true (no change needed)
    ynew, sigma = alpha_table_interpolate(p, p.alpha)
    abs(ynew) < sigma && return true  # already converged
    # Bisect to find zero of log(ω(α)) = 0
    R    = 1.0 + p.alpha * sqrt(summet) * p.rate / gradl
    alf1 = p.alpha * R
    alf2 = p.alpha / R
    alf1 = max(alf1, p.alpha); alf2 = min(alf2, p.alpha)
    y1, _ = alpha_table_interpolate(p, alf1)
    if y1 > 0
        p.alpha = alf1; return false
    end
    y2, _ = alpha_table_interpolate(p, alf2)
    if y2 <= 0
        p.alpha = alf2; return false
    end
    r = 1.0
    while r > p.tol
        r /= 2
        p.alpha = (alf1 + alf2) / 2
        ynew, _ = alpha_table_interpolate(p, p.alpha)
        if ynew <= 0; alf1 = p.alpha; y1 = ynew
        else;         alf2 = p.alpha; y2 = ynew
        end
    end
    (y1 - y2 != 0) && (p.alpha = alf1 + (alf2 - alf1) * y1 / (y1 - y2))
    return false
end

# ─────────────────────────────────────────────────────────────────────────────
# Beta control  (MEMLB)
# Chooses β ≥ α/min(1,rate) so that step distance ≤ rate
# ─────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
# mackay_alpha_update!  –  MacKay fixed-point alpha update
# ─────────────────────────────────────────────────────────────────────────────

"""
    mackay_alpha_update!(p, omega; damp=p.mackay_damp) → Bool

MacKay fixed-point update for the regularisation parameter α, operating on
the MemSys4 ω estimate directly.

**Derivation.**  MemSys4 defines ω as

    ω = Good · scale² · aim / (−2 · α · S)

At the fixed point ω = 1, the MacKay condition α = Good·scale²·aim/(−2S)
is satisfied.  Substituting, α_new = α_old · ω.  This avoids the biased
Lanczos estimate of Good (which is systematically low at small CG depth M),
because Good cancels between numerator and denominator of ω.

The update is damped by exponent `damp ∈ (0, 1]`:

    α_new = α_old · ω^damp

`damp = 1` is the full Newton step; `damp = 0.5` (default) is a geometric
half-step in log-space that prevents overshoot.  The fixed point is ω = 1
regardless of `damp`.

Returns `true` if α converged (relative log-change ≤ `p.utol`), `false`
otherwise — the same sense as `lcode3` returned by `alpha_update!`.
"""
function mackay_alpha_update!(p::MaximENTParams, omega::Float64;
                               damp::Float64 = p.mackay_damp)
    alpha_old = p.alpha
    p.alpha   = alpha_old * omega^damp
    return abs(log(p.alpha) - log(alpha_old)) <= p.utol
end

function beta_lower_bound(p::MaximENTParams, d0::Float64)
    bmin = sqrt(d0)
    if bmin > p.alpha / min(1.0, p.rate)
        return bmin, false   # distance penalty active
    else
        return p.alpha / min(1.0, p.rate), true
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Image update within positivity limits (MEMUPD / VUPDT1)
# h_new = max(h + δh,  h_old/5,  FLOOR)
# ─────────────────────────────────────────────────────────────────────────────
function image_update!(s::MaximENTState, p::MaximENTParams, delta_h::Vector{Float64})
    if p.methd1 in (1, 3, 4)   # positive-only entropy
        @. s.h = max(s.h + delta_h, s.h / 5, FLOOR)
    else                        # bipolar (pos/neg)
        s.h .+= delta_h
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Hidden-space CG  (MEMCGY)
# Solves  (β I + A) y = b,  b = √[metric] · (−∇L),  A = Memtr∘Memop
# Updates h ← clamp(h + √[metric]·y)
# Returns (hlow, hhigh, dist, lcodeb, lcodec)
# ─────────────────────────────────────────────────────────────────────────────
function cg_solve!(s::MaximENTState, p::MaximENTParams,
                 summet::Float64, linfwd!, linadj!)
    b    = s.gradL      # <4> will hold b then y then δh
    w2   = s.gradS      # <2> hidden workspace
    w23  = s.d_w1       # <23>
    w25  = s.d_w2       # <25>
    w26  = s.d_w3       # <26>
    w27  = s.d_w4       # <27>
    w28  = s.d_w5       # <28>
    fill!(w23, 0.0); fill!(w27, 0.0)

    # b ← √[metric] · b  (into area <4>)
    b .*= s.sqrt_metric
    temp = dot(b, b)
    if temp < TINY
        return 0.0, 0.0, 0.0, true, false
    end

    scale  = sqrt(temp)
    gamma0 = 1.0           # = temp / scale^2
    gamma  = gamma0

    # Beta control: β ≥ α such that ‖y‖ ≤ rate·√summet
    beta, lcodeb = beta_lower_bound(p, temp / (p.rate^2 * summet))

    # First forward transform: w28 = Memop(b) / scale
    _apply_fwd!(s, p, b, w28, linfwd!)
    w28 ./= scale
    fill!(w26, 0.0)

    delta  = dot(w28, w28)
    phi_b  = beta;  eps_b = 1.0
    delb   = phi_b + delta
    QB     = 1.0 / delb
    V0     = 1.0;  Y0 = V0 / delb
    eps_ab = gamma; phi_ab = 0.0; delab = 0.0; QAB = 0.0

    icode = 3   # default: iteration limit
    for L in 1:LMAX
        # Convergence check
        beta * (1.0 + p.utol) * QB + QAB >= gamma0 && (icode = 0; break)
        delta <= TINY  && (icode = 4; break)

        tmp = 1.0/eps_ab - delab * eps_ab / gamma
        phi_ab += beta / (eps_ab * delta * eps_ab) + tmp * gamma * tmp

        # Update data-space gradient accumulator: w26 -= w28 / delta
        @. w26 -= w28 / delta

        # Backward transform: w2 = Memtr(w26)
        _apply_adj!(s, p, w26, w2, linadj!)

        # w2 += b/scale  and  gamma = ‖w2‖²
        @. w2 += b / scale
        gamma = dot(w2, w2)

        delab  = phi_ab + (gamma / eps_ab) / eps_ab
        eps_ab = gamma / (eps_ab * delab)
        QAB   += 1.0 / delab

        beta * (1.0 + p.utol) * QB + QAB >= gamma0 && (icode = 0; break)
        gamma <= gamma0 * p.tol && (icode = 1; break)

        tmp2   = 1.0 / eps_b
        eps_b  = delta / (eps_b * delb)
        phi_b += beta / (eps_b * gamma * eps_b) +
                 (1.0/eps_b - tmp2) * delta * (1.0/eps_b - tmp2)

        # Forward transform: w25 = Memop(w2)
        _apply_fwd!(s, p, w2, w25, linfwd!)

        # w28 += w25 / gamma  and  delta = ‖w28‖²
        @. w28 += w25 / gamma
        delta = dot(w28, w28)

        delb  = phi_b + (delta / eps_b) / eps_b
        QB   += 1.0 / delb

        tmp   = 1.0 / (eps_b * gamma)
        V0   += tmp
        @. w27 += w26 * tmp

        tmp  = 1.0 / delb
        Y0  += V0 * tmp
        @. w23 += w27 * tmp
    end

    # Reconstruct solution: b (area <4>) ← y = b*Y0 + Memtr(scale·w23)
    # Fortran: MEMSMA(ST,4,Y0,2,4) means <4> := <4>*Y0 + <2>
    # i.e. b = b*Y0 + w2  (NOT b += Y0*w2)
    w23 .*= scale
    _apply_adj!(s, p, w23, w2, linadj!)
    @. b = b * Y0 + w2

    YY   = dot(b, b)
    hlow  = scale^2 * QB / 2.0
    hhigh = scale^2 * (gamma0 - QAB) / (2.0 * beta)
    dist  = sqrt(YY / max(summet, EPS))

    # δh = √[metric] · y; update h
    b .*= s.sqrt_metric
    image_update!(s, p, b)

    return hlow, hhigh, dist, lcodeb, (icode <= 2)
end


# ─────────────────────────────────────────────────────────────────────────────
# Additional transform variants (FLAG=FALSE for MESCAL)
# ─────────────────────────────────────────────────────────────────────────────

# MEMOP FLAG=FALSE: generate mock data (fwd!, no metric)
# Fortran: CALL MEMOP(ST,5,25,.FALSE.)  →  <25> = MemEx(ICF * <5>)
# With ICF=identity: out_d = memex(h)
function _compute_mockdata!(s::MaximENTState, p::MaximENTParams, out_d, fwd!)
    fwd!(s.h, out_d)
    p.ntrans += 1
end

# MEMTR FLAG=FALSE: backward transform without metric scaling
# Fortran: CALL MEMTR(ST,24,4,.FALSE.)  →  <4> = TrICF * Tropus * [acc] * <24>
# Used in MESCAL to get gradL without sqrt_metric factor
function _compute_gradient!(s::MaximENTState, p::MaximENTParams, in_d, out_h, linadj!)
    linadj!(in_d .* s.acc_vec, out_h)
    p.ntrans += 1
end

# ─────────────────────────────────────────────────────────────────────────────
# Data-space conjugate gradient  (MEMCGP)
# Builds Lanczos bidiagonalisation scalars gam[], del[] for evidence calc.
# Operator:  A = Memop ∘ Memtr   (data space)
#
# Input:  s.d_w3 (<26>) = random vector b  (modified in-place)
# Uses:   s.gradS (<2>), s.d_w2 (<25>), s.d_w5 (<28>)
# Returns: (gam[0..M], del[0..N], M, N, lcodec)
#          gam and del are 1-indexed in Julia (Fortran 0-indexed)
# ─────────────────────────────────────────────────────────────────────────────
function lanczos_bidiag!(s::MaximENTState, p::MaximENTParams, linfwd!, linadj!)
    b   = s.d_w3   # <26>  input random vector, modified in-place
    w2  = s.gradS  # <2>   hidden-space workspace
    w25 = s.d_w2   # <25>  data-space workspace
    w28 = s.d_w5   # <28>  conjugate accumulator (data space)

    gam = zeros(LMAX+1)
    del = zeros(LMAX+1)
    M   = 0
    N   = 0

    temp = dot(b, b)
    if temp < TINY
        return gam, del, 0, 0, false
    end

    scale  = sqrt(temp)
    gamma0 = 1.0
    gamma  = gamma0
    gam[1] = scale            # gam[0] in Fortran = sqrt(gamma)*scale = scale

    b ./= scale               # Normalise b (area <26>) in-place

    eps_ab = gamma0
    phi_ab = 0.0
    delab  = 0.0
    QAB    = 0.0
    eps_b  = 1.0
    QB     = 0.0

    # Initial pass: w28 = b/gamma  (=b since gamma=1),  w2 = Memtr(w28)
    # delta = ||w2||^2,  w25 = Memop(w2)
    copy!(w28, b)           # MESMUL(26, 1/gamma, 28) with gamma=1
    _apply_adj!(s, p, w28, w2, linadj!)   # w2 = sqrt_metric * Tropus * [acc] * w28
    delta  = dot(w2, w2)
    _apply_fwd!(s, p, w2, w25, linfwd!)     # w25 = [acc] * Opus * sqrt_metric * w2
    del[1] = sqrt(delta) / scale      # del[0] in Fortran
    phi_b  = p.alpha / gamma          # = alpha

    icode = 3
    for L in 1:LMAX
        # ── QB update ──
        delb = phi_b + (delta / eps_b) / eps_b
        QB  += 1.0 / delb

        # ── Convergence checks ──
        p.utol * QAB + QAB + p.alpha * QB >= gamma0 && (icode = 0; break)
        delta <= TINY                                 && (icode = 4; break)

        # ── For L>1, recompute w25 = Memop(w2) ──
        L > 1 && _apply_fwd!(s, p, w2, w25, linfwd!)

        # ── QAB update scalars ──
        tmp    = 1.0/eps_ab - delab * eps_ab / gamma
        phi_ab += p.alpha / (eps_ab * delta * eps_ab) + tmp * gamma * tmp

        # ── MEMSMD(25, -1/delta, 26, gamma):
        #    b += (-1/delta)*w25;  gamma = dot(b,b) ──
        @. b  -= w25 / delta
        gamma  = dot(b, b)

        delab  = phi_ab + (gamma / eps_ab) / eps_ab
        M      = L
        gam[M+1] = sqrt(gamma) * scale   # gam[M] in Fortran
        eps_ab = gamma / (eps_ab * delab)
        QAB   += 1.0 / delab

        # ── Convergence checks (again after QAB update) ──
        p.utol * QAB + QAB + p.alpha * QB >= gamma0 && (icode = 0; break)
        gamma <= gamma0 * p.tol                      && (icode = 1; break)

        # ── MEMSMA(26, 1/gamma, 28, 28):  w28 += (1/gamma)*b ──
        @. w28 += b / gamma

        # ── eps_b update ──
        tmp2   = 1.0 / eps_b
        eps_b  = delta / (eps_b * delb)
        phi_b += p.alpha / (eps_b * gamma * eps_b) +
                 (1.0/eps_b - tmp2) * delta * (1.0/eps_b - tmp2)

        # ── New backward transform: w2 = Memtr(w28), delta = ||w2||^2 ──
        _apply_adj!(s, p, w28, w2, linadj!)
        delta = dot(w2, w2)
        N     = L
        del[N+1] = sqrt(delta) / scale   # del[N] in Fortran
    end

    lcodec = (icode <= 2)
    return gam, del, M, N, lcodec
end

# ─────────────────────────────────────────────────────────────────────────────
# Hidden-space CG for quantification  (MEMCGH)
# Estimates Q = b^T (αI + A)^{-1} b,  A = Memtr ∘ Memop  (hidden space)
# Optionally builds the output vector t = (αI + A)^{-1/2} b  (if flag=true)
#
# Input:  s.gradL (<4>) = b  (used as workspace, output if flag=false or t)
# Uses:   s.gradS (<2>), s.d_w2 (<25>), s.d_w3 (<26>), s.d_w4 (<27>), s.d_w5 (<28>)
# coeff[0..n] needed only if flag=true
# Returns: (qlo, qhi, lcode, n, gam[0..n], del[0..n])
# ─────────────────────────────────────────────────────────────────────────────
function cg_hessian!(s::MaximENTState, p::MaximENTParams, alpha::Float64, flag::Bool,
                 coeff::Vector{Float64}, linfwd!, linadj!)
    b   = s.gradL   # <4>  input/output
    w2  = s.gradS   # <2>
    w25 = s.d_w2    # <25>
    w26 = s.d_w3    # <26>
    w27 = s.d_w4    # <27>
    w28 = s.d_w5    # <28>

    gam = zeros(LMAX+1)
    del = zeros(LMAX+1)
    N   = 0

    temp = dot(b, b)
    if temp < TINY
        fill!(b, 0.0)
        return 0.0, 0.0, true, 0, gam, del
    end

    scale  = sqrt(temp)
    gamma0 = 1.0
    gamma  = gamma0
    gam[1] = scale              # gam[0] Fortran

    # Normalise b → w2 (Fortran: MESMUL(4, 1/scale, 2))
    @. w2 = b / scale
    fill!(w26, 0.0)
    fill!(w28, 0.0)
    T0 = flag ? coeff[1] : 0.0  # coeff[0] in Fortran (1-indexed here)
    flag && (fill!(w27, 0.0))

    eps_b  = 1.0
    phi_b  = alpha / gamma
    eps_ab = gamma0
    phi_ab = 0.0
    delab  = 0.0
    QB     = 0.0
    QAB    = 0.0

    icode = 3
    L     = 0
    while true
        # ── Forward transform and MEMSMD:
        #    w28 += (1/gamma)*Memop(w2);  delta = dot(w28,w28) ──
        _apply_fwd!(s, p, w2, w25, linfwd!)    # w25 = Memop(w2)
        @. w28 += w25 / gamma
        delta = dot(w28, w28)
        N     = L
        del[N+1] = sqrt(delta) / scale   # del[L] Fortran

        L += 1

        # ── QB update ──
        delb = phi_b + (delta / eps_b) / eps_b
        QB  += 1.0 / delb

        # ── Convergence checks ──
        (1.0 + p.utol) * alpha * QB + QAB >= gamma0 && (icode = 0; break)
        delta <= TINY                                 && (icode = 4; break)
        L >  LMAX                                     && (icode = 3; break)

        # ── QAB scalars ──
        tmp    = 1.0/eps_ab - delab * eps_ab / gamma
        phi_ab += alpha / (eps_ab * delta * eps_ab) + tmp * gamma * tmp

        # ── MEMSMA(28, -1/delta, 26, 26):  w26 += (-1/delta)*w28 ──
        @. w26 -= w28 / delta

        # ── Backward transform: w2 = Memtr(w26) ──
        _apply_adj!(s, p, w26, w2, linadj!)

        # ── MEMSMD(4, 1/scale, 2, gamma):  w2 += b/scale;  gamma = dot(w2,w2) ──
        @. w2 += b / scale
        gamma = dot(w2, w2)
        gam[L+1] = sqrt(gamma) * scale   # gam[L] Fortran

        delab  = phi_ab + (gamma / eps_ab) / eps_ab
        eps_ab = gamma / (eps_ab * delab)
        QAB   += 1.0 / delab

        # ── Convergence checks ──
        (1.0 + p.utol) * alpha * QB + QAB >= gamma0 && (icode = 0; break)
        gamma <= gamma0 * p.tol                      && (icode = 1; break)

        # ── eps_b update ──
        tmp2   = 1.0 / eps_b
        eps_b  = delta / (eps_b * delb)
        phi_b += alpha / (eps_b * gamma * eps_b) +
                 (1.0/eps_b - tmp2) * delta * (1.0/eps_b - tmp2)

        if flag
            # Accumulate: w27 += coeff[L]*w26
            c = L+1 <= length(coeff) ? coeff[L+1] : 0.0
            @. w27 += w26 * c
            T0 += c
        end
    end

    qlo = scale^2 * QB
    qhi = scale^2 * (gamma0 - QAB) / alpha
    lcode = (icode <= 2)

    if flag
        # Build output vector in b (<4>): t = Memtr(w27) + (T0/scale)*b
        _apply_adj!(s, p, w27, w2, linadj!)    # w2 = Memtr(w27)
        T0_s = T0 / scale
        @. w2 += T0_s * b                   # MEMSMD(4, T0/scale, 2, temp)
        temp2 = dot(w2, w2)
        w2  .*= scale                        # MESMUL(2, scale, 2)
        copy!(b, w2)                         # result in area <4>
        lcode = abs(temp2 - qlo) <= p.utol * qlo
    end

    return qlo, qhi, lcode, N, gam, del
end

# ─────────────────────────────────────────────────────────────────────────────
# Eigenstructure of symmetric matrix from bidiagonal L  (replaces MEMDET)
#
# Returns (λ, vecs) where:
#   λ     : eigenvalues of L·Lᵀ (uplo=:L) or Lᵀ·L (uplo=:U), descending
#   vecs  : first m rows of the corresponding eigenvectors
#
# Key identity:  if B = U·Σ·Vᵀ  then
#   eigenvalues of B·Bᵀ = Σ²,  eigenvectors = columns of U
#   eigenvalues of Bᵀ·B = Σ²,  eigenvectors = columns of V
#
# Using svd(Bidiagonal) hits LAPACK's dedicated bidiagonal SVD (dbdsqr),
# which is more accurate than forming L·Lᵀ explicitly.
# ─────────────────────────────────────────────────────────────────────────────
function _bidiag_svd(d::AbstractVector, e::AbstractVector, m::Int, uplo::Symbol)
    F = svd(Bidiagonal(d, e, :L))   # sign of e doesn't affect σ² or |U/V|²
    λ = F.S .^ 2                    # eigenvalues, already descending
    vecs = (uplo === :U ? F.V : F.U)[1:m, :]   # Lᵀ·L → V,  L·Lᵀ → U
    return λ, vecs
end

# ─────────────────────────────────────────────────────────────────────────────
# Random vector generation  (MEVRND)
# NCORR = 0 : random signs (±1) using s.d_w3 (<26>)
# NCORR > 0 : Gaussian with inter-sample correlation
# ─────────────────────────────────────────────────────────────────────────────
function random_vector!(s::MaximENTState, p::MaximENTParams, area_m::Vector{Float64}, ncorr::Int)
    if ncorr <= 0
        # Random signs: generate normal, extract sign
        randn!(p.rng, area_m)
        @. area_m = sign(area_m)
    elseif ncorr == 1
        randn!(p.rng, area_m)
    else
        # Correlated Gaussian across calls: overlap = NCORR samples
        fac = 1.0 / sqrt(Float64(ncorr))
        randn!(p.rng, area_m)
        area_m .*= fac
        saved = copy(p.rng)
        tmp   = similar(area_m)
        for _ in 2:ncorr
            randn!(p.rng, tmp)
            area_m .+= tmp .* fac
        end
        # Restore so next call overlaps by (ncorr-1) samples
        copy!(p.rng, saved)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# alpha*good for Classic MaxEnt at alpha=∞  (MEPRB0)
# Called from MESCAL when MEMRUN=1, methods 1 or 2.
# Returns alfgd = E[||Memtr(rand)||^2]  (trace estimate)
# ─────────────────────────────────────────────────────────────────────────────
function trace_estimate_init!(s::MaximENTState, p::MaximENTParams, linfwd!, linadj!)
    w26 = s.d_w3   # <26>  random vector workspace
    w2  = s.gradS  # <2>   hidden output

    alfgd = 0.0
    for _ in 1:p.nrand
        random_vector!(s, p, w26, 0)          # random signs
        _apply_adj!(s, p, w26, w2, linadj!)  # w2 = sqrt_metric * Tropus * [acc] * w26
        alfgd += dot(w2, w2)
    end
    return alfgd / p.nrand
end

# ─────────────────────────────────────────────────────────────────────────────
# Evidence estimation: Good d.o.f. and log-determinant  (MEPROB)
# Returns (glow, ghigh, gdev, dlow, dhigh, ddev, lcodeg)
# ─────────────────────────────────────────────────────────────────────────────
function evidence_estimate!(s::MaximENTState, p::MaximENTParams, alpha::Float64, utol::Float64,
                 linfwd!, linadj!)
    w26 = s.d_w3   # <26>  random input vector for each trial

    glow  = 0.0; ghigh = 0.0; gdev = 0.0
    dlow  = 0.0; dhigh = 0.0; ddev = 0.0
    weight = 0.0
    lcodeg = false

    # Workspace for lower-bound trace accumulation only
    tol2 = sqrt(p.tol)

    ritz_λ_all = Float64[]
    ritz_w_all = Float64[]

    for i in 1:p.nrand
        # Generate random sign vector in area <26>
        random_vector!(s, p, w26, 0)

        # Data-space Lanczos bidiagonalization
        gam, del, M, N, lg = lanczos_bidiag!(s, p, linfwd!, linadj!)
        w = lg ? 1.0 : 1.0e-6   # pre-compute trial weight

        goodlo = 0.0; goodev = 0.0; detllo = 0.0; detdev = 0.0
        goodhi = 0.0; detlhi = 0.0

        if M > 0
            # ── Lower bound: eigenstructure of Lᵀ·L, first eigenvector row only ──
            d_lo    = gam[1:M] .* del[1:M]           # L diagonal
            e_lo    = M > 1 ? .-gam[2:M] .* del[1:M-1] : Float64[]   # L subdiagonal
            trace   = sum(d_lo)
            λlo, Ulo = _bidiag_svd(d_lo, e_lo, 1, :U)  # Lᵀ·L eigenvectors

            xs     = λlo[1:M] ./ alpha
            ys     = Ulo[1, 1:M] .^ 2 ./ alpha
            log1xs = log1p.(xs)
            goodlo = sum(@. ys / (1 + xs))
            goodev = sum(@. ys / (1 + xs) * (xs / (1 + xs)))
            detllo = sum(@. ifelse(xs > tol2, ys * log1xs / xs, ys))
            detdev = sum(@. ifelse(xs > tol2, ys * log1xs^2 / xs, ys * xs))

            x0sq   = (gam[1] * del[1])^2
            goodlo *= x0sq;  detllo *= x0sq
            goodev  = 2x0sq * goodev
            detdev  = 2x0sq * detdev
            # Accumulate Ritz pairs now that x0sq and w are both known
            append!(ritz_λ_all, λlo[1:M])
            append!(ritz_w_all, x0sq .* Ulo[1, 1:M] .^ 2 .* w)
        end

        # ── Upper bound: augmented bidiagonal with trace-cap tail ──
        n_up   = M + 1
        d_up   = [gam[j] * del[j] for j in 1:n_up]
        e_up   = n_up > 1 ? [.-gam[j+1] * del[j] for j in 1:n_up-1] : Float64[]
        d_up[n_up] = (M > 0) ? trace / p.tol : 0.0   # regularising tail
        λup, Uup = _bidiag_svd(d_up, e_up, 1, :L)    # L·Lᵀ eigenvectors

        xs     = λup ./ alpha
        vs     = Uup[1, :] .^ 2
        log1xs = log1p.(xs)
        goodhi = sum(@. vs * xs / (1 + xs))
        detlhi = sum(@. vs * ifelse(xs > tol2, log1xs, xs))
        goodhi *= gam[1]^2
        detlhi *= gam[1]^2

        # ── Accumulate (with weight 1 if CG converged, SMALL if not) ──
        glow  += goodlo * w
        ghigh += goodhi * w
        gdev  += goodev * w
        dlow  += detllo * w
        dhigh += detlhi * w
        ddev  += detdev * w
        weight += w
        lcodeg = lcodeg || lg
    end

    if weight > 0
        glow  /= weight; ghigh /= weight; gdev = sqrt(gdev) / weight
        dlow  /= weight; dhigh /= weight; ddev = sqrt(ddev) / weight
        ritz_w_all ./= weight
    end
    return glow, ghigh, gdev, dlow, dhigh, ddev, lcodeg, ritz_λ_all, ritz_w_all
end

# ─────────────────────────────────────────────────────────────────────────────
# MESCAL: compute scalar statistics for current h
# Fills s.sqrt_metric, s.gradS, s.gradL, s.residuals
# Returns (S, test, chisq, scale, omega, var, summet, gradl,
#          lcodeh, lcodeo, lcodet, lcodeg,
#          plow, phigh, pdev, glow, ghigh, gdev)
# ─────────────────────────────────────────────────────────────────────────────
function maxent_scalars!(s::MaximENTState, p::MaximENTParams, memrun::Int,
                 linfwd!, linadj!, fwd!)
    DIV(x, y) = x / max(y, EPS)

    # ── Step 1: mock data →  residuals ──
    _compute_mockdata!(s, p, s.d_w2, fwd!)    # <25> = MemEx(h)

    local alhood::Float64, data_cnt::Float64, units::Float64
    if p.methd2 == 1
        alhood, data_cnt, units = likelihood_gradient!(s, p)
    else
        alhood, data_cnt, units = poisson_gradient!(s, p)
    end
    units -= data_cnt * PILOG

    # ── Step 2: gradL = Tropus([acc] * residuals)  (no metric) ──
    _compute_gradient!(s, p, s.residuals, s.gradL, linadj!)

    # ── Step 3: entropy → sqrt_metric, gradS, S, summet ──
    S, summet = entropy_gradient!(s, p)

    if memrun == 1
        S = 0.0
        fill!(s.gradS, 0.0)
    end

    # ── Step 4: TEST ──
    test   = gradient_angle(s)
    lcodet = (test <= 1.0)

    # ── Step 5: GRADL scalar = ‖√[metric] · (−∇L)‖ ──
    alf2s = mapreduce((l, m) -> (l * m)^2, +, s.gradL, s.sqrt_metric)
    gradl = sqrt(alf2s)

    # ── Step 6: SCALE and CHISQ ──
    local alfs::Float64
    if memrun == 1
        alf2s  = -alf2s / 2.0
        alfs   = 0.0
    else
        alfs   = p.alpha * S
    end

    if p.methd0 == 2
        p.scale = sqrt(2.0 * (alhood - alfs) / data_cnt)
    else
        p.scale = 1.0
    end
    chisq = DIV(2.0 * alhood, p.scale^2)

    # ── Step 7: omega and evidence ──
    lcodeh = true
    lcodeg = true
    glow   = 0.0; ghigh = 0.0; gdev = 0.0
    dlow   = 0.0; dhigh = 0.0; ddev = 0.0
    plow   = 0.0; phigh = 0.0; pdev = 0.0
    alfgd  = 0.0
    omega  = 0.0
    var    = 0.0

    if p.methd0 in (1, 2, 3)
        if memrun == 1
            alfgd = trace_estimate_init!(s, p, linfwd!, linadj!)
        else
            glow, ghigh, gdev, dlow, dhigh, ddev, lcodeg, rλ, rw =
                evidence_estimate!(s, p, p.alpha, p.utol, linfwd!, linadj!)
            s.ritz_λ = rλ
            s.ritz_w = rw
            lcodeh = p.hhigh / (p.scale^2) <= p.utol * glow
        end

        scl2 = p.scale^2
        plow  = units - data_cnt * log(p.scale) - dhigh/2.0 +
                alfs/scl2 - chisq/2.0
        phigh = units - data_cnt * log(p.scale) - dlow/2.0  +
                alfs/scl2 - chisq/2.0
        probl = (plow + phigh) / 2.0
        pdev  = ddev / 2.0

        good = (glow + ghigh) / 2.0

        if p.methd0 == 3
            omega = (memrun == 1) ? 0.0 : DIV(p.aim, p.alpha)
            var   = p.utol^2 / 12.0
        else
            if memrun == 1
                omega = DIV(alfgd * scl2 * p.aim, -2.0 * alf2s)
                var   = DIV(test, 1.0 - test) + p.utol^2 / 12.0
            elseif memrun == 2
                omega = DIV(good * scl2 * p.aim, -2.0 * alfs)
                var   = DIV(test, 1.0 - test) + p.utol^2 / 12.0
            end
        end
    elseif p.methd0 == 4
        omega = DIV(data_cnt * p.aim, chisq)
        var   = DIV(test, 1.0 - test) + p.utol^2 / 12.0
    end

    lcodeo = abs(omega - 1.0) <= p.utol

    # ── Step 8: restore gradS (recompute entropy, or zero for memrun==1) ──
    if memrun == 1
        S = 0.0
        fill!(s.gradS, 0.0)
    else
        S, summet = entropy_gradient!(s, p)   # restores s.gradS and s.sqrt_metric
    end

    return (S=S, test=test, chisq=chisq, scale=p.scale,
            omega=omega, var=var, summet=summet, gradl=gradl,
            lcodeh=lcodeh, lcodeo=lcodeo, lcodet=lcodet, lcodeg=lcodeg,
            plow=plow, phigh=phigh, pdev=pdev,
            glow=glow, ghigh=ghigh, gdev=gdev, alhood=alhood)
end

# ─────────────────────────────────────────────────────────────────────────────
# Ritz-value-based MacKay alpha update  (ritz_alpha path)
#
# Uses cached Ritz pairs (s.ritz_λ, s.ritz_w) from the last evidence_estimate! call to
# solve the MacKay fixed-point equation
#
#   Good(α) = −2αS / (scale²·aim)
#
# via bisection, without any new operator calls.  The lower-bound Good estimate
# is  Good(α) ≈ Σⱼ ritz_w[j] / (α + ritz_λ[j]).
#
# Returns true when |log α_new − log α_old| ≤ utol.
# ─────────────────────────────────────────────────────────────────────────────
function ritz_alpha_update!(p::MaximENTParams, s::MaximENTState, S::Float64)
    isempty(s.ritz_λ) && return true   # no Ritz data yet
    abs(S) < EPS      && return true   # entropy ≈ 0, can't solve

    # Lower-bound Good as function of α
    good(α) = sum(w / (α + λ) for (w, λ) in zip(s.ritz_w, s.ritz_λ))
    rhs(α)  = -2.0 * α * S / (p.scale^2 * p.aim)
    f(α)    = good(α) - rhs(α)

    # Bisect: f decreases monotonically (Good↓, rhs↑), unique root
    α_lo = EPS
    α_hi = maximum(s.ritz_λ) * 1.0e8
    # Extend upper bracket if needed
    f_hi = f(α_hi)
    while f_hi >= 0.0
        α_hi *= 1.0e4
        f_hi  = f(α_hi)
        α_hi > 1.0e40 && break
    end
    f(α_lo) <= 0.0 && return true   # degenerate: Good ≤ rhs even at α ≈ 0

    for _ in 1:64
        α_mid = sqrt(α_lo * α_hi)
        f(α_mid) > 0.0 ? (α_lo = α_mid) : (α_hi = α_mid)
        (α_hi / max(α_lo, EPS) - 1.0) < 1.0e-6 && break
    end

    α_old   = p.alpha
    p.alpha = sqrt(α_lo * α_hi)
    return abs(log(p.alpha) - log(α_old)) <= p.utol
end

# ─────────────────────────────────────────────────────────────────────────────
# MEM4: one MaxEnt iterate
#
# memrun:  1 = start (initialise h from model)
#          2 = continue without alpha re-init
#          3 = continue with alpha re-init
#          4 = evaluate only (no step)
#
# Returns a NamedTuple with all output statistics.
# ─────────────────────────────────────────────────────────────────────────────
function maxent_step!(s::MaximENTState, p::MaximENTParams, memrun::Int,
               linfwd!, linadj!, fwd!)
    @assert 1 <= memrun <= 4 "memrun must be 1..4"

    # ── Initialise h for MEMRUN=1 ──
    if memrun == 1
        copy!(s.h, s.model)           # MEHSET: h := model
        p.methd2 == 2 && fill!(s.acc_vec, 0.0)  # MEZERO(22) for Poisson
        p.istat  = 0
        p.ntrans = 0
    else
        # copy h (area 5) from previous result (already in s.h)
        # (Fortran: MECOPY(1,5) → already done by previous iterate exit)
    end

    # ── Compute statistics for current h ──
    st = maxent_scalars!(s, p, memrun, linfwd!, linadj!, fwd!)

    lcode2 = (p.istat ÷ 4 % 2 == 0)   # previous beta-control flag
    p.istat = 0
    st.lcodeh || (p.istat += 1)
    st.lcodeo || (p.istat += 2)
    st.lcodet || (p.istat += 16)
    st.lcodeg || (p.istat += 64)

    omega  = st.omega
    lcode0 = st.lcodeh
    lcode1 = st.lcodeo
    lcode4 = st.lcodet
    lcode6 = st.lcodeg

    # ── Override omega on startup ──
    if memrun == 1 && omega >= 1.0
        omega = 0.5
    end

    lcode3 = true   # alpha control diagnostic
    lcode5 = true   # CG diagnostic
    hlow   = 0.0
    hhigh  = 0.0
    dist   = 0.0

    if memrun != 4
        # ── Alpha control ──
        if memrun == 1
            lcode3 = alpha_init!(p, st.summet, st.gradl, omega, st.var)
        elseif p.ritz_alpha
            lcode3 = ritz_alpha_update!(p, s, st.S)
        elseif p.mackay_alpha
            lcode3 = mackay_alpha_update!(p, st.omega; damp=p.mackay_damp)
        elseif memrun == 2
            lcode3 = alpha_update!(p, st.summet, omega, st.var, st.gradl, false,
                            lcode2, lcode4)
        elseif memrun == 3
            lcode3 = alpha_update!(p, st.summet, omega, st.var, st.gradl, true,
                            lcode2, lcode4)
        end

        # ── b = alpha*gradS - gradL  (Fortran: MEMSMA(2, -alpha, 4, 4)) ──
        # area4 = (-alpha)*area2 + area4  = alpha*gradS - gradL
        @. s.gradL = -p.alpha * s.gradS + s.gradL

        # ── Hidden-space CG step: updates h ──
        hlow, hhigh, dist, lcodeb, lcode5 =
            cg_solve!(s, p, st.summet, linfwd!, linadj!)

        p.hhigh = hhigh

        lcodeb || (p.istat += 4)
        lcode3 || (p.istat += 8)
        lcode5 || (p.istat += 32)
    end

    # ── Update visible image f = ICF * h  (ICF = identity) ──
    copy!(s.f, s.h)   # MEMICF(1,11): f = ICF * h

    ntrnsx = p.ntrans

    return (S=st.S, test=st.test, chisq=st.chisq, scale=st.scale,
            omega=omega, alpha=p.alpha,
            plow=st.plow, phigh=st.phigh, pdev=st.pdev,
            glow=st.glow, ghigh=st.ghigh, gdev=st.gdev,
            istat=p.istat, ntrans=ntrnsx,
            hlow=hlow, hhigh=hhigh, dist=dist, alhood=st.alhood)
end

# ─────────────────────────────────────────────────────────────────────────────
# reconstruct!: outer loop calling maxent_step! repeatedly
#
# Usage:
#   p  = MaximENTParams(methd=[1,1,1], aim=1.0, ...)
#   s  = MaximENTState(nhid, ndat)
#   s.model .= prior_image
#   s.data  .= observed_data
#   s.acc_vec .= 1.0 ./ sigma_data
#   reconstruct!(s, p, 50, linfwd!, linadj!, fwd!; verbose=true)
#   result_image = s.f   (or s.h)
#
# linfwd!(in_h, out_d)    : visible → data transform  (linear part)
# linadj!(in_d, out_h)  : data → visible transpose
# fwd!(in_h, out_d)   : full forward model (may be nonlinear)
#
# Returns a vector of per-iteration named tuples.
# ─────────────────────────────────────────────────────────────────────────────
function reconstruct!(s::MaximENTState, p::MaximENTParams, maxiter::Int,
                      linfwd!, linadj!, fwd!;
                      verbose::Bool = false,
                      ndata::Int = count(>(0.0), s.acc_vec))
    history = []
    for k in 1:maxiter
        memrun = (k == 1) ? 1 : 2
        result = maxent_step!(s, p, memrun, linfwd!, linadj!, fwd!)
        push!(history, result)

        if verbose
            icode_str = string(result.istat, base=2, pad=7)
            rchi2 = result.chisq / max(ndata, 1)
            @printf("  Iter %3d  S=%10.4g  χ²/N=%8.4f  ω=%8.5f  α=%10.4g  code=%s\n",
                    k, result.S, rchi2, result.omega, result.alpha, icode_str)
        end

        # Converged: omega ≈ 1 and no status flags set
        if result.omega >= (1.0 - p.utol) && result.istat == 0
            verbose && println("  Converged at iteration $k")
            break
        end
    end
    return history
end

