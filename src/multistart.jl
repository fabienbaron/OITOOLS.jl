# Many Levenberg-Marquardt fits of one model from many starting points, and the distinct minima
# they reach — PMOIRED's `gridFit` / `analyseGrid` (`pmoired/oimodels.py`): starts on a grid,
# uniform or normal draws per parameter, conditions that drop starts before fitting, and the
# same rules for setting fits aside and for grouping the rest into minima. An explored parameter
# that is not free is held at its start's value in that fit (PMOIRED's grid parameter in
# `doNotFit`): a χ² profile over it. Here the fits run on threads (one compiled model per thread,
# or per start when parameters are held — compiling again is parsing, the generated code is
# reused), start i draws from `Xoshiro(seed + i)` whatever the thread count, and a component on
# an orbit has its minima merged across whole periods of MJD0.

using Random, Printf, LinearAlgebra, ForwardDiff

"""
    MultiStartMinimum

One distinct minimum of a [`multistart_fit`](@ref): `fit`, the best fit that reached it (an
`LsqFitResult`), `best` the index of its start, and `starts`, the indices of the starts that
ended there.
"""
struct MultiStartMinimum
    fit    :: LsqFitResult
    starts :: Vector{Int}
    best   :: Int
end

"""
    MultiStartResult

What [`multistart_fit`](@ref) returns.

- `minima` — the distinct minima, lowest χ²ᵣ first ([`MultiStartMinimum`](@ref));
  `minima[1].fit` is the global one;
- `starts` — the starting point of every fit (the free parameters then the held ones × starts),
  `list_free_params`, `held` (explored parameters that are not free, held at their start's
  value in each fit) and `explored` (the parameters drawn or gridded);
- `fits` — every fit, `nothing` where it failed; `status` — `:minimum` (grouped into a
  minimum), `:failed` (an error, or χ² not finite), `:no_uncertainties`, `:did_not_move` (the
  explored parameters moved less than 1σ in all), `:chi2` (above the best χ²ᵣ + `deltachi2`);
  `errors` the message of each failure;
- `nrejected` — starts dropped by `constrain` before fitting; `nclamped` — starts moved inside
  the bounds; `seed`, `elapsed` (s), `nthreads`.
"""
struct MultiStartResult
    minima           :: Vector{MultiStartMinimum}
    list_free_params :: Vector{String}
    held             :: Vector{String}
    explored         :: Vector{String}
    starts           :: Matrix{Float64}
    fits             :: Vector{Union{Nothing,LsqFitResult}}
    status           :: Vector{Symbol}
    errors           :: Vector{String}
    nrejected        :: Int
    nclamped         :: Int
    seed             :: UInt64
    elapsed          :: Float64
    nthreads         :: Int
end

function Base.show(io::IO, ::MIME"text/plain", r::MultiStartResult)
    n = length(r.fits)
    nmin = count(==(:minimum), r.status)
    println(io, "MultiStartResult: $n fits on $(r.nthreads) thread(s) in ",
            @sprintf("%.2f s", r.elapsed), r.nrejected > 0 ? " ($(r.nrejected) starts outside `constrain`)" : "")
    println(io, "  fits taken into account: $nmin / $n")
    for (s, what) in ((:failed, "failed"), (:no_uncertainties, "have no uncertainties"),
                      (:did_not_move, "did not move by 1σ from their start"),
                      (:chi2, "are above χ²ᵣ min + deltachi2"))
        k = count(==(s), r.status)
        k > 0 && println(io, "  $k $what")
    end
    if !isempty(r.minima)
        @printf(io, "  distinct minima: %d / %d [~%.1f starts per minimum]\n", length(r.minima),
                nmin, nmin / length(r.minima))
        length(r.minima) > nmin / 2 &&
            println(io, "  too many distinct minima: the grid too coarse or too few starts; the global minimum may not be reliable")
        length(r.minima) < nmin / 4 &&
            println(io, "  few distinct minima: the grid could be coarser, or fewer starts used")
        b = r.minima[1].fit
        @printf(io, "  best: χ²ᵣ %.4f, reached from %d start(s)\n", b.chi2r, length(r.minima[1].starts))
        for (k, v, s) in zip(b.list_free_params, b.x_opt, b.stderror)
            @printf(io, "    %-20s %12.6g ± %.3g\n", k, v, s)
        end
    end
end

# A parameter → its specification, from a Dict (sorted by name) or a vector of pairs (in order).
_ms_pairs(d::AbstractDict) = [String(k) => v for (k, v) in sort(collect(d); by = first)]
_ms_pairs(v::AbstractVector) = [String(first(p)) => last(p) for p in v]
_ms_pairs(::Nothing) = Pair{String,Any}[]

# PMOIRED's grid axis: n = int((max − min)/step + 1) values from min to max.
function _ms_axis(name, spec)
    length(spec) == 3 || throw(ArgumentError("grid $name: (min, max, step), not $(repr(spec))"))
    lo, hi, st = Float64.(spec)
    st > 0 && hi >= lo || throw(ArgumentError("grid $name: needs min ≤ max and step > 0"))
    n = floor(Int, (hi - lo) / st + 1)
    return n == 1 ? [lo] : collect(range(lo, hi; length = n))
end

# The starting points (free parameters × starts), every start inside `constrain`: grid nodes in
# PMOIRED's order (the first gridded parameter varies slowest), each with its uniform and normal
# draws from `Xoshiro(seed + i)`.
function _ms_starts(x0, free, grid, urand, nrand, nstarts, seed)
    idx(n) = (i = findfirst(==(n), free); i === nothing ?
              throw(ArgumentError("$n is explored but is not a parameter of the model")) : i)
    names = [first.(grid); first.(urand); first.(nrand)]
    length(unique(names)) == length(names) ||
        throw(ArgumentError("a parameter is explored once: in grid, rand or randn"))
    for (n, s) in urand
        length(s) == 2 && s[1] <= s[2] || throw(ArgumentError("rand $n: (min, max), not $(repr(s))"))
    end
    for (n, s) in nrand
        length(s) == 2 && s[2] >= 0 || throw(ArgumentError("randn $n: (mean, σ ≥ 0), not $(repr(s))"))
    end
    axes = [_ms_axis(n, s) for (n, s) in grid]
    nodes = isempty(axes) ? [Float64[]] :
            vec([reverse(collect(t)) for t in Iterators.product(reverse(axes)...)])
    N = isempty(axes) ? something(nstarts, 0) : length(nodes)
    isempty(axes) && (nstarts === nothing || nstarts < 1) &&
        throw(ArgumentError("without a grid, give nstarts"))
    gi, ui, ni = idx.(first.(grid)), idx.(first.(urand)), idx.(first.(nrand))
    X = repeat(x0, 1, N)
    for i in 1:N
        isempty(axes) || (X[gi, i] .= nodes[i])
        rng = Xoshiro(seed + i)
        for (k, (_, (lo, hi))) in zip(ui, urand)
            X[k, i] = lo + rand(rng) * (hi - lo)
        end
        for (k, (_, (μ, σ))) in zip(ni, nrand)
            X[k, i] = μ + randn(rng) * σ
        end
    end
    return X, names
end

# Where a component on an orbit has its MJD0 and P among the resolver's values, and the one free
# parameter MJD0 moves with (`(i_MJD0, i_P, k, coefficient)`), or `nothing`: minima whose MJD0
# differ by whole periods are the same orbit.
function _ms_period_shift(fm, model_dict, x)
    cn = try _orbit_component(model_dict) catch; return nothing end
    iM = findfirst(==("$cn,orb MJD0"), fm.all_names); iP = findfirst(==("$cn,orb P"), fm.all_names)
    (iM === nothing || iP === nothing) && return nothing
    g = ForwardDiff.gradient(p -> fm.resolver.fn(p)[iM], x)
    nz = findall(!iszero, g)
    length(nz) == 1 || return nothing
    return (iM, iP, nz[1], g[nz[1]])
end

# A minimum's parameters by name: the fitted free ones and the held ones of its best start.
_ms_values(r, m) = Dict{String,Float64}([zip(r.list_free_params, m.fit.x_opt)...;
                                         zip(r.held, r.starts[(length(r.list_free_params) + 1):end, m.best])...])

# PMOIRED's `analyseGrid`: fits set aside, then the rest grouped into distinct minima. `X` has
# the free parameters' starts then the held ones'; `gstep` is the grid step of each held
# parameter (NaN for a drawn one).
function _ms_analyse(fits, X, free, explored, deltachi2, shift, fm, gstep)
    N = length(fits)
    nf = length(free)
    status = fill(:minimum, N)
    chi2 = [f === nothing ? NaN : f.chi2r for f in fits]
    fin = filter(isfinite, chi2)
    chi2min = isempty(fin) ? NaN : minimum(fin)
    ie = [k for k in (findfirst(==(n), free) for n in explored) if k !== nothing]
    for i in 1:N
        f = fits[i]
        if f === nothing || !isfinite(chi2[i])
            status[i] = :failed
        elseif deltachi2 !== nothing && chi2[i] > chi2min + deltachi2
            status[i] = :chi2
        elseif !any(s -> isfinite(s) && s > 0, f.stderror)
            status[i] = :no_uncertainties
        else
            # the explored parameters moved by less than 1σ in all
            # (when only held parameters are explored, there is nothing to have moved: unlike
            # PMOIRED, which would set every fit aside, the test is not made)
            mv = 0.0
            for k in ie
                s = f.stderror[k]
                isfinite(s) && s > 0 && (mv += (f.x_opt[k] - X[k, i])^2 / s^2)
            end
            !isempty(ie) && mv < 1 && (status[i] = :did_not_move)
        end
    end
    # MJD0 and P of each kept fit, for the period merge
    resolved = shift === nothing ? nothing :
               Dict(i => (fm.resolver.fn(fits[i].x_opt)[shift[1]], fm.resolver.fn(fits[i].x_opt)[shift[2]])
                    for i in 1:N if status[i] === :minimum)
    # the free parameters, then the gridded held ones with σ = step/2 (PMOIRED's: a held drawn
    # parameter takes no part in the distance)
    hg = [h for h in eachindex(gstep) if isfinite(gstep[h])]
    vec_b(i) = [fits[i].x_opt; X[nf .+ hg, i]]
    vec_s(i) = [fits[i].stderror; gstep[hg] ./ 2]
    keep = [i for i in 1:N if status[i] === :minimum]
    remaining = copy(keep)
    minima = MultiStartMinimum[]
    for i in keep
        i in remaining || continue
        bi, si = vec_b(i), vec_s(i)
        nd = length(bi)
        grp = Int[]
        for j in remaining
            bj = vec_b(j)
            if shift !== nothing
                (Mi, Pi), (Mj, _) = resolved[i], resolved[j]
                n = round((Mj - Mi) / Pi)
                n != 0 && (bj = copy(bj); bj[shift[3]] -= n * Pi / shift[4])
            end
            sj = vec_s(j)
            # PMOIRED's distance: the mean over the fitted parameters (NaN terms left out) of
            # (bᵢ − bⱼ)²/(σᵢσⱼ), divided once more by their number
            s, c = 0.0, 0
            for k in 1:nd
                t = (bi[k] - bj[k])^2 / (si[k] * sj[k])
                isnan(t) || (s += t; c += 1)
            end
            c > 0 && (s / c) / nd < 1 && push!(grp, j)
        end
        isempty(grp) && (grp = [i])
        best = grp[argmin([chi2[j] for j in grp])]
        push!(minima, MultiStartMinimum(fits[best], sort(grp), best))
        filter!(j -> j ∉ grp, remaining)
    end
    sort!(minima; by = m -> m.fit.chi2r)
    return minima, status
end

"""
    multistart_fit(model_dict, list_free_params, data; grid, rand, randn, nstarts, constrain,
                   lb, ub, weights, constraints, rv, deltachi2, seed, threaded, maxIter, verbose)
        -> MultiStartResult

Levenberg-Marquardt fits ([`fit_model_lsqfit`](@ref)) from many starting points, and the
distinct minima they reach: PMOIRED's `gridFit`.

The parameters explored, each in one of
- `grid` — `"p" => (min, max, step)`: every node, PMOIRED's count ⌊(max − min)/step + 1⌋ from min
  to max; several grids make every combination (the first varying slowest);
- `rand` — `"p" => (min, max)`: uniform; `randn` — `"p" => (mean, σ)`: normal; one draw per start;

as a Dict (taken by name) or a vector of pairs (in order). With a grid the number of starts is
its node count (draws added to each node); without one, `nstarts`. The other free parameters
start from `model_dict`. An explored parameter that is not free (a fixed number in
`model_dict`) is held at the start's value in that fit — a grid over held parameters is a χ²
profile, every node's fit of the rest. A start outside `lb`/`ub` is moved onto the bound.

`constrain` drops starts before fitting: conditions in the form of `constraints`
(`("sqrt(\$B,x^2 + \$B,y^2)", "<=", 12)`, [`ModelConstraint`](@ref)), held when not violated.
`constraints`, `lb`, `ub`, `weights`, `rv` and `maxIter` apply to every fit, as in
[`fit_model_lsqfit`](@ref).

The fits are then sorted as PMOIRED's `analyseGrid` does: set aside when they failed, have no
uncertainties, moved their explored parameters by less than 1σ in all, or end above the best
χ²ᵣ + `deltachi2` (when given); the rest grouped — two fits are one minimum when the mean over
the free parameters of (bᵢ − bⱼ)²/(σᵢσⱼ), divided by their number, is below 1 (a gridded held
parameter counts, with σ half its step) — each minimum kept as its lowest χ²ᵣ. For a component on an orbit, MJD0 values whole periods apart count as
the same.

Threads (`threaded = true`): one compiled model per thread; start i draws its values from
`Xoshiro(seed + i)`, so the starts and results do not depend on the thread count.
"""
function multistart_fit(model_dict::Dict{String}, list_free_params::Vector{String}, data::OIdata;
                        grid = nothing, rand = nothing, randn = nothing,
                        nstarts::Union{Nothing,Integer} = nothing, constrain = [],
                        lb = Dict{String,Float64}(), ub = Dict{String,Float64}(),
                        weights = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0], constraints = [],
                        rv::Union{Nothing,RVData} = nothing, deltachi2 = nothing,
                        seed::Integer = 1, threaded::Bool = true, maxIter = 200,
                        verbose::Bool = true)
    t0 = time()
    free = list_free_params
    nf = length(free)
    gp, up, np_ = _ms_pairs(grid), _ms_pairs(rand), _ms_pairs(randn)
    held = [n for n in [first.(gp); first.(up); first.(np_)] if n ∉ free]
    for n in held
        get(model_dict, n, nothing) isa Real ||
            throw(ArgumentError("$n is explored but is neither free nor a fixed number in the model"))
    end
    x0 = Float64[model_dict[p] for p in free]
    X, explored = _ms_starts([x0; Float64[model_dict[h] for h in held]], [free; held], gp, up, np_,
                             nstarts, UInt64(seed))
    gstep = [(i = findfirst(p -> first(p) == h, gp); i === nothing ? NaN : Float64(last(gp[i])[3]))
             for h in held]

    # one model with the constraints and the start conditions as extra expressions
    cons, conds = parse_constraints(constraints), parse_constraints(constrain)
    augmented = copy(model_dict)
    ckeys = augment_constraints!(augmented, cons)
    skeys = String[]
    for (i, c) in enumerate(conds)
        k = "__start_$(i)__"; augmented[k] = constraint_expr(c); push!(skeys, k)
    end
    nB = size(data.uv, 2)
    # the model with the held parameters at `hv` (the same generated code whatever the values)
    function compile(hv = Float64[])
        md = augmented
        isempty(held) || (md = copy(augmented); for (h, v) in zip(held, hv); md[h] = v; end)
        return parse_model(md, free; nB_workspace = nB)
    end
    fm = compile()
    cidx = [findfirst(==(k), fm.all_names) for k in ckeys]
    sidx = [findfirst(==(k), fm.all_names) for k in skeys]
    hvals(s) = X[(nf + 1):end, s]
    model_at(s) = isempty(held) ? fm : compile(hvals(s))

    # starts outside `constrain` dropped; inside the bounds
    ok = [isempty(conds) ||
          (m = model_at(s); all(iszero(_violation(c.op, v)) for (c, i) in zip(conds, sidx)
                                for v in _as_penalty_vec(m.resolver.fn(X[1:nf, s])[i])))
          for s in axes(X, 2)]
    nrejected = count(!, ok)
    X = X[:, ok]
    N = size(X, 2)
    N > 0 || throw(ArgumentError("every start is outside `constrain`"))
    lo = Float64[get(lb, p, -Inf) for p in free]; hi = Float64[get(ub, p, Inf) for p in free]
    nclamped = 0
    for s in 1:N
        c = clamp.(view(X, 1:nf, s), lo, hi)
        c == view(X, 1:nf, s) || (nclamped += 1; X[1:nf, s] .= c)
    end
    rvix = rv === nothing ? nothing : _rv_index(fm, _orbit_component(model_dict))

    # the fits, threads taking the next start as they free up
    nw = threaded ? min(Threads.nthreads(), N) : 1
    # one model per thread; with held parameters, one per start, compiled under a lock (the
    # generated-function cache is shared)
    models = isempty(held) ? [w == 1 ? fm : compile() for w in 1:nw] : nothing
    lk = ReentrantLock()
    fits = Vector{Union{Nothing,LsqFitResult}}(nothing, N)
    errors = fill("", N)
    next = Threads.Atomic{Int}(1)
    nblas = BLAS.get_num_threads()
    nw > 1 && BLAS.set_num_threads(1)
    try
        work(w) = while (i = Threads.atomic_add!(next, 1)) <= N
            try
                m = isempty(held) ? models[w] : lock(() -> compile(hvals(i)), lk)
                fits[i] = _lsqfit_compiled(m, free, X[1:nf, i], data, lo, hi; weights, cons,
                                           constraint_idx = cidx, rv, rvix, maxIter,
                                           verb = false, warn_cov = false)
            catch e
                e isa InterruptException && rethrow()
                errors[i] = sprint(showerror, e)
            end
        end
        if nw == 1
            work(1)
        else
            @sync for w in 1:nw
                Threads.@spawn work(w)
            end
        end
    finally
        BLAS.set_num_threads(nblas)
    end

    # (an orbit whose period or periastron depends on a held parameter is not merged across
    # periods: the shift is applied to a free one)
    cn = try _orbit_component(model_dict) catch; "" end
    uses_held(k) = (v = get(model_dict, k, nothing); k in held ||
                    (v isa AbstractString && any(r -> r in held, extract_refs(v))))
    shift = !isempty(cn) && (uses_held("$cn,orb MJD0") || uses_held("$cn,orb P")) ? nothing :
            _ms_period_shift(fm, model_dict, x0)
    minima, status = _ms_analyse(fits, X, free, explored, deltachi2, shift, fm, gstep)
    r = MultiStartResult(minima, copy(free), held, explored, X, fits, status, errors, nrejected,
                         nclamped, UInt64(seed), time() - t0, nw)
    verbose && show(stdout, MIME"text/plain"(), r)
    return r
end
