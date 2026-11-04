# The free-form profile by variational inference (Modeling › Model Optimization › Free-form
# profile): `fit_ld_profile_vi` on a disc of the model with its learned prior, bagged, run on a
# worker like a fit; its I(μ) is drawn as the posterior median with the band of the samples
# under the fits table. Needs VarInf loaded.

# Its canvas, built by `gui()` before the window exists, and the last run.
const LDVI_CANVAS = Ref{Any}(nothing)
const LDVI_RESULT = Ref{Any}(nothing)
const LDVI_START  = Ref{Any}(nothing)        # (μ, I/I(1), I/⟨I⟩) of the profile it started from
const LDVI_DRAW   = Ref((:centre, :mu))      # the normalisation and the abscissa on screen

"""
    build_ldvi_view(fig, ax) -> NamedTuple

The band, the posterior median and the starting profile, empty. **Call before the window is
shown.**
"""
function build_ldvi_view(fig, ax)
    for it in (:scrollzoom, :dragpan, :rectanglezoom, :limitreset)
        Makie.deregister_interaction!(ax, it)
    end
    μ  = Makie.Observable(Float32[0, 1])
    lo = Makie.Observable(Float32[0, 0]); hi = Makie.Observable(Float32[0, 0])
    m  = Makie.Observable(Float32[0, 0])
    μ0 = Makie.Observable(Float32[0, 1]); i0 = Makie.Observable(Float32[0, 0])
    Makie.band!(ax, μ, lo, hi; color = (:royalblue, 0.3))
    Makie.lines!(ax, μ, m; color = :royalblue, linewidth = 2)
    Makie.lines!(ax, μ0, i0; color = :gray30, linestyle = :dash, linewidth = 1.2)
    # No title: the panel's own line above the plot names the run.
    ax.xlabel = "μ"; ax.ylabel = "I(μ)/I(centre)"
    return (; figure = fig, axis = ax, μ, lo, hi, m, μ0, i0)
end

"""
    update_ldvi_view!(panel, r, start, norm, axis = :mu) -> panel

`norm` is `:centre` (I/I(1), the centre at 1) or `:mean` (I/⟨I⟩); `axis` is `:mu` (I(μ)) or
`:r` (I(r), r/R = √(1 − μ²)). The abscissa spans [0, 1]; the ordinate fits the band.
Assignment only, as on every live canvas.
"""
function update_ldvi_view!(p, r, start, norm::Symbol, axis::Symbol = :mu)
    c = norm === :centre
    x(μ) = axis === :r ? sqrt.(max.(0.0, 1 .- μ .^ 2)) : μ
    p.μ[]  = Float32.(x(r.μ))
    p.lo[] = Float32.(c ? r.I_lo : r.F_lo); p.hi[] = Float32.(c ? r.I_hi : r.F_hi)
    p.m[]  = Float32.(c ? r.I_median : r.F_median)
    if start !== nothing
        p.μ0[] = Float32.(x(start[1])); p.i0[] = Float32.(c ? start[2] : start[3])
    end
    v = axis === :r ? "r" : "μ"
    p.axis.xlabel = axis === :r ? "r/R" : "μ"
    p.axis.ylabel = c ? "I($v)/I(centre)" : "I($v)/⟨I⟩"
    p.axis.limits[] = ((0.0, 1.0), nothing)
    Makie.reset_limits!(p.axis)
    return p
end

"""
    shell_ldvi_discs() -> String

The model's discs, by name, one per line: any law, or a free-form disc.
"""
function shell_ldvi_discs()
    m = _model()
    names = try OITOOLS._component_names(m.dict) catch; String[] end
    return join([n for n in names if !isempty(OITOOLS._ld_law_keys(n, m.dict))], "\n")
end

"""
    shell_ldvi_start(component, ngrid, method, ndraws, v2, t3amp, t3phi, scale_errors,
                     bootstrap) -> String

Start `fit_ld_profile_vi` on the current dataset, on a worker, with its default (learned)
prior; `shell_job_poll` follows it. `method` is `geovi` or `mgvi`; closure phases are not
used. `bootstrap` is the number of replicates of bagged VI, 0 for none.
Returns `""` or `!` and a reason.
"""
function shell_ldvi_start(component::AbstractString, ngrid, method::AbstractString, ndraws,
                          v2::Bool, t3amp::Bool, t3phi::Bool, scale_errors::Bool = false,
                          bootstrap = 0)
    sh = _shell()
    sh.job === nothing || return "! something is already running"
    vi_available() || return "! the free-form profile needs VarInf: start Julia with `using VarInf`"
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"
    md = copy(_model(sh).dict)
    cn = String(component)
    isempty(OITOOLS._ld_law_keys(cn, md)) && return "! " * cn * " is not a disc"
    String(method) in ("geovi", "mgvi") || return "! unknown method: " * String(method)
    n = round(Int, Float64(ngrid))
    5 <= n <= 64 || return "! between 5 and 64 nodes"
    nb = round(Int, Float64(bootstrap))
    (nb == 0 || 2 <= nb <= 1000) || return "! bootstrap: 0 or at least 2 replicates"
    # Closure phases are left out: a symmetric disc's are 0 or 180°, which V² already says,
    # and their wrapped residuals are steps the samplers cannot follow.
    weights = fitting_weights(; v2, t3amp, t3phi = false, visamp = false, visphi = false,
                              flux = false, diffvis = false)
    all(iszero, weights) && return "! the free-form profile needs V² or T3amp"
    t3phi && console!(sh, "  closure phases left out of the free-form profile (0 or 180° for a disc)")
    kw = (; ngrid = n, method = Symbol(method), n_draws = round(Int, Float64(ndraws)), weights,
            scale_errors, bootstrap = nb)
    μ = ld_nodes(n)
    I0 = ld_profile(md, cn, μ)
    g0 = OITOOLS._ld_table(:piecewise, n, 0.0).g0
    LDVI_START[] = (μ, I0, I0 ./ (2 * sum(g0 .* I0)))
    data = e.data[1, 1]
    console!(sh, "fit_ld_profile_vi(model, " * _literal(cn) * ", data; ngrid = $n, " *
                 "method = :" * String(method) * ", n_draws = $(kw.n_draws), weights = $(weights), " *
                 "scale_errors = $(scale_errors), bootstrap = $nb)";
             kind = :cmd)
    return start_job!(sh, :ldvi, stop -> fit_ld_profile_vi(md, cn, data; kw..., verb = true))
end

"What a finished run leaves: the result kept, the view drawn, a line for the panel."
function finish_ldvi!(sh::ShellState, res)
    if !res.ok
        msg = "! free-form profile failed: " * _cause(res.err)
        sh.fitlog = _collapse_cr(res.output) * "\n" * msg
        return console!(sh, msg)
    end
    sh.fitlog = _collapse_cr(res.output)
    r = res.result
    LDVI_RESULT[] = r
    LDVI_CANVAS[] === nothing || update_ldvi_view!(LDVI_CANVAS[], r, LDVI_START[], LDVI_DRAW[]...)
    line = @sprintf("%s: diameter %.4f ± %.4f mas, χ²r %.3f, %d samples (%.1f s)", r.component,
                    r.diameter_median, r.diameter_std, r.chi2r, size(r.I_samples, 2), r.seconds)
    r.error_scale == 1 || (line *= @sprintf(", error bars × %.2f", r.error_scale))
    r.bootstrap > 0 && (line *= @sprintf(", bagged over %d replicates (θ ± %.4f on the full data alone)",
                                         r.bootstrap, r.single.diameter_std))
    r.centre === :map && (line *= "; the variational mean diverged, drawn about the maximum")
    console!(sh, "  " * line)
    return line
end

"""
    shell_ldvi_draw(norm, axis = "mu") -> String

Redraw the last run normalised by the centre (`centre`) or by the mean intensity (`mean`),
against μ (`mu`) or r/R (`r`). The choice is kept for the next run.
"""
function shell_ldvi_draw(norm::AbstractString, axis::AbstractString = "mu")
    LDVI_DRAW[] = (String(norm) == "mean" ? :mean : :centre, String(axis) == "r" ? :r : :mu)
    r = LDVI_RESULT[]
    r === nothing && return "! no free-form profile yet"
    LDVI_CANVAS[] === nothing || update_ldvi_view!(LDVI_CANVAS[], r, LDVI_START[], LDVI_DRAW[]...)
    return ""
end

"""
    shell_ldvi_summary() -> String

The last run as `component\\tdiameter\\tσ\\tchi2r\\tmethod\\tprior\\tsamples\\terror scale\\treplicates`,
or `""`; replicates is 0 without bagging.
"""
function shell_ldvi_summary()
    r = LDVI_RESULT[]
    r === nothing && return ""
    return join((r.component, @sprintf("%.5f", r.diameter_median),
                 @sprintf("%.5f", r.diameter_std),
                 @sprintf("%.3f", r.chi2r), String(r.method), String(r.prior),
                 string(size(r.I_samples, 2)), @sprintf("%.3f", r.error_scale), string(r.bootstrap)), '\t')
end

"""
    shell_ldvi_adopt() -> String

Put the last run's profile into the model: the disc becomes a `ldpiecewise` disc on the run's
grid at the variational mean, its nodes fixed, its diameter free if it was.
"""
function shell_ldvi_adopt()
    sh = _shell(); m = _model(sh)
    r = LDVI_RESULT[]
    r === nothing && return "! no free-form profile yet"
    cn = r.component
    old = OITOOLS._ld_law_keys(cn, m.dict)
    isempty(old) && return "! the model has no disc called " * cn
    wasfree = (cn * "," * old[1]) in m.free
    for s in old
        k = cn * "," * s
        delete!(m.dict, k); delete!(m.lb, k); delete!(m.ub, k); filter!(!=(k), m.free)
    end
    for k in r.free
        m.dict[k] = r.model[k]
    end
    wasfree && push!(m.free, r.free[1])
    console!(sh, "  " * cn * ": the free-form profile adopted, " * string(length(r.μ)) *
                 " nodes (fixed), diameter " * @sprintf("%.5f", Float64(r.model[r.free[1]])))
    return ""
end

"""
    shell_ldvi_save(path) -> String

The last run's profile as CSV: μ, r/R, then I/I(1) and I/⟨I⟩, each as median, mean, σ, 16th
and 84th percentiles; the diameter (median ± σ) in a comment line.
"""
function shell_ldvi_save(path::AbstractString)
    r = LDVI_RESULT[]
    r === nothing && return "! no free-form profile yet"
    p = String(path)
    try
        open(p, "w") do io
            @printf(io, "# %s: diameter %.6f ± %.6f mas; %s, %s prior, %d samples; χ²r %.4f; error bars × %.4f; bagged over %d replicates\n",
                    r.component, r.diameter_median, r.diameter_std, r.method, r.prior,
                    size(r.I_samples, 2), r.chi2r, r.error_scale, r.bootstrap)
            println(io, "mu,r_over_R,I,I_mean,I_std,I_16,I_84,I_over_mean,I_over_mean_mean," *
                        "I_over_mean_std,I_over_mean_16,I_over_mean_84")
            for k in eachindex(r.μ)
                @printf(io, "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", r.μ[k],
                        sqrt(max(0.0, 1 - r.μ[k]^2)), r.I_median[k], r.I_mean[k], r.I_std[k],
                        r.I_lo[k], r.I_hi[k], r.F_median[k], r.F_mean[k], r.F_std[k], r.F_lo[k],
                        r.F_hi[k])
            end
        end
    catch err
        return "! " * _cause(err)
    end
    console!(_shell(), "  free-form profile written to " * p)
    return ""
end
