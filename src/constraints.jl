# constraints.jl
#
# Constraints on a parametric fit: relations between parameters that bounds cannot express.
#
# Load order (within OITOOLS module):
#   include("parse_model.jl")        # FlatModel, the resolver, default_bounds
#   include("constraints.jl")        # this file
#   include("fit_model.jl")          # the fitters, which apply what is defined here
#
# Why this is not bounds
# ----------------------
# `lb`/`ub` describe an axis-aligned box, so they can say "diamout is between 0 and 40" but not
# "diamout is larger than diamin". Every fitter here takes the box natively — NLopt, LsqFit and
# UltraNest all consume it — and none of them takes a relation as a bound.
#
# How a constraint is evaluated
# -----------------------------
# The resolver is the only thing that can evaluate an expression over parameter names, so each
# constraint becomes an extra parameter of the compiled model, holding the signed distance
#
#     d = lhs - rhs
#
# with both sides free to name parameters or be whole expressions. `augment_constraints!` adds
# them; the fitters look up where they landed in `all_names`. This is the mechanism `priors`
# already uses.
#
# How a constraint is enforced — two mechanisms, and the difference matters
# ------------------------------------------------------------------------
# `fit_model` hands them to **NLopt as real nonlinear constraints**, so `diamout > diamin`
# holds at the optimum. `constrained_opt` and `add_nlopt_constraints!` are that path.
#
# `fit_model_lsqfit` and `fit_model_ultranest` have no such machinery, and use a **one-sided
# quadratic penalty** on the normalised violation:
#
#     r = d / |tol|
#     penalty += (the part of r that violates the operator)^2
#
# which is what PMOIRED's `prior` list does, so a model converts between the two packages
# without changing meaning. `tol` is the softness scale, not a threshold: the penalty reaches 1
# when the constraint is violated by `tol`, so a smaller `tol` is stiffer. A satisfied
# constraint contributes exactly zero and exactly zero gradient, and so cannot pull a fit that
# already obeys it.
#
# A penalty is soft, and soft is not a detail. Measured on `demos/data/AlphaCenA.oifits`,
# capping a uniform-disc diameter 0.5 mas below its best fit costs Δχ² ≈ 1.5e6; a penalty stiff
# enough to win against that would wreck the conditioning, and one that is not simply loses.
# The fit then settles where the two gradients balance — neither the constrained answer nor the
# unconstrained one. That is why `fit_model` does not use the penalty, and why the two other
# fitters say in their docstrings that theirs is soft.

using ForwardDiff, NLopt

"""
Default tolerance for a constraint, in the units of the quantity being constrained.

It means two related things, one per enforcement mechanism. To NLopt it is the feasibility
tolerance: how far outside the constraint a point may sit and still count as satisfied. To the
penalty used by the least-squares and nested-sampling fitters it is the stiffness: a violation
of `tol` costs one unit of penalty.

`1e-3` is a thousandth of a milliarcsecond on an angular size, and a tenth of a percent on a
flux fraction — below anything interferometry resolves, in both readings.
"""
const DEFAULT_CONSTRAINT_TOL = 1e-3

const _CONSTRAINT_OPS = Dict(
    "<"  => :<,  "<=" => :≤, "≤"  => :≤,
    ">"  => :>,  ">=" => :≥, "≥"  => :≥,
    "="  => :(=), "==" => :(=),
)

"""
    ModelConstraint(lhs, op, rhs; tol = DEFAULT_CONSTRAINT_TOL)

A relation between model parameters, enforced during fitting.

`lhs` and `rhs` are parameter names (`"disk,diamout"`), expressions over parameter names
(`"2 * star,ud"`), or — for `rhs` — a plain number. `op` is one of `"<"`, `"<="`, `">"`,
`">="`, `"="`, given as a string or a symbol. `tol` is how much violation counts as none; see
[`DEFAULT_CONSTRAINT_TOL`](@ref).

[`fit_model`](@ref) hands these to NLopt as real nonlinear constraints, so they hold at the
optimum. [`fit_model_lsqfit`](@ref) and [`fit_model_ultranest`](@ref) have no such machinery
and apply a one-sided quadratic penalty instead, which is soft: a steep enough χ² can overrule
it, and `tol` then sets the stiffness rather than a feasibility tolerance.

```julia
ModelConstraint("disk,diamout", ">", "disk,diamin")   # a ring must have an outside
ModelConstraint("star,f", "+", 0.0)                   # ── invalid: `+` is not a relation
ModelConstraint("star,f + disk,f", "=", 1.0)          # fluxes sum to one
ModelConstraint("star,ud", "<", 4.0; tol = 0.01)      # deliberately loose
```

Globals — model keys with no comma — must be written with a `\\\$`, as `"\\\$PA"`, exactly as in
a model expression. A bare name with no comma and no arithmetic is taken to be one and gets
its `\\\$` added.

See also [`fit_model`](@ref) and [`read_model_file`](@ref).
"""
struct ModelConstraint
    lhs :: String
    op  :: Symbol
    rhs :: Union{Float64,String}
    tol :: Float64

    function ModelConstraint(lhs, op, rhs; tol = DEFAULT_CONSTRAINT_TOL)
        o = op isa Symbol ? op : get(_CONSTRAINT_OPS, String(op)) do
            error("Unknown constraint operator '$op'; expected one of " *
                  join(sort(collect(keys(_CONSTRAINT_OPS))), ", "))
        end
        o = get(_CONSTRAINT_OPS, string(o), o)
        o in (:<, :≤, :>, :≥, :(=)) ||
            error("Unknown constraint operator '$op'; expected one of " *
                  join(sort(collect(keys(_CONSTRAINT_OPS))), ", "))
        tol = Float64(tol)
        tol > 0 || error("Constraint tolerance must be positive, got $tol")
        new(String(lhs), o, rhs isa Number ? Float64(rhs) : String(rhs), tol)
    end
end

_op_string(op::Symbol) = op === :≤ ? "<=" : op === :≥ ? ">=" : op === :(=) ? "=" : string(op)

function Base.show(io::IO, c::ModelConstraint)
    rhs = c.rhs isa Float64 ? string(c.rhs) : c.rhs
    print(io, "ModelConstraint(", c.lhs, " ", _op_string(c.op), " ", rhs,
              ", tol=", c.tol, ")")
end

"""
    parse_constraints(specs) -> Vector{ModelConstraint}

Accept constraints in any of the forms a caller is likely to have them in.

A `ModelConstraint` passes through. A tuple is read positionally as
`(lhs, op, rhs)` or `(lhs, op, rhs, tol)` — the PMOIRED `prior` layout, so a converted PMOIRED
model needs no reshaping. A `Dict` or `NamedTuple` is read by the keys `param`/`lhs`, `op`,
`value`/`rhs` and `tol`, which is the shape a TOML file parses to.
"""
parse_constraints(specs) = ModelConstraint[_as_constraint(s) for s in specs]
parse_constraints(::Nothing) = ModelConstraint[]

_as_constraint(c::ModelConstraint) = c

function _as_constraint(t::Union{Tuple,AbstractVector})
    length(t) == 3 && return ModelConstraint(t[1], t[2], t[3])
    length(t) == 4 && return ModelConstraint(t[1], t[2], t[3]; tol = t[4])
    error("A constraint tuple must be (lhs, op, rhs) or (lhs, op, rhs, tol), got $(length(t)) entries")
end

function _as_constraint(d::AbstractDict)
    get2(a, b) = haskey(d, a) ? d[a] : haskey(d, b) ? d[b] :
                 error("A constraint needs a '$a' (or '$b'); got keys $(collect(keys(d)))")
    lhs = get2("param", "lhs")
    rhs = get2("value", "rhs")
    haskey(d, "op") || error("A constraint needs an 'op'; got keys $(collect(keys(d)))")
    ModelConstraint(lhs, d["op"], rhs; tol = get(d, "tol", DEFAULT_CONSTRAINT_TOL))
end

# NamedTuples index by Symbol, so they take their own accessor rather than stringly-typing the
# one above.
function _as_constraint(d::NamedTuple)
    lhs = haskey(d, :param) ? d.param : haskey(d, :lhs) ? d.lhs :
          error("A constraint needs a 'param' (or 'lhs'); got fields $(keys(d))")
    rhs = haskey(d, :value) ? d.value : haskey(d, :rhs) ? d.rhs :
          error("A constraint needs a 'value' (or 'rhs'); got fields $(keys(d))")
    haskey(d, :op) || error("A constraint needs an 'op'; got fields $(keys(d))")
    ModelConstraint(lhs, d.op, rhs; tol = get(d, :tol, DEFAULT_CONSTRAINT_TOL))
end


# ─────────────────────────────────────────────────────────────────────────────
# Turning a constraint into something the resolver can evaluate
# ─────────────────────────────────────────────────────────────────────────────

"""
Reference one side of a constraint the way a model expression would.

`"star,ud"` becomes `"\\\$star,ud"`, `"2*star,ud + disk,f"` gets a `\\\$` on each parameter, and a
bare `"PA"` is read as a global and gets one too. A side that already contains a `\\\$` is left
exactly as written, so a caller who knows the expression language keeps full control.
"""
function _constraint_ref(side::AbstractString)
    s = String(side)
    occursin('$', s) && return s
    d = _dollarify(s)
    # A lone identifier has no comma, so `_dollarify` left it alone; in a model dict that is a
    # global, and referencing one takes the same `$` as anything else.
    occursin(r"^\s*[A-Za-z_][A-Za-z0-9_]*\s*$", d) ? "\$" * strip(d) : d
end

"The expression whose value is `lhs - rhs`, i.e. the signed distance from the constraint."
function constraint_expr(c::ModelConstraint)
    rhs = c.rhs isa Float64 ? repr(c.rhs) : _constraint_ref(c.rhs)
    "(" * _constraint_ref(c.lhs) * ") - (" * rhs * ")"
end

"""
The part of a normalised distance `r` that violates `op`, and zero when it does not.

`min`/`max` rather than a branch so the result carries `r`'s type: these run on ForwardDiff
duals during the gradient pass.
"""
_violation(op::Symbol, r) =
    op === :> || op === :≥ ? min(r, zero(r)) :
    op === :< || op === :≤ ? max(r, zero(r)) :
    r                                          # :(=) — both directions cost

# A constraint on a chromatic parameter is a constraint at every wavelength: the resolver
# broadcasts, so the value arrives as a vector and every element is penalised. Scalars are
# wrapped rather than special-cased so both go through one code path, including under
# ForwardDiff, where the element type is a dual and must be preserved.
_as_penalty_vec(v::Number)        = [v]
_as_penalty_vec(v::AbstractArray) = vec(v)


# ─────────────────────────────────────────────────────────────────────────────
# The penalties, in the three shapes the fitters need
# ─────────────────────────────────────────────────────────────────────────────

"""
    constraint_penalty(x, constraints, model, idx) -> Float64

Total penalty at `x`. `idx[i]` is the position of constraint `i`'s difference expression in
`model.all_names`; [`augment_constraints!`](@ref) produces both.
"""
function constraint_penalty(x, constraints, model, idx)
    isempty(constraints) && return 0.0
    pen = 0.0
    vals = model.resolver.fn(x)
    for (i, c) in enumerate(constraints)
        r = _as_penalty_vec(vals[idx[i]]) ./ c.tol
        pen += sum(abs2, _violation.(c.op, r))
    end
    return pen
end

"""
    constraint_penalty_and_grad!(g, x, constraints, model, idx) -> Float64

Total penalty at `x`, accumulating its gradient into `g`.

A satisfied constraint adds nothing to either, and is skipped before the Jacobian is taken —
which is what keeps constraints affordable, since in a converged fit most of them are inactive.
"""
function constraint_penalty_and_grad!(g, x, constraints, model, idx)
    isempty(constraints) && return 0.0
    pen  = 0.0
    vals = model.resolver.fn(x)
    for (i, c) in enumerate(constraints)
        r = _as_penalty_vec(vals[idx[i]]) ./ c.tol
        v = _violation.(c.op, r)
        all(iszero, v) && continue
        pen += sum(abs2, v)
        J = ForwardDiff.jacobian(x_ -> _as_penalty_vec(model.resolver.fn(x_)[idx[i]]), x)
        g .+= (2 / c.tol) .* (J' * v)
    end
    return pen
end

"""
    constraint_residuals(x, constraints, model, idx) -> Vector{Float64}

The violations as residuals, for a least-squares fitter that sums their squares itself.

One entry per constraint for a scalar model, more when a chromatic parameter makes the
resolver broadcast. The count is fixed for a given model, which is what `curve_fit` requires
of its residual vector.
"""
function constraint_residuals(x, constraints, model, idx)
    isempty(constraints) && return Float64[]
    vals = model.resolver.fn(x)
    return reduce(vcat, (_violation.(c.op, _as_penalty_vec(vals[idx[i]]) ./ c.tol)
                         for (i, c) in enumerate(constraints)))
end

"""
    constraint_jacobian(x, constraints, model, idx) -> Matrix{Float64}

Rows of ∂residual/∂x matching [`constraint_residuals`](@ref), stacked in the same order.
"""
function constraint_jacobian(x, constraints, model, idx)
    isempty(constraints) && return zeros(0, length(x))
    blocks = map(enumerate(constraints)) do (i, c)
        r = _as_penalty_vec(model.resolver.fn(x)[idx[i]]) ./ c.tol
        J = ForwardDiff.jacobian(x_ -> _as_penalty_vec(model.resolver.fn(x_)[idx[i]]), x)
        # Zero where the constraint is satisfied, matching the residual's own flat region.
        active = [!iszero(_violation(c.op, rj)) for rj in r]
        (active ./ c.tol) .* J
    end
    return reduce(vcat, blocks)
end

"""
    augment_constraints!(augmented, constraints) -> Vector{String}

Add each constraint's difference expression to a model dict under a reserved key, and return
those keys.

The resolver is the only thing that knows how to evaluate an expression over model parameters,
so a constraint is evaluated by making it a parameter of the model being compiled. The keys are
named so they cannot collide with a real parameter, and the caller drops them from any model it
reports. This is the mechanism `priors` already uses.
"""
function augment_constraints!(augmented::AbstractDict, constraints)
    keys_ = String[]
    for (i, c) in enumerate(constraints)
        k = "__constraint_$(i)__"
        augmented[k] = constraint_expr(c)
        push!(keys_, k)
    end
    return keys_
end

"""
    check_constraints(constraints, model_dict; verb = true) -> Vector{Bool}

Report which constraints the model satisfies as it stands, without fitting anything.

Useful before a long run: a starting model that violates a constraint is not an error — the
penalty will push it back — but it is worth knowing about, because it may equally mean the
constraint was written backwards.
"""
function check_constraints(constraints, model_dict::AbstractDict; verb::Bool = true)
    cs = parse_constraints(constraints)
    isempty(cs) && return Bool[]
    augmented = Dict{String,Any}(String(k) => v for (k, v) in model_dict)
    ckeys = augment_constraints!(augmented, cs)
    ok = Bool[]
    for (c, k) in zip(cs, ckeys)
        d = _resolve_numeric(k, augmented)
        v = _violation(c.op, first(_as_penalty_vec(d)) / c.tol)
        push!(ok, iszero(v))
        if verb && !iszero(v)
            @warn "Constraint is violated by the starting model" constraint = c distance = d
        end
    end
    return ok
end


# ─────────────────────────────────────────────────────────────────────────────
# Hard constraints, for the fitter that can take them
# ─────────────────────────────────────────────────────────────────────────────
#
# NLopt implements nonlinear constraints properly, and a properly-implemented constraint is
# worth far more than a penalty here. Measured on `demos/data/AlphaCenA.oifits`,
# capping a uniform-disc diameter 0.5 mas below its best fit costs Δχ² ≈ 1.5e6 — so a penalty
# stiff enough to bind would have to out-weigh that, and one that does not simply loses. The
# fit settles wherever the two gradients balance, which is neither the constrained optimum nor
# the unconstrained one, and nothing in the result says so.
#
# NLopt's constraint machinery has no such failure mode: `star,ud < cap` lands on `cap`.
#
# Not every algorithm takes constraints, and the package default (`:LD_LBFGS`) does not. Rather
# than overriding a caller's choice of algorithm — derivative-free versus gradient is a
# deliberate decision — such a method is wrapped in `:AUGLAG`, NLopt's augmented-Lagrangian
# driver, which solves a sequence of unconstrained problems with the chosen algorithm inside.

"Algorithms taking nonlinear inequality constraints without a wrapper."
const _NLOPT_INEQUALITY = (:LD_MMA, :LD_CCSAQ, :LD_SLSQP, :LN_COBYLA, :GN_ISRES, :GN_AGS)

"Algorithms also taking nonlinear equality constraints — a much shorter list."
const _NLOPT_EQUALITY = (:LD_SLSQP, :GN_ISRES)

"""
    nlopt_handles_constraints(method, constraints) -> Bool

Whether this NLopt algorithm can take these constraints directly, rather than through `AUGLAG`.
"""
function nlopt_handles_constraints(method::Symbol, constraints)
    method in _NLOPT_INEQUALITY || return false
    any(c -> c.op === :(=), constraints) ? method in _NLOPT_EQUALITY : true
end

"""
    constrained_opt(method, n, constraints; ftol_rel, xtol_rel, maxeval) -> (Opt, wrapped)

An NLopt `Opt` able to carry `constraints`, and whether it needed the `AUGLAG` wrapper.

The tolerances are set on the returned optimiser and, when wrapping, on the inner one too:
`AUGLAG` has no stopping criterion of its own to hand down, so an inner optimiser left at
NLopt's defaults would never terminate.
"""
function constrained_opt(method::Symbol, n::Integer, constraints;
                         ftol_rel = 1e-8, xtol_rel = 1e-6, maxeval = 2000)
    if isempty(constraints) || nlopt_handles_constraints(method, constraints)
        opt = Opt(method, n)
        NLopt.ftol_rel!(opt, ftol_rel); NLopt.xtol_rel!(opt, xtol_rel); NLopt.maxeval!(opt, maxeval)
        return opt, false
    end
    inner = Opt(method, n)
    NLopt.ftol_rel!(inner, ftol_rel); NLopt.xtol_rel!(inner, xtol_rel); NLopt.maxeval!(inner, maxeval)
    opt = Opt(:AUGLAG, n)
    NLopt.local_optimizer!(opt, inner)          # NLopt copies it, so `inner` need not outlive this
    NLopt.ftol_rel!(opt, ftol_rel); NLopt.xtol_rel!(opt, xtol_rel); NLopt.maxeval!(opt, maxeval)
    return opt, true
end

"""
    add_nlopt_constraints!(opt, constraints, model, idx, x0)

Register each constraint on an NLopt optimiser, in NLopt's own `c(x) ≤ 0` / `c(x) = 0` form.

`>` and `>=` are registered negated, since NLopt knows only one direction. Each constraint's
`tol` becomes NLopt's feasibility tolerance for it, which is the same quantity it means in a
penalty — the violation that counts as none.

A chromatic parameter makes the resolver return a vector, and a constraint on it holds at every
wavelength; each element is registered separately. `x0` is evaluated once to learn how many
there are.
"""
function add_nlopt_constraints!(opt, constraints, model, idx, x0)
    for (i, c) in enumerate(constraints)
        s = (c.op === :> || c.op === :≥) ? -1.0 : 1.0
        n_out = length(_as_penalty_vec(model.resolver.fn(x0)[idx[i]]))
        for j in 1:n_out
            f = let s = s, ii = idx[i], jj = j, model = model
                function (x::Vector, g::Vector)
                    if length(g) > 0
                        J = ForwardDiff.jacobian(
                                x_ -> _as_penalty_vec(model.resolver.fn(x_)[ii]), x)
                        g .= s .* view(J, jj, :)
                    end
                    return s * _as_penalty_vec(model.resolver.fn(x)[ii])[jj]
                end
            end
            c.op === :(=) ? NLopt.equality_constraint!(opt, f, c.tol) :
                            NLopt.inequality_constraint!(opt, f, c.tol)
        end
    end
    return opt
end

"""
    constraint_violations(x, constraints, model, idx) -> Vector{Float64}

How far each constraint is from being satisfied at `x`, in its own `tol` units, and zero where
it is satisfied.

The check worth running on a converged fit: a hard constraint that the optimiser could not
satisfy is a fact about the fit, and it should be reported rather than assumed away.
"""
function constraint_violations(x, constraints, model, idx)
    isempty(constraints) && return Float64[]
    vals = model.resolver.fn(x)
    return [maximum(abs, _violation.(c.op, _as_penalty_vec(vals[idx[i]]) ./ c.tol))
            for (i, c) in enumerate(constraints)]
end
