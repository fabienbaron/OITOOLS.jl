# Constraints and TOML model files.
#
# The assertion that matters is §"a constraint the data fights is still obeyed". A soft penalty
# passes an easy test — put the constraint somewhere the fit was going anyway and it looks
# enforced — and fails the real one. On this dataset, capping the diameter 0.5 mas below its
# best fit costs Δχ² ≈ 1.5e6, so anything short of a real constraint settles part-way and
# reports a number that is neither the constrained nor the unconstrained answer.
#
# Nothing here hard-codes a fitted value. The unconstrained fit is run first and the constraint
# is placed relative to its result, so the test states a relationship between two fits rather
# than a number that a change of data or optimiser would invalidate.

using OITOOLS, Test, LinearAlgebra

@testset "constraints and model files" begin

data = readoifits(joinpath(@__DIR__, "..", "demos", "data", "AlphaCenA.oifits");
                  warn = false, verbose = false)[1, 1]

UD    = Dict{String,Any}("star,ud" => 8.0, "star,f" => 1.0)
UDLB  = Dict("star,ud" => 0.0)
UDUB  = Dict("star,ud" => 20.0)

@testset "operators are accepted as written, and only real relations" begin
    for (s, sym) in ("<" => :<, "<=" => :≤, "≤" => :≤,
                     ">" => :>, ">=" => :≥, "≥" => :≥,
                     "=" => :(=), "==" => :(=))
        @test ModelConstraint("a,ud", s, 1.0).op === sym
        @test ModelConstraint("a,ud", Symbol(s), 1.0).op === sym   # symbols too
    end
    @test_throws ErrorException ModelConstraint("a,ud", "+", 1.0)
    @test_throws ErrorException ModelConstraint("a,ud", "!=", 1.0)
    # tol is a scale, so it cannot be zero or negative
    @test_throws ErrorException ModelConstraint("a,ud", "<", 1.0; tol = 0.0)
    @test_throws ErrorException ModelConstraint("a,ud", "<", 1.0; tol = -1.0)
end

@testset "each side is referenced the way a model expression would be" begin
    ce(c) = OITOOLS.constraint_expr(c)
    @test ce(ModelConstraint("disk,diamout", ">", "disk,diamin")) ==
          "(\$disk,diamout) - (\$disk,diamin)"
    @test ce(ModelConstraint("star,ud", "<", 4.0))    == "(\$star,ud) - (4.0)"
    @test ce(ModelConstraint("star,f + disk,f", "=", 1.0)) == "(\$star,f + \$disk,f) - (1.0)"
    # a bare name is a global, and referencing one takes the same `$`
    @test ce(ModelConstraint("PA", "<", 90.0))        == "(\$PA) - (90.0)"
    # a side the caller already wrote in the expression language is left exactly as written
    @test ce(ModelConstraint("\$star,ud * 2", ">", 1.0)) == "(\$star,ud * 2) - (1.0)"
end

@testset "constraints arrive in whichever shape the caller has them" begin
    # the PMOIRED `prior` tuple layout, so a converted model needs no reshaping
    cs = parse_constraints([("a,ud", "<", 3.0), ("a,ud", ">", "b,ud", 1e-4)])
    @test length(cs) == 2
    @test cs[1].rhs === 3.0 && cs[1].tol == DEFAULT_CONSTRAINT_TOL
    @test cs[2].rhs == "b,ud" && cs[2].tol == 1e-4
    # the shape a TOML [[constraints]] entry parses to
    d = parse_constraints([Dict("param" => "a,ud", "op" => ">", "value" => 0.5, "tol" => 0.01)])
    @test only(d) == ModelConstraint("a,ud", ">", 0.5; tol = 0.01)
    @test only(parse_constraints([(param = "a,ud", op = ">", value = 0.5)])).lhs == "a,ud"
    @test isempty(parse_constraints([]))
    @test_throws ErrorException parse_constraints([("a,ud", "<")])
end

@testset "the penalty is one-sided, and its gradient matches finite differences" begin
    md   = Dict{String,Any}("star,ud" => 8.0, "star,f" => 1.0)
    cons = parse_constraints([ModelConstraint("star,ud", "<", 7.0)])
    aug  = copy(md)
    keys_ = OITOOLS.augment_constraints!(aug, cons)
    fm   = parse_model(aug, ["star,ud"]; nB_workspace = 1)
    idx  = [findfirst(==(k), fm.all_names) for k in keys_]

    # satisfied costs exactly zero — a constraint the fit already obeys must not pull it
    @test OITOOLS.constraint_penalty([6.5], cons, fm, idx) == 0.0
    g = zeros(1)
    @test OITOOLS.constraint_penalty_and_grad!(g, [6.5], cons, fm, idx) == 0.0
    @test g == [0.0]

    # violated costs (violation/tol)^2
    @test OITOOLS.constraint_penalty([7.5], cons, fm, idx) ≈ (0.5 / 1e-3)^2

    for x in ([8.0], [7.5], [7.2])
        g = zeros(1)
        OITOOLS.constraint_penalty_and_grad!(g, x, cons, fm, idx)
        h  = 1e-7
        fd = (OITOOLS.constraint_penalty(x .+ h, cons, fm, idx) -
              OITOOLS.constraint_penalty(x .- h, cons, fm, idx)) / (2h)
        @test g[1] ≈ fd rtol = 1e-5
    end

    # the least-squares form squares to the same penalty, and its Jacobian matches
    r = OITOOLS.constraint_residuals([7.5], cons, fm, idx)
    @test sum(abs2, r) ≈ OITOOLS.constraint_penalty([7.5], cons, fm, idx)
    J = OITOOLS.constraint_jacobian([7.5], cons, fm, idx)
    @test size(J) == (length(r), 1)
    g75 = zeros(1)
    OITOOLS.constraint_penalty_and_grad!(g75, [7.5], cons, fm, idx)
    @test 2 * (J' * r)[1] ≈ g75[1] rtol = 1e-8
end

@testset "check_constraints reports the starting model without fitting" begin
    ok = check_constraints([ModelConstraint("r,fwhmout", ">", "r,fwhmin")],
                           Dict{String,Any}("r,fwhmin" => 2.0, "r,fwhmout" => 5.0, "r,f" => 1.0))
    @test ok == [true]
    bad = @test_logs (:warn,) match_mode = :any check_constraints(
              [ModelConstraint("r,fwhmout", ">", "r,fwhmin")],
              Dict{String,Any}("r,fwhmin" => 5.0, "r,fwhmout" => 2.0, "r,f" => 1.0))
    @test bad == [false]
    @test isempty(check_constraints([], UD))
end

@testset "a constraint the data fights is still obeyed" begin
    # The point of the whole file. Δχ² across this cap is ~1.5e6, which any soft penalty of a
    # sane stiffness loses to.
    free = ["star,ud"]
    r0   = fit_model(UD, free, data; lb = UDLB, ub = UDUB)
    cap  = r0.x_opt[1] - 0.5

    for method in (:LD_LBFGS,        # the default: takes no constraints, so AUGLAG wraps it
                   :LN_NELDERMEAD,   # derivative-free, also wrapped
                   :LD_MMA,          # takes them natively
                   :LN_COBYLA)       # takes them natively, derivative-free
        r = fit_model(UD, free, data; lb = UDLB, ub = UDUB, method,
                      constraints = [ModelConstraint("star,ud", "<", cap)])
        @testset "$method" begin
            @test r.x_opt[1] <= cap + 1e-3       # inside, to the constraint's own tolerance
            @test r.x_opt[1] ≈ cap atol = 1e-3   # and pressed against it, not part-way
        end
    end
end

@testset "a constraint that is not binding changes nothing" begin
    free = ["star,ud"]
    r0 = fit_model(UD, free, data; lb = UDLB, ub = UDUB)
    r1 = fit_model(UD, free, data; lb = UDLB, ub = UDUB,
                   constraints = [ModelConstraint("star,ud", "<", r0.x_opt[1] + 5)])
    @test r1.x_opt ≈ r0.x_opt
    # and the reported chi2 is the data chi2, not a penalised one
    @test r1.chi2 ≈ r0.chi2
end

@testset "a relation between two parameters holds at the optimum" begin
    # A Gaussian ring whose outer FWHM starts smaller than its inner one: the model is
    # unphysical as given, and the constraint is what makes the answer meaningful.
    md   = Dict{String,Any}("r,fwhmin" => 3.0, "r,fwhmout" => 2.0, "r,f" => 1.0)
    free = ["r,fwhmin", "r,fwhmout"]
    r = fit_model(md, free, data;
                  lb = Dict("r,fwhmin" => 0.1, "r,fwhmout" => 0.1),
                  ub = Dict("r,fwhmin" => 20.0, "r,fwhmout" => 20.0),
                  constraints = [ModelConstraint("r,fwhmout", ">", "r,fwhmin")])
    @test r.x_opt[2] > r.x_opt[1]
end

@testset "an equality constraint ties parameters together" begin
    md   = Dict{String,Any}("a,ud" => 3.0, "a,f" => 0.5, "b,ud" => 8.0, "b,f" => 0.3)
    free = ["a,ud", "a,f", "b,ud", "b,f"]
    r = fit_model(md, free, data;
                  lb = Dict(p => 0.0 for p in free),
                  ub = Dict("a,ud" => 20.0, "a,f" => 1.0, "b,ud" => 20.0, "b,f" => 1.0),
                  constraints = [ModelConstraint("a,f + b,f", "=", 1.0; tol = 1e-6)])
    @test r.x_opt[2] + r.x_opt[4] ≈ 1.0 atol = 1e-4
end

@testset "every fitter takes constraints, and says so when it cannot" begin
    free = ["star,ud"]
    c = [ModelConstraint("star,ud", "<", 7.0)]
    # Levenberg-Marquardt: the constraint enters as residual rows
    rl = fit_model_lsqfit(UD, free, data; lb = UDLB, ub = UDUB, constraints = c)
    @test rl.x_opt[1] < fit_model_lsqfit(UD, free, data; lb = UDLB, ub = UDUB).x_opt[1]
    # and the reported chi2 counts data rows only, so it stays comparable
    @test rl.chi2 ≈ OITOOLS.chi2_flat(rl.model, rl.x_opt, data;
                              weights = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]) rtol = 1e-8

    # a compiled FlatModel cannot grow the constraint expressions, and must say so rather than
    # accept the argument and drop it
    fm = parse_model(UD, free; nB_workspace = size(data.uv, 2))
    @test_throws ErrorException fit_model(fm, [8.0], data; constraints = c)
    @test_throws ErrorException fit_model_lsqfit(fm, [8.0], data; constraints = c)
end

@testset "TOML model files round-trip everything a fit needs" begin
    path = joinpath(mktempdir(), "binary.toml")
    md = Dict{String,Any}("star,ud"    => 6.5,
                          "star,f"     => 0.8,
                          "disk,fwhm"  => "\$star,ud * 3",
                          "disk,f"     => 0.2,
                          "disk,pa"    => 45.0,
                          "PA"         => 30.0)
    free = ["star,ud", "disk,pa"]
    lb, ub = default_bounds(md, free)
    cons = [ModelConstraint("disk,fwhm", ">", "star,ud"; tol = 1e-3),
            ModelConstraint("star,f + disk,f", "=", 1.0)]
    priors = [("star,ud", 6.0, 0.5)]

    write_model_file(path, md; free, lb, ub, constraints = cons, priors, name = "test binary")
    m = read_model_file(path)

    @test m.model == md
    @test m.free  == free
    @test m.name  == "test binary"
    @test m.constraints == cons
    @test m.priors == priors
    for p in free
        @test m.lb[p] == lb[p]
        @test m.ub[p] == ub[p]
    end

    # what read gives back is what a fitter takes
    r = fit_model(m.model, m.free, data; m.lb, m.ub, m.constraints, m.priors)
    @test length(r.x_opt) == length(free)

    # writing is deterministic, so two files of the same model compare equal
    p2 = joinpath(mktempdir(), "binary.toml")
    write_model_file(p2, md; free, lb, ub, constraints = cons, priors, name = "test binary")
    @test read(path, String) == read(p2, String)
end

@testset "bounds that constrain nothing are not written" begin
    path = joinpath(mktempdir(), "m.toml")
    md = Dict{String,Any}("star,ud" => 6.5, "star,f" => 1.0)
    write_model_file(path, md; free = ["star,ud"],
                     lb = Dict("star,ud" => -Inf), ub = Dict("star,ud" => Inf))
    @test !occursin("[bounds]", read(path, String))
    m = read_model_file(path)
    @test isempty(m.lb) && isempty(m.ub)
    # a half-open bound is a real bound and survives
    write_model_file(path, md; free = ["star,ud"],
                     lb = Dict("star,ud" => 0.0), ub = Dict("star,ud" => Inf))
    m = read_model_file(path)
    @test m.lb["star,ud"] == 0.0 && m.ub["star,ud"] == Inf
end

@testset "a malformed model file is described, not silently accepted" begin
    dir = mktempdir()
    write(joinpath(dir, "nomodel.toml"), "free = []\n")
    @test_throws ErrorException read_model_file(joinpath(dir, "nomodel.toml"))
    @test_throws ErrorException read_model_file(joinpath(dir, "absent.toml"))

    # a free parameter must exist and be numeric: the resolver raises on a string there, and
    # the file is the last place that can still say which key is wrong
    write(joinpath(dir, "ghost.toml"), """
          free = ["star,ud"]
          [model]
          "star,f" = 1.0
          """)
    @test_throws ErrorException read_model_file(joinpath(dir, "ghost.toml"))
    write(joinpath(dir, "expr.toml"), """
          free = ["star,ud"]
          [model]
          "star,ud" = "\$other,ud"
          """)
    @test_throws ErrorException read_model_file(joinpath(dir, "expr.toml"))
    @test_throws ErrorException write_model_file(joinpath(dir, "x.toml"),
                                    Dict{String,Any}("a,ud" => "\$b,ud"); free = ["a,ud"])

    write(joinpath(dir, "bounds.toml"), """
          [model]
          "star,ud" = 1.0
          [bounds]
          "star,ud" = [5.0]
          """)
    @test_throws ErrorException read_model_file(joinpath(dir, "bounds.toml"))
    write(joinpath(dir, "empty.toml"), """
          [model]
          "star,ud" = 1.0
          [bounds]
          "star,ud" = [5.0, 1.0]
          """)
    @test_throws ErrorException read_model_file(joinpath(dir, "empty.toml"))
end

@testset "a hand-written file reads without the writer having produced it" begin
    dir = mktempdir(); path = joinpath(dir, "hand.toml")
    write(path, """
          free = ["ring,fwhmin", "ring,fwhmout"]

          [model]
          "ring,fwhmin"  = 1.5
          "ring,fwhmout" = 2.5
          "ring,f"       = 1.0

          [bounds]
          "ring,fwhmin"  = [0.1, 10.0]
          "ring,fwhmout" = [0.2, 20.0]

          [[constraints]]
          param = "ring,fwhmout"
          op    = ">"
          value = "ring,fwhmin"
          tol   = 0.001
          """)
    m = read_model_file(path)
    @test m.free == ["ring,fwhmin", "ring,fwhmout"]
    @test m.model["ring,fwhmin"] == 1.5
    @test m.lb["ring,fwhmout"] == 0.2
    @test only(m.constraints) == ModelConstraint("ring,fwhmout", ">", "ring,fwhmin")
    @test isempty(m.priors)
end

end
