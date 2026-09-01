# The Model perspective's data layer.
#
# Assertions are against the OITOOLS parser rather than against expectations written here: the
# panel's whole purpose is to report what the parser did, so a test that encoded its own idea
# of the answer could agree with itself while disagreeing with the model being fitted.

@testset "Model data layer" begin

    md = Dict{String,Any}(
        "star,ud"   => 6.5,
        "star,f"    => 0.8,
        "disk,fwhm" => "\$star,ud * 3",     # derived
        "disk,f"    => 0.2,
        "disk,pa"   => 45.0,
        "PA"        => 30.0,                # a global: bare key, no component
    )
    free = ["star,ud", "disk,pa"]

    @testset "every key becomes exactly one row" begin
        rows = model_rows(md, free)
        @test length(rows) == length(md)
        @test Set(r.key for r in rows) == Set(keys(md))
    end

    @testset "the three states are exclusive and correctly assigned" begin
        rows = model_rows(md, free)
        by = Dict(r.key => r for r in rows)
        @test by["star,ud"].mode   === PARAM_FREE
        @test by["disk,pa"].mode   === PARAM_FREE
        @test by["disk,fwhm"].mode === PARAM_EXPR
        @test by["star,f"].mode    === PARAM_FIXED
        @test by["PA"].mode        === PARAM_FIXED
        # free and derived cannot coexist: the resolver raises on a string in list_free_params,
        # so no row may be both
        @test !any(r -> r.mode === PARAM_FREE && !isempty(r.expr), rows)
    end

    @testset "fit indices match list_free_params positions" begin
        rows = model_rows(md, free)
        v = free_parameter_vector(rows)
        @test [r.key for r in v] == free          # row i IS x[i]
        @test [r.fitindex for r in v] == 1:length(free)
        @test all(r.fitindex == 0 for r in rows if r.mode !== PARAM_FREE)
    end

    @testset "derived values are resolved, not left blank" begin
        by = Dict(r.key => r for r in model_rows(md, free))
        @test by["disk,fwhm"].value ≈ 6.5 * 3     # $star,ud * 3
        @test by["disk,fwhm"].expr == "\$star,ud * 3"
    end

    @testset "bounds come from the parser's own defaults" begin
        rows = model_rows(md, free)
        by   = Dict(r.key => r for r in rows)
        lb, ub = default_bounds(md, free)
        @test by["star,ud"].lb == lb["star,ud"]
        @test by["star,ud"].ub == ub["star,ud"]
        @test by["disk,pa"].lb == lb["disk,pa"]
        # supplied bounds win over the defaults
        rows2 = model_rows(md, free; lb = Dict("star,ud" => 5.0), ub = Dict("star,ud" => 8.0))
        b2 = Dict(r.key => r for r in rows2)["star,ud"]
        @test (b2.lb, b2.ub) == (5.0, 8.0)
    end

    @testset "at-bound flags a pinned parameter, not a small one" begin
        # default_bounds gives a diameter 0..1000, so a span-relative test would call every
        # small diameter pinned. What matters is a value sitting ON its bound.
        rows = model_rows(Dict{String,Any}("a,ud" => 6.5, "b,ud" => 0.0, "b,f" => 1.0),
                          ["a,ud", "b,ud"])
        by = Dict(r.key => r for r in rows)
        @test by["a,ud"].atbound == false
        @test by["b,ud"].atbound == true
    end

    @testset "ordering: globals first, then geometry key before its modifiers" begin
        rows = model_rows(md, free)
        @test first(rows).component == GLOBAL_COMPONENT
        pos  = Dict(r.key => i for (i, r) in enumerate(rows))
        @test pos["star,ud"]   < pos["star,f"]      # geometry key leads its component
        @test pos["disk,fwhm"] < pos["disk,pa"]
    end

    @testset "inspection reports what the parser understood" begin
        insp = model_inspection(md)
        byname = Dict(c.name => c for c in insp.components)
        @test Set(keys(byname)) == Set(OITOOLS._component_names(md))
        @test byname["star"].kind         === OITOOLS._identify_kind("star", md)
        @test byname["star"].geometry_key == "ud"     # the key that decided the kind
        @test byname["disk"].geometry_key == "fwhm"
        @test insp.globals == ["PA"]
    end

    @testset "the silent misparse is made visible" begin
        # `gaussian_ring` is not a geometry key. With `fwhm` also present the component is a
        # plain Gaussian and `gaussian_ring` an inert parameter — no error, no warning at fit
        # time. Two shipped demos had exactly this.
        bad = Dict{String,Any}("ring,gaussian_ring" => 3.0, "ring,fwhm" => 5.0, "ring,f" => 1.0)
        insp = model_inspection(bad)
        c = only(insp.components)
        @test c.kind === :gaussian                    # NOT a ring
        @test c.geometry_key == "fwhm"                # and this is why
        @test "ring,gaussian_ring" in c.unrecognised  # named, in the component that owns it
        @test "ring,gaussian_ring" in insp.unrecognised
    end

    @testset "a clean model reports nothing unrecognised" begin
        @test isempty(model_inspection(md).unrecognised)
    end

    @testset "broadcasting is reported, because it changes every parameter" begin
        # One chromatic expression makes the resolver broadcast ALL of them, so every
        # parameter silently becomes a per-uv-point vector.
        @test model_inspection(md).broadcasting == false
        chromatic = merge(md, Dict{String,Any}("disk,f" => "0.2 * (\$WL / 1.6)^2"))
        @test model_inspection(chromatic).broadcasting == true
    end

    @testset "chromatic rows resolve at a wavelength, or not at all" begin
        # A chromatic expression has no single value: the resolver broadcasts it to one entry
        # per uv point. Without a wavelength the row cannot be resolved and must say so with a
        # NaN -- the failure to avoid is reporting a number that is not the parameter.
        chromatic = Dict{String,Any}("s,ud" => 2.0, "s,f" => "0.7 * (\$WL/1.6e-6)^-4")
        r0 = only(filter(r -> r.key == "s,f", model_rows(chromatic, String[])))
        @test r0.mode === PARAM_EXPR
        @test isnan(r0.value)

        r1 = only(filter(r -> r.key == "s,f", model_rows(chromatic, String[]; wl = 1.5e-6)))
        @test r1.value ≈ 0.7 * (1.5 / 1.6)^-4
        r2 = only(filter(r -> r.key == "s,f", model_rows(chromatic, String[]; wl = 1.6e-6)))
        @test r2.value ≈ 0.7

        # The substitution must not leak into the model: WL is a display value, and a model
        # that gained a global called WL would be a different model.
        @test !haskey(chromatic, "WL")
        @test length(model_rows(chromatic, String[]; wl = 1.6e-6)) == length(chromatic)
    end

    @testset "malformed input is described, not thrown at the caller" begin
        # A string named in list_free_params is malformed; the row must say so rather than
        # pretend a value. The GUI's job is to prevent this state, and to be legible when a
        # model arrives in it from a file.
        rows = model_rows(Dict{String,Any}("s,ud" => "\$other", "s,f" => 1.0), ["s,ud"])
        r = only(filter(r -> r.key == "s,ud", rows))
        @test r.mode === PARAM_FREE
        @test isnan(r.value)
    end
end
