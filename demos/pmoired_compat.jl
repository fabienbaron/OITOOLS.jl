"""
pmoired_compat.jl

Utility to convert a PMOIRED-style Python dict literal (as a Julia String)
into an equivalent Julia Dict literal string.

Transformation rules applied:
  { ... }           ->  Dict( ... )
  'key': value      ->  "key" => value         (outside string literals)
  :                 ->  =>                      (outside string literals, single colon only)
  'plain string'    ->  "plain string"
  'expr with \$ref' ->  raw"expr with \$ref"   (raw string prevents interpolation)
  True / False      ->  true / false
  None              ->  nothing

Usage:
    julia_str = pmoired_to_julia(python_str)
    param     = eval(Meta.parse(julia_str))   # only at parse_model time, not in hot path
"""

"""
    pmoired_to_julia(s) -> String

Convert a PMOIRED Python dict literal string to an equivalent Julia Dict
literal string.  Handles nested dicts, expression strings with `\$` references,
single-quoted strings, Python boolean/None literals, and list literals.

Example
-------
```python
# Python / PMOIRED:
param = {
    'star,ud':    3.2,
    'star,f':     0.7,
    'ring,udout': '\$star,ud * 8',
    'ring,f':     '1 - \$star,f',
    'ring,incl':  30.0,
}
fitOnly = ['star,ud', 'star,f', 'ring,udout', 'ring,f']
```

```julia
# Generated Julia:
param = Dict(
    "star,ud"    => 3.2,
    "star,f"     => 0.7,
    "ring,udout" => raw"\$star,ud * 8",
    "ring,f"     => raw"1 - \$star,f",
    "ring,incl"  => 30.0,
)
fitOnly = ["star,ud", "star,f", "ring,udout", "ring,f"]
```
"""
function pmoired_to_julia(s::AbstractString)::String
    buf   = IOBuffer()
    i     = firstindex(s)
    last  = lastindex(s)

    # ------------------------------------------------------------------ helpers
    # Advance iterator safely
    nextidx(s, i) = i > last ? i : nextind(s, i)

    # Read a quoted string (single or double), return (content, next_i)
    function read_string(s, start, delim)
        content = IOBuffer()
        i = nextidx(s, start)          # skip opening quote
        while i <= last
            c = s[i]
            if c == '\\'               # escape sequence — preserve verbatim
                write(content, c)
                i = nextidx(s, i)
                if i <= last
                    write(content, s[i])
                    i = nextidx(s, i)
                end
            elseif c == delim          # closing quote
                i = nextidx(s, i)
                break
            else
                write(content, c)
                i = nextidx(s, i)
            end
        end
        return String(take!(content)), i
    end

    # ------------------------------------------------------------------ main loop
    while i <= last
        c = s[i]

        # ---- single-quoted string ----------------------------------------
        if c == '\''
            content, i = read_string(s, i, '\'')
            if occursin('$', content)
                write(buf, "raw\"", content, "\"")
            else
                write(buf, '"', content, '"')
            end

        # ---- double-quoted string ----------------------------------------
        elseif c == '"'
            content, i = read_string(s, i, '"')
            if occursin('$', content)
                write(buf, "raw\"", content, "\"")
            else
                write(buf, '"', content, '"')
            end

        # ---- Python hash comment — preserve as Julia comment -------------
        elseif c == '#'
            while i <= last && s[i] != '\n'
                write(buf, s[i])
                i = nextidx(s, i)
            end

        # ---- { -> Dict( --------------------------------------------------
        elseif c == '{'
            write(buf, "Dict(")
            i = nextidx(s, i)

        # ---- } -> ) ------------------------------------------------------
        elseif c == '}'
            write(buf, ')')
            i = nextidx(s, i)

        # ---- : -> =>  (key-value separator; avoid :: ) -------------------
        elseif c == ':'
            next_i = nextidx(s, i)
            if next_i <= last && s[next_i] == ':'   # :: — pass through
                write(buf, "::")
                i = nextidx(s, next_i)
            else
                write(buf, " =>")
                i = next_i
            end

        # ---- everything else passes through verbatim --------------------
        else
            write(buf, c)
            i = nextidx(s, i)
        end
    end

    result = String(take!(buf))

    # ---- Python keywords → Julia keywords --------------------------------
    # Use word-boundary anchors to avoid substring replacement
    result = replace(result, r"\bTrue\b"  => "true")
    result = replace(result, r"\bFalse\b" => "false")
    result = replace(result, r"\bNone\b"  => "nothing")

    return result
end


"""
    pmoired_to_julia_file(infile, outfile)

Read a Python/PMOIRED notebook snippet from `infile`, convert it, write to `outfile`.
Lines that look like `param = {...}` or `fitOnly = [...]` are converted;
other lines are passed through.
"""
function pmoired_to_julia_file(infile::AbstractString, outfile::AbstractString)
    open(outfile, "w") do out
        for line in eachline(infile)
            println(out, pmoired_to_julia(line))
        end
    end
end


# ===========================================================================
# Self-tests
# ===========================================================================
function _test_pmoired_to_julia()
    println("Running pmoired_to_julia tests...\n")
    pass = 0; fail = 0

    function check(input, expected, label)
        got = pmoired_to_julia(input)
        if got == expected
            println("  PASS  $label")
            pass += 1
        else
            println("  FAIL  $label")
            println("    input:    $input")
            println("    expected: $expected")
            println("    got:      $got")
            fail += 1
        end
    end

    # -- basic key-value ---------------------------------------------------
    check("{'ud': 3.2}",
          "Dict(\"ud\" => 3.2)",
          "simple ud key")

    # -- multi-component comma key -----------------------------------------
    check("{'star,ud': 0.2, 'star,f': 0.7}",
          "Dict(\"star,ud\" => 0.2, \"star,f\" => 0.7)",
          "comma in key name")

    # -- plain string value ------------------------------------------------
    check("{'type': 'ring'}",
          "Dict(\"type\" => \"ring\")",
          "plain string value")

    # -- expression with $ — must become raw"..." --------------------------
    check("{'ring,f': '1 - \$star,f'}",
          "Dict(\"ring,f\" => raw\"1 - \$star,f\")",
          "dollar expression -> raw string")

    # -- expression with multiple $ refs -----------------------------------
    check("{'ring,udout': '\$star,ud * \$scale'}",
          "Dict(\"ring,udout\" => raw\"\$star,ud * \$scale\")",
          "multiple dollar refs in expression")

    # -- Python True/False/None --------------------------------------------
    check("{'flag': True, 'other': False, 'val': None}",
          "Dict(\"flag\" => true, \"other\" => false, \"val\" => nothing)",
          "Python bool/None keywords")

    # -- list literal (fitOnly) --------------------------------------------
    check("['star,ud', 'ring,f']",
          "[\"star,ud\", \"ring,f\"]",
          "list of strings (fitOnly)")

    # -- nested dict -------------------------------------------------------
    check("{'a': {'b': 1}}",
          "Dict(\"a\" => Dict(\"b\" => 1))",
          "nested dict")

    # -- scientific notation and negative numbers --------------------------
    check("{'wl0': 1.6e-6, 'flux': -0.3}",
          "Dict(\"wl0\" => 1.6e-6, \"flux\" => -0.3)",
          "scientific notation and negatives")

    # -- escape sequence in string value -----------------------------------
    check("""{'note': 'it\\'s fine'}""",
          """Dict("note" => "it\\'s fine")""",
          "escaped single quote inside single-quoted string")

    # -- comment passthrough -----------------------------------------------
    check("# a comment",
          "# a comment",
          "comment passthrough")

    # -- whitespace preserved in keys, values  ----------------------------
    check("{ 'ud' : 3.2 }",
          "Dict( \"ud\"  => 3.2 )",
          "extra whitespace around : preserved")

    println("\n$pass passed, $fail failed.")
    return fail == 0
end


# Uncomment to run tests when this file is included:
# _test_pmoired_to_julia()
