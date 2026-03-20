# example_model_fitting_pmoired_conversion.jl
#
# Demonstrates converting PMOIRED Python dict literals to Julia using
# pmoired_to_julia() from OITOOLS.
#
# Usage:
#   include("example_model_fitting_pmoired_conversion.jl")

using OITOOLS

# ── Basic key-value ──────────────────────────────────────────────────────────
println("=== PMOIRED → Julia conversion examples ===\n")

input = "{'ud': 3.2}"
println("Input:  ", input)
println("Output: ", pmoired_to_julia(input))
println()

# ── Multi-component model with comma keys ────────────────────────────────────
input = "{'star,ud': 0.2, 'star,f': 0.7}"
println("Input:  ", input)
println("Output: ", pmoired_to_julia(input))
println()

# ── Expression strings with \$ references ────────────────────────────────────
input = "{'ring,udout': '\$star,ud * 8', 'ring,f': '1 - \$star,f'}"
println("Input:  ", input)
println("Output: ", pmoired_to_julia(input))
println()

# ── Python True/False/None ───────────────────────────────────────────────────
input = "{'flag': True, 'other': False, 'val': None}"
println("Input:  ", input)
println("Output: ", pmoired_to_julia(input))
println()

# ── fitOnly list ─────────────────────────────────────────────────────────────
input = "['star,ud', 'ring,f']"
println("Input:  ", input)
println("Output: ", pmoired_to_julia(input))
println()

# ── Nested dict ──────────────────────────────────────────────────────────────
input = "{'a': {'b': 1}}"
println("Input:  ", input)
println("Output: ", pmoired_to_julia(input))
println()

# ── Full PMOIRED model conversion ───────────────────────────────────────────
println("=== Full model example ===\n")

pmoired_model = """
{
    'star,ud':    3.2,
    'star,f':     0.7,
    'ring,udout': '\$star,ud * 8',
    'ring,f':     '1 - \$star,f',
    'ring,incl':  30.0,
}
"""

julia_str = pmoired_to_julia(pmoired_model)
println("PMOIRED input:")
println(pmoired_model)
println("Julia output:")
println(julia_str)

# Parse into an actual Julia Dict
model_dict = pmoired_to_dict(pmoired_model)
println("Parsed Dict:")
for (k, v) in sort(collect(model_dict))
    println("  \"$k\" => $v")
end
