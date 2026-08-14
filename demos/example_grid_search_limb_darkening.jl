##############################################################################
#  example_grid_search_limb_darkening.jl
#
#  Diameter vs linear limb-darkening coefficient χ² map for α Centauri A
#  (VLTI/PIONIER, Kervella et al. 2017, A&A 597, A137 — `aa29505-16.pdf`).
#
#  The point of the map is the *shape* of the valley: θ and u are strongly and
#  positively correlated, which is why a limb-darkened diameter is always less
#  well determined than a uniform-disc one from the same data.
#
#  Run headless with:  MPLBACKEND=Agg julia --project=.. <this file>
##############################################################################

using OITOOLS, PythonPlot, Printf

# The OIFITS carries the raw PIONIER wavelength scale; Kervella et al. multiply it by
# γ = 1.00481 (their Eq. 1) before quoting diameters, so multiply a fitted diameter by γ
# to compare with their Table 3. See example_model_fitting_limb_darkening.jl.
const GAMMA = 1.00481

data = readoifits(joinpath(@__DIR__, "data", "AlphaCenA.oifits"); warn=false, verbose=false)[1,1]

diameters = range(7.5, 9.5, step=0.005)   # mas
ld        = range(0.0, 0.5, step=0.002)   # linear limb-darkening coefficient u

model = Dict{String,Any}(
    "star,ldlin" => 8.0,
    "star,u"     => 0.3,
    "star,f"     => 1.0,
)
free_params = ["star,ldlin", "star,u"]
params = dict_to_model(model, free_params)

chi2 = [model_to_chi2(params, [d, u], data; weights=[1.0,0,0]) for d in diameters, u in ld]

res = findmin(chi2)
d_best, u_best = diameters[res[2][1]], ld[res[2][2]]
@printf("Grid minimum:  χ² = %.2f  (χ²r = %.3f)   θ = %.4f mas (×γ = %.4f)   u = %.4f\n",
        res[1], res[1]/data.nv2, d_best, GAMMA*d_best, u_best)
@printf("Kervella+2017 Table 3, linear (u free):  θ = 8.458 ± 0.005 ± 0.035 mas,  u = 0.1761 ± 0.0062\n")

# ── Plot ─────────────────────────────────────────────────────────────────────
# χ² spans three orders of magnitude across this box, so the raw values would render as
# a single bright pixel on black. Plot Δχ² relative to the minimum on a log scale.
#
# Two panels, because the two things worth seeing live at different scales: the
# degeneracy valley is only visible across the whole box, and the confidence region is
# ~0.02 mas wide, i.e. a few pixels there.
dchi2 = chi2 .- res[1]
u_pub, d_pub = 0.1761, 8.458/GAMMA          # published values on the OIFITS' own scale

fig, axs = subplots(1, 2, figsize=(12, 5))

for (ax, (dlim, ulim, ttl)) in zip(axs,
        ((extrema(diameters), extrema(ld), "full grid — the θ/u degeneracy"),
         ((8.36, 8.48), (0.12, 0.24),      "zoom — 1σ, 2σ, 3σ joint confidence")))
    im = ax.imshow(log10.(dchi2 .+ 1), cmap=ColorMap("gist_heat_r"), interpolation="none",
                   aspect="auto", origin="lower",
                   extent=[first(ld), last(ld), first(diameters), last(diameters)])
    # Δχ² = 2.30, 6.17, 11.8 are the 68.3 / 95.4 / 99.7% levels for two free parameters
    ax.contour(collect(ld), collect(diameters), dchi2,
               levels=[2.30, 6.17, 11.8], colors="black", linewidths=0.8)
    ax.plot([u_best], [d_best], "+", color="black",     markersize=11, label="grid minimum")
    ax.plot([u_pub],  [d_pub],  "x", color="tab:blue",  markersize=11, label="Kervella+2017 (θ/γ)")
    ax.set_xlim(ulim...); ax.set_ylim(dlim...)
    ax.set_xlabel("Linear limb-darkening coefficient  u")
    ax.set_ylabel("Diameter  (mas)")
    ax.set_title(ttl, fontsize=10)
    ax.legend(loc="lower right", fontsize=8)
    fig.colorbar(im, ax=ax, aspect=40, label=L"$\log_{10}(\Delta\chi^2 + 1)$")
end

suptitle(L"$\alpha$ Cen A — $\chi^2$ map, linear limb-darkened disc")
tight_layout()

savefig(joinpath(@__DIR__, "alphacenA_ldlin_chi2map.png"), dpi=120)
println("Wrote alphacenA_ldlin_chi2map.png")
