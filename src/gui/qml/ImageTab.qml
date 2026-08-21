// OITOOLS — the Image perspective: turn one dataset into one reconstruction.
//
// The engine is the primary control, and it is not a cosmetic choice: the seven engines take
// their regularisation in five incompatible forms, three of them return a posterior rather
// than an image, and two of them are optional dependencies that may not be installed. So the
// panel below the engine selector CHANGES SHAPE with the engine instead of offering one
// flattened list of knobs that would be wrong for every engine but VMLMB.
//
// This file is the view only. Every point where Julia will be called is marked `// wire:`;
// nothing here computes, loads or reconstructs anything.

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml            // supplies `Julia`, through which the reconstruction is run
import Makie            // MakieArea, for this perspective's own canvas

Item {
    id: root

    // ── scaling, supplied by the window ───────────────────────────────────────
    //
    // Same contract as Main.qml: the window has already worked out the screen's scale, and a
    // tab that re-derived it would drift away from the chrome around it. Every pixel length
    // goes through dp() and every explicit point size through pt().
    // Colormap for the reconstruction. The list comes from Julia so the buttons and
    // `set_colormap!` cannot disagree about the names.
    property var colormapNames: []
    property string colormap: "viridis"

    // Names come from Julia so the buttons and `set_colormap!` cannot disagree.
    Component.onCompleted:
        root.colormapNames = Julia.shell_image_colormaps()
                                  .split("\n")
                                  .filter(function (t) { return t.length > 0 })

    property real uiScale:   1.0
    property real fontScale: 1.0
    property real baseFontPt: 11

    function dp(px)     { return Math.round(px * uiScale) }   // pixel lengths
    function pt(points) { return points * fontScale }         // explicit point sizes

    // ── image geometry (§5.1) ─────────────────────────────────────────────────
    property int  nx: 128
    property real pixsize: 0.2                                 // mas per pixel
    readonly property real fov: nx * pixsize                    // mas across the image
    property bool scaleRegs: true                               // scale_regularizers on pixsize change
    property string ftMode: "nfft"                              // setup_ft mode

    // Gaussian, not Dirac. A single lit pixel is the least committal image there is and the
    // worst possible start for an unregularised gradient descent: measured on 2004-data1 it
    // begins at chi2r 1.4e10 and ends at 1.2e6 with the flux blown up to 4234, where a
    // Gaussian start on the same data reaches chi2r 312 with the flux near 1.
    property string startImage: "gaussian"                      // dirac | gaussian | fits
    property string startImagePath: ""
    property real   startFwhm: -1                               // mas; negative ⇒ chosen from the FOV
    property int    startSeed: 1                                // for the random start

    // ── what the file actually contains (§1.8) ────────────────────────────────
    //
    // Julia pushes these from the valid-point counts after a load. They are the ONLY source of
    // "the file has none"; the two imaging exclusions below are a different thing entirely.
    property bool haveCvis:    false
    property bool haveV2:      false
    property bool haveT3amp:   false
    property bool haveT3phi:   false
    property bool haveDiffvis: false
    property bool haveFlux:    false

    // ── which observables the reconstruction uses ─────────────────────────────
    readonly property bool useCvis:    cvisBox.checked
    readonly property bool useV2:      v2Box.checked
    readonly property bool useT3amp:   t3ampBox.checked
    readonly property bool useT3phi:   t3phiBox.checked
    readonly property bool useDiffvis: diffvisBox.checked
    readonly property bool useFlux:    fluxBox.checked
    readonly property bool anyObservable: useV2 || useT3amp || useT3phi || useDiffvis

    // reconstruct()'s weights vector is three long, [V², T3amp, T3φ], and there is no fourth
    // slot to put a complex visibility or a flux in. That is a property of the engine and it
    // never changes with the file, so it gets its own wording: a user who sees the same grey
    // box for "your file has no T3amp" and "no imaging engine can use cvis" concludes the file
    // is broken.
    readonly property string engineExclusion:
        "not supported by any imaging engine — reconstruct() weights are [V², T3amp, T3φ] only.\n" +
        "This is permanent, and says nothing about what this file contains."

    // ── engine (§5.2) ─────────────────────────────────────────────────────────
    property string engine: "vmlmb"
    // Pigeons and OIVI are weak dependencies: OIVI is unregistered and needs a sibling VarInf
    // clone, so the entries stay visible-but-disabled with the reason on them rather than
    // vanishing, which would leave no trace of the capability existing at all.
    property bool pigeonsAvailable: false
    property bool oiviAvailable:    false

    // The three SQUEEZE rows share a sampler, so almost every branch wants "is this SQUEEZE"
    // rather than a list of keys. Naming it once is what keeps a new variant from having to be
    // added in six places.
    readonly property bool isSqueeze:
        engine === "squeeze" || engine === "tempering" || engine === "squeeze_sparco"

    // The SPARCO row IS the annealing sampler with a parametric component sampled beside the
    // image, so the panel follows the engine instead of a separate tick box.
    readonly property bool squeezeUseSparco: engine === "squeeze_sparco"

    readonly property bool stochastic: isSqueeze || engine === "vi"

    function engineReasonFor(key) {
        if (key === "tempering" && !pigeonsAvailable) return "needs Pigeons"
        if (key === "vi"        && !oiviAvailable)    return "needs OIVI"
        return ""
    }
    function engineAvailableFor(key) { return engineReasonFor(key) === "" }

    function engineIndexOf(key) {
        for (var i = 0; i < engineModel.count; ++i)
            if (engineModel.get(i).key === key) return i
        return 0
    }
    // Keep the combo on whatever Julia last set, so the engine can be driven from a script.
    onEngineChanged: engineBox.currentIndex = engineIndexOf(engine)

    // Positivity is a bound, not a term in the objective, and the reason it is unarguable
    // differs by engine — hiding it in the regulariser list would misrepresent both facts.
    // ── regularisation state, one block per engine shape (§5.3) ───────────────
    //
    // reconstruct() and reconstruct_hybrid() take ["name", μ, extra...] specs, so they share
    // one editable list; the other three take something that is not a spec list at all, and
    // pretending otherwise is exactly the flattening §5.3 warns against.
    property alias vmlmbRegs:   regModel          // rows of ["name", μ, extra...]
    property alias squeezeRegs: squeezeRegModel   // SQUEEZE's own six-name vocabulary

    // Only these names exist. Offering a free-text field would let a typo through, and an
    // unknown name is an error deep inside the engine rather than a message here.
    readonly property var specNames: [
        "centering", "tv", "tvsq", "l1l2", "l1l2w", "l1hyp", "l2sq",
        "compactness", "radialvar", "entropy", "support",
        "transspectral_tv", "transspectral_tvsq", "transspectral_structnorm",
        "transspectral_l1l2", "transspectral_grouptv", "transspectral_poly"
    ]
    // l1l2 cannot be evaluated without its α, so that row grows a field the others do not
    // have; the rest of these extras are optional and default inside the engine.
    function specExtraLabel(name) {
        switch (name) {
        case "l1l2":               return "α"
        case "transspectral_l1l2": return "δ"
        case "tv":                 return "ϵ"
        case "compactness":        return "w"
        case "transspectral_poly": return "degree"
        default:                   return ""
        }
    }
    function specExtraRequired(name) { return name === "l1l2" || name === "transspectral_l1l2" }
    function specNeedsImage(name)    { return name === "support" }
    function specNeedsMatrices(name) { return name === "radialvar" }

    // BSMEM: the entropy prior is a spec entry (["mem", image]) but the stopping criterion is
    // a 4-element Int vector, which is not a regulariser in any sense the list could express.
    property var  bsmemMethod: [4, 1, 1, 2]
    property bool bsmemUsePrior: false
    property string bsmemPriorPath: ""
    property int  bsmemMaxiter: 200

    // BSDMM: proximal blocks selected by weight, with the operator chosen by name. A zero
    // weight is how a block is switched off, so there is no separate enable flag.
    property real   bsdmmMuTv: 0.0
    property real   bsdmmMuGroup: 0.0
    // Non-zero by default. Centering is what pins the image against the translation degeneracy
    // of the bispectrum, and ADMM needs at least one proximal block to split across, so a run
    // with it at zero is either unpinned or has nothing to solve. 1e2 is what the shipped BSDMM
    // demos use alongside mu_reg = 1e5.
    property real   bsdmmMuCen: 1e2
    property string bsdmmRegType: "tv"          // tv | l2smooth
    property string bsdmmGroupType: "sparsity"  // sparsity | tv
    property int    bsdmmMaxiter: 200

    // VMLMB / SPARCO run controls
    property int vmlmbMaxiter: 200
    property int sparcoRounds: 5

    // SQUEEZE (§5.4a)
    property int  nelements: 0            // 0 ⇒ default_nelements(data)
    property int  niter: 200
    property int  nchains: 4
    property real tmin: 1.0
    property real fAnywhere: 0.1
    property real fCopycat: 0.2
    property real centMult: 1.0
    property bool autoCentering: true
    property string squeezePriorPath: ""

    // Parallel tempering, on the Pigeons side
    property int nreplicas: 10
    property int nrounds: 10

    // SPARCO parametric part, feeding either FlatModel (hybrid) or SqueezeSparco (annealing)
    property real sparcoFStar: 0.5
    property real sparcoUd: 0.0
    property real sparcoEnvIndx: 0.0
    property real sparcoLambda0: 1.65e-6      // metres
    property real sparcoFBg: 0.0
    property real sparcoBgIndx: 0.0
    property bool sparcoFreeFStar: true
    property bool sparcoFreeUd: false
    property bool sparcoFreeEnvIndx: false
    property bool sparcoFreeFBg: false
    property bool sparcoFreeBgIndx: false

    // VI (OIVI): there are no regularisers at all — the correlated-field prior IS the
    // regularisation, so that panel shows prior hyperparameters where the others show μ.
    property string viEngine: "mgvi"          // map | mgvi | geovi | hybrid
    property int    viSamples: 8
    property real   cfSlopeMean: -3.0
    property real   cfSlopeStd: 0.5
    property real   cfFluctMean: 1.0
    property real   cfFluctStd: 0.5
    property real   cfFlexMean: 1.0
    property real   cfFlexStd: 0.5
    property bool   cfUseIwp: false
    property real   cfAspMean: 1.0
    property real   cfAspStd: 0.5
    property real   cfOffsetMean: 0.0
    property real   cfOffsetStd: 0.1
    property bool   cfUseOffset: false

    // ── run state ─────────────────────────────────────────────────────────────
    property bool   running: false
    property real   progress: -1          // < 0 ⇒ unknown, shown as an indeterminate bar
    property string statusText: "idle"
    property string roundText: ""         // the stochastic engines report per round, not per iteration

    readonly property string blockedReason:
          !engineAvailableFor(engine)                        ? engineReasonFor(engine)
        : !anyObservable                                     ? "no observables selected"
        : (nx <= 0 || pixsize <= 0)                          ? "image geometry is not set"
        : (startImage === "fits" && startImagePath === "")   ? "no starting image chosen"
        : ""
    readonly property bool canRun: !running && blockedReason === ""

    // ── results (§5.5) ────────────────────────────────────────────────────────
    property bool   hasResult: false
    property bool   resultIsPosterior: false   // set by Julia; only the stochastic engines produce one
    property string resultMode: "mean"         // mean | sigma | sample
    property int    sampleCount: 0
    property int    sampleIndex: 0

    // ── diagnostics (§5.5) ────────────────────────────────────────────────────
    //
    // Thresholds live here, once, from demos/pigeons_diagnostics.md. Julia pushes raw numbers
    // and nothing else; deciding red/amber/green in Julia would put the rules in two places.
    property bool hasDiagnostics: false
    property int  diagRestarts: 0
    property real diagLambda: 0
    property real diagSwapMin: 0
    property real diagSwapMean: 0
    property real diagRhoMax: 0
    property real diagRhoMean: 0
    property string diagRhoTrend: "→"
    property real diagExplorerMin: 0
    property real diagExplorerMean: 0
    property real diagLogZ: 0
    property real diagLogZDelta: 0
    property var  diagSwapPerRung: []          // one acceptance per ladder rung, for the bottleneck view
    // annealing
    property real diagAcceptMin: 0
    property real diagAcceptMean: 0
    property real diagTemperature: 0
    property real diagChi2rBest: 0
    property real diagFlatChi2r: 0
    property int  diagBestChain: 0
    // VI
    property real diagElbo: 0
    property real diagElboDelta: 0

    function fmt(v, d) { return isNaN(v) ? "—" : Number(v).toFixed(d === undefined ? 3 : d) }

    readonly property var diagRows:
          engine === "tempering" ? [
            { label: "restarts",
              value: hasDiagnostics ? String(diagRestarts) : "—",
              state: !hasDiagnostics ? "unknown" : (diagRestarts > 0 ? "ok" : "bad"),
              tip: "round trips between the coldest and hottest chain. > 0 is necessary for reliable results; 0 means the chains are stuck." },
            { label: "Λ",
              value: hasDiagnostics ? fmt(diagLambda, 2) + "  → suggests " + Math.max(2, Math.round(2 * diagLambda)) + " chains" : "—",
              state: !hasDiagnostics ? "unknown" : (diagLambda >= nreplicas - 1 ? "bad" : "ok"),
              tip: "communication barrier. Set n_chains ≈ 2Λ; Λ ≈ n_chains − 1 means the ladder is saturated." },
            { label: "swap α",
              value: hasDiagnostics ? "min " + fmt(diagSwapMin, 2) + " · mean " + fmt(diagSwapMean, 2) : "—",
              state: !hasDiagnostics ? "unknown"
                   : (diagSwapMin <= 0 ? "bad"
                   : ((diagSwapMean < 0.3 || diagSwapMean > 0.5 || diagSwapMin < 0.3) ? "warn" : "ok")),
              tip: "swap acceptance, wanted in 0.3–0.5. min(α) = 0 marks one specific bottleneck rung — see the ladder below." },
            { label: "autocorr |ρ|",
              value: hasDiagnostics ? "max " + fmt(diagRhoMax, 2) + " · mean " + fmt(diagRhoMean, 2) + "  " + diagRhoTrend : "—",
              state: !hasDiagnostics ? "unknown" : (diagRhoMax > 0.9 ? "warn" : "ok"),
              tip: "near 1 means consecutive samples are the same sample: more scans are needed before the ensemble is independent." },
            { label: "explorer α",
              value: hasDiagnostics ? "min " + fmt(diagExplorerMin, 2) + " · mean " + fmt(diagExplorerMean, 2) : "—",
              state: !hasDiagnostics ? "unknown" : (diagExplorerMin < 0.3 ? "warn" : "ok"),
              tip: "AutoMALA self-tunes toward 0.5–0.6; below 0.3 the explorer is not moving." },
            { label: "log(Z₁/Z₀)",
              value: hasDiagnostics ? fmt(diagLogZ, 3) + "  (Δ " + fmt(diagLogZDelta, 3) + " this round)" : "—",
              state: !hasDiagnostics ? "unknown" : (Math.abs(diagLogZDelta) > 0.1 ? "warn" : "ok"),
              tip: "must stabilise across rounds. Once stable it also gives model comparison between two reconstructions." }
          ]
        : isSqueeze ? [
            { label: "chain acceptance",
              value: hasDiagnostics ? "min " + fmt(diagAcceptMin, 2) + " · mean " + fmt(diagAcceptMean, 2) : "—",
              state: !hasDiagnostics ? "unknown" : (diagAcceptMin < 0.05 ? "bad" : (diagAcceptMin < 0.2 ? "warn" : "ok")),
              tip: "collapsing acceptance is annealing's version of a stuck ladder: the chain has stopped exploring." },
            { label: "temperature",
              value: hasDiagnostics ? fmt(diagTemperature, 3) + "  (tmin " + fmt(tmin, 3) + ")" : "—",
              state: "unknown",
              tip: "current point on the schedule, against the floor asked for." },
            { label: "χ²r best chain",
              value: hasDiagnostics ? fmt(diagChi2rBest, 3) + "  (flat image scores " + fmt(diagFlatChi2r, 3) + ")" : "—",
              state: !hasDiagnostics ? "unknown" : (diagChi2rBest < diagFlatChi2r ? "ok" : "bad"),
              tip: "flat_chi2r is the free baseline: a reconstruction that does not beat a flat image has found nothing." },
            { label: "best chain",
              value: hasDiagnostics ? String(diagBestChain) + " of " + String(nchains) : "—",
              state: "unknown",
              tip: "which chain the returned image came from; the others are the ensemble behind the σ map." }
          ]
        : engine === "vi" ? [
            { label: "ELBO",
              value: hasDiagnostics ? fmt(diagElbo, 3) + "  (Δ " + fmt(diagElboDelta, 3) + ")" : "—",
              state: !hasDiagnostics ? "unknown" : (Math.abs(diagElboDelta) > 0.1 ? "warn" : "ok"),
              tip: "the variational bound must settle before the posterior means anything." },
            { label: "posterior samples",
              value: hasDiagnostics ? String(sampleCount) : "—",
              state: "unknown",
              tip: "how many samples back the σ map and the scrubber." }
          ]
        : []

    function lightColor(state) {
        switch (state) {
        case "ok":   return "#2e7d32"
        case "warn": return "#ef6c00"
        case "bad":  return "#c62828"
        default:     return "#bdbdbd"
        }
    }



    // ── observable boxes follow the file, not the user's last session ─────────
    //
    // Available comes back ticked, absent comes back unticked: "ticked but absent" is the
    // state that must be unrepresentable. Unticking something present is legitimate and
    // common — V²-only, or dropping a suspect T3amp — so nothing here prevents it.
    // wire: after a dataset loads Julia sets the have* flags from the VALID-POINT counts
    //       (nv2, nt3amp, nt3phi, ndiffphase — not array lengths), which resets these boxes;
    //       on Run it turns the ticks back into weights = [V², T3amp, T3φ] and the
    //       use_diffphases kwarg, which is where diffvis goes since it is not a weight.
    function resetObservables() {
        cvisBox.checked    = false
        v2Box.checked      = haveV2
        t3ampBox.checked   = haveT3amp
        t3phiBox.checked   = haveT3phi
        diffvisBox.checked = haveDiffvis
        fluxBox.checked    = false
    }
    onHaveV2Changed:      resetObservables()
    onHaveT3ampChanged:   resetObservables()
    onHaveT3phiChanged:   resetObservables()
    onHaveDiffvisChanged: resetObservables()

    // ── running a reconstruction ─────────────────────────────────────────────
    //
    // The call holds Qt's thread for as long as the optimiser takes, so the status has to be
    // on screen BEFORE it starts — hence the timer. Without it the window freezes with the
    // button still reading "Run" and nothing saying why.
    // The VMLMB regulariser rows as `name,mu[,extra];name,mu` — one string, because that is
    // what crosses into Julia cleanly. Only this engine's rows are sent: the other engines
    // take their regularisation in shapes that are not a named list at all.
    // The iteration count the chosen engine actually reads. They are not the same control:
    // VMLMB counts gradient steps, BSMEM MEM iterations, ADMM splittings and SQUEEZE sweeps,
    // and a shared "maxiter" would silently mean four different amounts of work.
    function engineMaxiter() {
        if (engine === "bsmem") return bsmemMaxiter
        if (engine === "bsdmm") return bsdmmMaxiter
        if (isSqueeze)          return niter
        return vmlmbMaxiter
    }

    // Everything the chosen engine needs beyond the shared geometry and observables, as
    // `key\tvalue` lines. Only the current engine's settings are sent: an engine that received
    // another's keys would ignore them, which is fine, but the console line that echoes them
    // would then describe settings that had no effect.
    function engineOptions() {
        var o = []
        function put(k, v) { o.push(k + "\t" + v) }

        if (engine === "bsmem") {
            put("method1", bsmemMethod[0]); put("method2", bsmemMethod[1])
            put("method3", bsmemMethod[2]); put("method4", bsmemMethod[3])
            put("maxiter", bsmemMaxiter)
            if (bsmemUsePrior && bsmemPriorPath.length > 0) put("prior", bsmemPriorPath)

        } else if (engine === "bsdmm") {
            put("mu_reg", bsdmmMuTv); put("mu_cen", bsdmmMuCen)
            put("reg_type", bsdmmRegType); put("maxiter", bsdmmMaxiter)

        } else if (engine === "sparco") {
            put("lambda0", sparcoLambda0); put("f_star", sparcoFStar)
            put("f_bg", sparcoFBg); put("env_indx", sparcoEnvIndx)
            put("ud", sparcoUd); put("rounds", sparcoRounds)

        } else if (isSqueeze) {
            put("nelements", nelements); put("niter", niter); put("nchains", nchains)
            put("tmin", tmin); put("f_anywhere", fAnywhere); put("f_copycat", fCopycat)
            put("cent_mult", centMult); put("auto_centering", autoCentering ? 1 : 0)
            put("seed", startSeed)
            if (squeezePriorPath.length > 0) put("prior", squeezePriorPath)
            if (engine === "squeeze_sparco") {
                put("lambda0", sparcoLambda0); put("f_star", sparcoFStar)
                put("ud", sparcoUd); put("env_indx", sparcoEnvIndx)
                put("f_bg", sparcoFBg); put("bg_indx", sparcoBgIndx)
                put("free_f_star", sparcoFreeFStar ? 1 : 0)
                put("free_ud", sparcoFreeUd ? 1 : 0)
                put("free_env_indx", sparcoFreeEnvIndx ? 1 : 0)
                put("free_f_bg", sparcoFreeFBg ? 1 : 0)
                put("free_bg_indx", sparcoFreeBgIndx ? 1 : 0)
            }
        }
        return o.join("\n")
    }

    function regulariserSpec() {
        // SQUEEZE's names are spelled like oichi2.jl's but act on the integer histogram, so it
        // has its own closed list and its own model behind it.
        if (isSqueeze) {
            var sq = []
            for (var s = 0; s < squeezeRegModel.count; ++s) {
                var q = squeezeRegModel.get(s)
                if (!q.use) continue
                sq.push(q.name + "," + q.lambda)
            }
            return sq.join(";")
        }
        // BSMEM ignores the spec list -- its entropy mode is the `method` vector, and its prior
        // travels as an option -- so sending it regularisers would describe a run it did not do.
        if (engine !== "vmlmb" && engine !== "sparco") return ""
        var out = []
        for (var i = 0; i < regModel.count; ++i) {
            var r = regModel.get(i)
            if (!r.name || r.name.length === 0) continue
            if (r.on === false) continue        // unticked: left out of this run, not deleted
            var row = r.name + "," + r.mu
            // `l1l2` is the one name that REQUIRES its second argument; sending it without
            // one would error inside reconstruct rather than here.
            if (r.extra !== undefined && String(r.extra).length > 0) row += "," + r.extra
            out.push(row)
        }
        return out.join(";")
    }

    // ── the Fourier plan ─────────────────────────────────────────────────────
    //
    // Rebuilt whenever nx, the pixel size or the transform changes, because those three ARE
    // the plan — a run with a stale one would transform against a geometry the panel no longer
    // shows. Debounced: building an NFFT plan is the expensive part of setting a
    // reconstruction up, and a spin box emits a value per click.
    property string planText: "plan: —"

    function refreshPlan() { planTimer.restart() }

    Timer {
        id: planTimer
        interval: 400; repeat: false
        onTriggered: {
            root.planText = Julia.shell_ft_setup(root.nx, root.pixsize, root.ftMode)
            root.consoleChanged()
        }
    }

    onNxChanged:      refreshPlan()
    onPixsizeChanged: refreshPlan()
    onFtModeChanged:  refreshPlan()

    // ── per-observable reduced chi2 ──────────────────────────────────────────
    //
    // One number for the whole fit hides which observable is being fitted and which is being
    // sacrificed, and that is usually the question. An unticked observable still gets a row:
    // its chi2 says how well the image predicts data it never saw.
    ListModel { id: breakdownModel }

    function refreshBreakdown() {
        breakdownModel.clear()
        var rows = Julia.shell_chi2_breakdown()
        if (rows.length === 0) return
        var lines = rows.split("\n")
        for (var i = 0; i < lines.length; ++i) {
            var f = lines[i].split("\t")
            if (f.length === 4)
                breakdownModel.append({ obs: f[0], chi2r: f[1], npts: f[2], used: f[3] === "1" })
        }
    }

    // Raised when Julia has written to the shared console and the window should re-read it.
    signal consoleChanged()

    property bool _continueRun: false

    function startRun(fromPrevious) {
        if (!canRun) return
        _continueRun = fromPrevious === true
        running = true
        statusText = _continueRun ? "continuing…" : "reconstructing…"
        runTimer.restart()
    }

    Timer {
        id: runTimer
        interval: 120; repeat: false
        onTriggered: {
            var line = Julia.shell_reconstruct(root.nx, root.pixsize, root.ftMode,
                                               root.startImage,
                                               root.useV2, root.useT3amp, root.useT3phi,
                                               root.engineMaxiter(),
                                               root.regulariserSpec(),
                                               root.startFwhm, root.startSeed,
                                               root._continueRun,
                                               root.engine, root.engineOptions())
            root.statusText = line
            root.hasResult = line.indexOf("chi2r") >= 0
            root.refreshBreakdown()
            root.running = false
            // The console pane belongs to the window, and the reconstruction wrote to it.
            root.consoleChanged()
            // Makie draws on demand here; without this the new image is not painted until
            // something else happens to invalidate the area.
            imageArea.update()
        }
    }

    // Which page of the regulariser stack an engine wants. Annealing and tempering share one:
    // tempering IS the annealing sampler with a Pigeons ladder around it.
    function regPageFor(key) {
        switch (key) {
        case "bsmem":     return 1
        case "bsdmm":     return 2
        case "squeeze":        return 3
        case "tempering":      return 3
        case "squeeze_sparco": return 3
        case "vi":        return 4
        case "sparco":    return 5
        default:          return 0
        }
    }

    // ── small reusable pieces ─────────────────────────────────────────────────

    // A disabled control receives no hover events, so its own ToolTip never appears — which is
    // precisely the control whose greying most needs explaining. This sits over the control as
    // an enabled sibling and exists only while there is something to say; with NoButton it
    // does not swallow the clicks of a control that is live.
    component ReasonTip: MouseArea {
        property string reason: ""
        anchors.fill: parent
        visible: reason.length > 0
        hoverEnabled: true
        acceptedButtons: Qt.NoButton
        ToolTip.visible: containsMouse
        ToolTip.text: reason
        ToolTip.delay: 400
    }

    // A numeric entry that does not fight the user: the field owns the text while it is being
    // edited and only reports a parsed, clamped value on commit, so a binding into `value`
    // survives typing and Julia can still push a suggested number in.
    component NumField: TextField {
        property real value: 0
        property real minimum: -1.0e12
        property real maximum:  1.0e12
        signal committed(real v)

        function display(x) {
            if (x === 0) return "0"
            var a = Math.abs(x)
            return (a >= 1.0e-4 && a < 1.0e6) ? String(Number(x.toPrecision(8))) : x.toExponential(3)
        }

        horizontalAlignment: TextInput.AlignRight
        selectByMouse: true
        inputMethodHints: Qt.ImhFormattedNumbersOnly
        validator: DoubleValidator {
            bottom: minimum; top: maximum
            notation: DoubleValidator.ScientificNotation
        }
        Component.onCompleted: text = display(value)
        onValueChanged: if (!activeFocus) text = display(value)
        onEditingFinished: {
            var v = parseFloat(text)
            if (isNaN(v)) { text = display(value); return }
            v = Math.min(maximum, Math.max(minimum, v))
            text = display(v)
            committed(v)
        }
    }

    // ── models ────────────────────────────────────────────────────────────────

    ListModel {
        id: engineModel
        ListElement { key: "vmlmb";     name: "VMLMB";                                 entry: "reconstruct" }
        ListElement { key: "bsmem";     name: "BSMEM";                                 entry: "reconstruct_bsmem" }
        ListElement { key: "bsdmm";     name: "BSDMM";                                 entry: "reconstruct_bsdmm" }
        ListElement { key: "sparco";    name: "SPARCO hybrid";                         entry: "reconstruct_hybrid" }
        // SQUEEZE appears three times because the three are used as three different
        // engines, not as one engine with settings: annealing and tempering are separate
        // entry points, and the SPARCO variant additionally samples a parametric component
        // alongside the image. Sharing one row and a checkbox would bury that.
        ListElement { key: "squeeze";        name: "Squeeze (Annealing)";
                      entry: "reconstruct_squeeze" }
        ListElement { key: "tempering";      name: "Squeeze (Tempering)";
                      entry: "reconstruct_squeeze_tempered" }
        ListElement { key: "squeeze_sparco"; name: "Squeeze (Annealing + SPARCO)";
                      entry: "reconstruct_squeeze  (model = SqueezeSparco)" }
        ListElement { key: "vi";        name: "Variational inference (OIVI)";          entry: "reconstruct_map / mgvi / geovi" }
    }

    // Empty on purpose. Positivity is already in force — it is VMLMB's `lower = 0`, not a
    // regulariser — so an empty list is an unregularised maximum-likelihood reconstruction,
    // which is a real one and the right thing to see before anything is added. A pre-filled
    // row would be applied by the first Run without the user having asked for it.
    ListModel {
        id: regModel
    }

    // SQUEEZE's names are spelled like oichi2.jl's but they act on the integer histogram, so
    // tv/entropy/compactness are different maths here and the list is closed at six.
    ListModel {
        id: squeezeRegModel
        ListElement { name: "l0";           use: false; lambda: 1.0; path: "" }
        ListElement { name: "tv";           use: false; lambda: 1.0; path: "" }
        ListElement { name: "entropy";      use: false; lambda: 1.0; path: "" }
        ListElement { name: "compactness";  use: false; lambda: 1.0; path: "" }
        ListElement { name: "centering";    use: true;  lambda: 1.0; path: "" }
        ListElement { name: "priorimage";   use: false; lambda: 1.0; path: "" }
    }

    // ── file pickers ──────────────────────────────────────────────────────────

    // FilePicker throughout: QtQuick.Dialogs.FileDialog leaves its window mapped on some
    // systems even after reporting itself closed -- see FilePicker.qml.
    FilePicker {
        id: startImageDialog
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Starting image (FITS)"
        filters: [{ label: "FITS images (*.fits *.fit)", patterns: "*.fits,*.fit" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) { root.startImagePath = path }
    }

    FilePicker {
        id: bsmemPriorDialog
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "MaxEnt prior image (FITS)"
        filters: [{ label: "FITS images (*.fits *.fit)", patterns: "*.fits,*.fit" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) { root.bsmemPriorPath = path }
    }

    FilePicker {
        id: squeezePriorDialog
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Squeeze prior image (FITS)"
        // A zero pixel here is a hard mask: no flux quantum may sit there.
        filters: [{ label: "FITS images (*.fits *.fit)", patterns: "*.fits,*.fit" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) { root.squeezePriorPath = path }
    }

    FilePicker {
        id: specImageDialog
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Regulariser image (FITS)"
        property int rowIndex: -1
        filters: [{ label: "FITS images (*.fits *.fit)", patterns: "*.fits,*.fit" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) {
            if (rowIndex >= 0) regModel.setProperty(rowIndex, "path", path)
        }
    }

    // ── the ["name", μ, extra...] spec editor, used by VMLMB and by SPARCO ─────
    //
    // One Component, two Loaders, one model: μ weights genuinely carry over between those two
    // engines, so editing the list on either page edits the same list.
    Component {
        id: specRegPanel

        ColumnLayout {
            spacing: dp(4)

            Repeater {
                model: regModel

                delegate: RowLayout {
                    id: specRow
                    required property int index
                    required property var model
                    Layout.fillWidth: true
                    spacing: dp(6)

                    // Off rather than removed. Trying a reconstruction with and without a term
                    // is the normal way to find out what it is doing, and deleting the row to
                    // do that loses the weight that was tuned to get there.
                    CheckBox {
                        checked: specRow.model.on === undefined ? true : specRow.model.on
                        ToolTip.visible: hovered
                        ToolTip.text: checked ? "in use — untick to leave it out of the run"
                                              : "not used by the next run"
                        onToggled: regModel.setProperty(specRow.index, "on", checked)
                    }

                    ComboBox {
                        Layout.preferredWidth: dp(180)
                        enabled: specRow.model.on === undefined ? true : specRow.model.on
                        model: root.specNames
                        currentIndex: root.specNames.indexOf(specRow.model.name)
                        onActivated: regModel.setProperty(specRow.index, "name", root.specNames[currentIndex])
                    }
                    Label { text: "μ" }
                    NumField {
                        Layout.preferredWidth: dp(80)
                        minimum: 0
                        value: specRow.model.mu
                        onCommitted: (v) => regModel.setProperty(specRow.index, "mu", v)
                    }
                    Label {
                        text: root.specExtraLabel(specRow.model.name)
                        visible: text.length > 0
                        color: root.specExtraRequired(specRow.model.name) ? "#000" : "#666"
                    }
                    NumField {
                        Layout.preferredWidth: dp(80)
                        visible: root.specExtraLabel(specRow.model.name).length > 0
                        value: specRow.model.extra
                        onCommitted: (v) => regModel.setProperty(specRow.index, "extra", v)
                    }
                    Button {
                        text: specRow.model.path === "" ? "prior…" : "prior ✓"
                        visible: root.specNeedsImage(specRow.model.name)
                        onClicked: { specImageDialog.rowIndex = specRow.index; specImageDialog.openAt("") }
                    }
                    Label {
                        // radialvar takes two matrices, H and G, which no numeric field can carry.
                        text: "H/G set in Julia"
                        color: "#888"; font.pointSize: pt(baseFontPt - 2)
                        visible: root.specNeedsMatrices(specRow.model.name)
                        // wire: Julia supplies the radialvar H and G matrices for this row
                    }
                    Item { Layout.fillWidth: true }
                    Button {
                        text: "−"
                        implicitWidth: dp(28)
                        implicitHeight: dp(24)
                        ToolTip.visible: hovered
                        ToolTip.text: "remove this regulariser"
                        onClicked: regModel.remove(specRow.index)
                    }
                }
            }

            RowLayout {
                Layout.fillWidth: true
                spacing: dp(8)
                Button {
                    text: "+ regularizer"
                    onClicked: regModel.append({ "on": true, "name": "tv", "mu": 1.0,
                                                 "extra": 1.0e-8, "path": "" })
                }
            }
        }
    }

    // ── the SPARCO parametric part, shared by the hybrid and by SqueezeSparco ──
    Component {
        id: sparcoPanel

        GridLayout {
            columns: 3
            columnSpacing: dp(8)
            rowSpacing: dp(4)

            Label { text: "f_star"; }
            NumField { Layout.preferredWidth: dp(90); minimum: 0; maximum: 1
                       value: root.sparcoFStar; onCommitted: (v) => root.sparcoFStar = v }
            CheckBox { text: "free"; checked: root.sparcoFreeFStar; onToggled: root.sparcoFreeFStar = checked }

            Label { text: "ud (mas)" }
            NumField { Layout.preferredWidth: dp(90); minimum: 0
                       value: root.sparcoUd; onCommitted: (v) => root.sparcoUd = v }
            CheckBox { text: "free"; checked: root.sparcoFreeUd; onToggled: root.sparcoFreeUd = checked }

            Label { text: "env_indx" }
            NumField { Layout.preferredWidth: dp(90)
                       value: root.sparcoEnvIndx; onCommitted: (v) => root.sparcoEnvIndx = v }
            CheckBox { text: "free"; checked: root.sparcoFreeEnvIndx; onToggled: root.sparcoFreeEnvIndx = checked }

            Label { text: "λ₀ (m)" }
            NumField { Layout.preferredWidth: dp(90); minimum: 1.0e-9
                       value: root.sparcoLambda0; onCommitted: (v) => root.sparcoLambda0 = v }
            Label { text: "always fixed"; color: "#888"; font.pointSize: pt(baseFontPt - 2) }

            Label { text: "f_bg" }
            NumField { Layout.preferredWidth: dp(90); minimum: 0; maximum: 1
                       value: root.sparcoFBg; onCommitted: (v) => root.sparcoFBg = v }
            CheckBox { text: "free"; checked: root.sparcoFreeFBg; onToggled: root.sparcoFreeFBg = checked }

            Label { text: "bg_indx" }
            NumField { Layout.preferredWidth: dp(90)
                       value: root.sparcoBgIndx; onCommitted: (v) => root.sparcoBgIndx = v }
            CheckBox { text: "free"; checked: root.sparcoFreeBgIndx; onToggled: root.sparcoFreeBgIndx = checked }

            // wire: these six become FlatModel params for reconstruct_hybrid, or
            // SqueezeSparco(; f_star, ud, env_indx, lambda0, f_bg, bg_indx, free) for annealing
        }
    }

    // ══════════════════════════════════════════════════════════════════════════
    //  layout
    // ══════════════════════════════════════════════════════════════════════════

    RowLayout {
        anchors.fill: parent
        spacing: 0

        // ── left: everything that defines the reconstruction ──────────────────
        ColumnLayout {
            Layout.preferredWidth: dp(430)
            Layout.minimumWidth: dp(340)
            Layout.fillHeight: true
            spacing: dp(6)

            ScrollView {
                id: setupScroll
                Layout.fillWidth: true
                Layout.fillHeight: true
                contentWidth: availableWidth
                clip: true

                ColumnLayout {
                    width: setupScroll.availableWidth
                    spacing: dp(8)

                    // ── setup (§5.1) ──────────────────────────────────────────
                    GroupBox {
                        title: "Setup"
                        Layout.fillWidth: true
                        Layout.margins: dp(6)

                        GridLayout {
                            anchors.fill: parent
                            columns: 3
                            columnSpacing: dp(8)
                            rowSpacing: dp(4)

                            Label { text: "nx (pixels)" }
                            SpinBox {
                                id: nxBox
                                from: 16; to: 2048; stepSize: 16
                                // Bound, not initialised: Julia pushes a data-driven nx after a
                                // load, and a box that only took its value at construction
                                // would show one number while the run used another.
                                value: root.nx
                                editable: true
                                onValueModified: root.nx = value
                            }
                            Label {
                                text: "FOV " + root.fov.toFixed(2) + " mas"
                                color: "#666"
                            }

                            Label { text: "pixel size (mas)" }
                            NumField {
                                id: pixField
                                Layout.preferredWidth: dp(100)
                                minimum: 1.0e-6
                                value: root.pixsize
                                // Julia pushes auto_pixsize(data) into root.pixsize after a
                                // load; this follows it.
                                onCommitted: (v) => {
                                    root.pixsize = v
                                    // wire: if scaleRegs, Julia calls scale_regularizers(regs, pixsize)
                                    //       so a tuned μ survives the resolution change
                                }
                            }
                            CheckBox {
                                id: scaleRegsBox
                                text: "scale regularisers"
                                checked: true
                                onToggled: root.scaleRegs = checked
                                ToolTip.visible: hovered
                                ToolTip.text: "rescale the μ weights through scale_regularizers() whenever the pixel size changes, so a tuned regularisation survives a resolution change"
                            }

                            Label { text: "Fourier transform" }
                            ComboBox {
                                id: ftBox
                                Layout.preferredWidth: dp(100)
                                model: ["nfft", "dft"]
                                onActivated: root.ftMode = currentText
                            }
                            Label {
                                // What a Run would actually transform with. The plan is what
                                // decides the representable field of view and resolution, and
                                // it is rebuilt silently whenever the geometry changes — so it
                                // is worth stating rather than leaving to be inferred.
                                text: root.planText
                                color: root.planText.indexOf("!") === 0 ? "#c62828" : "#888"
                                font.pointSize: pt(baseFontPt - 2)
                                elide: Text.ElideRight
                                Layout.fillWidth: true
                                ToolTip.visible: planHover.hovered && root.planText.length > 0
                                ToolTip.text: root.planText
                                HoverHandler { id: planHover }
                                //       (mode / size / channels / uv counts)
                            }

                            Label { text: "starting image" }
                            ComboBox {
                                id: startBox
                                Layout.preferredWidth: dp(100)
                                model: ["Dirac", "Gaussian", "Random", "FITS file"]
                                currentIndex: root.startImage === "gaussian" ? 1
                                            : root.startImage === "random"   ? 2
                                            : root.startImage === "fits"     ? 3 : 0
                                onActivated: root.startImage =
                                    ["dirac", "gaussian", "random", "fits"][currentIndex]
                            }
                            RowLayout {
                                spacing: dp(6)
                                Layout.fillWidth: true
                                NumField {
                                    Layout.preferredWidth: dp(80)
                                    visible: root.startImage === "gaussian"
                                    // Negative is the "choose it for me" sentinel. Zero is not:
                                    // zero IS a width, and a field showing 0 reads as a value
                                    // someone meant rather than as one nobody set.
                                    minimum: -1
                                    value: root.startFwhm
                                    onCommitted: (v) => root.startFwhm = v
                                }
                                Label {
                                    visible: root.startImage === "gaussian"
                                    text: root.startFwhm > 0 ? "mas FWHM" : "mas FWHM (−1 ⇒ from FOV)"
                                    color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                }
                                // A random start answers "does the answer depend on where I
                                // started", which a smooth one cannot. Its seed is part of the
                                // reconstruction, not a detail: without it the run cannot be
                                // repeated and an exported script does not reproduce its figure.
                                Label {
                                    visible: root.startImage === "random"
                                    text: "seed"
                                    color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                }
                                SpinBox {
                                    visible: root.startImage === "random"
                                    Layout.preferredWidth: dp(90)
                                    from: 0; to: 999999; editable: true
                                    value: root.startSeed
                                    onValueModified: root.startSeed = value
                                }
                                Button {
                                    visible: root.startImage === "fits"
                                    text: "Choose…"
                                    onClicked: startImageDialog.openAt("")
                                }
                                Label {
                                    visible: root.startImage === "fits"
                                    Layout.fillWidth: true
                                    text: root.startImagePath === "" ? "no file chosen" : root.startImagePath
                                    color: root.startImagePath === "" ? "#c62828" : "#666"
                                    elide: Text.ElideLeft
                                    font.pointSize: pt(baseFontPt - 2)
                                }
                                Item { Layout.fillWidth: true; visible: root.startImage === "dirac" }
                            }
                        }
                    }

                    // ── observables (§1.8) ────────────────────────────────────
                    GroupBox {
                        title: "Observables"
                        Layout.fillWidth: true
                        Layout.margins: dp(6)

                        ColumnLayout {
                            anchors.fill: parent
                            spacing: dp(4)

                            RowLayout {
                                Layout.fillWidth: true
                                spacing: dp(12)

                                // cvis and flux are dead for imaging whatever the file holds, so they
                                // are marked † and carry the engine's reason, not the file's.
                                Item {
                                    implicitWidth: cvisBox.implicitWidth
                                    implicitHeight: cvisBox.implicitHeight
                                    CheckBox { id: cvisBox; text: "cvis †"; checked: false; enabled: false }
                                    ReasonTip { reason: root.engineExclusion }
                                }
                                Item {
                                    implicitWidth: v2Box.implicitWidth
                                    implicitHeight: v2Box.implicitHeight
                                    CheckBox { id: v2Box; text: "v2"; checked: root.haveV2; enabled: root.haveV2 }
                                    ReasonTip { reason: root.haveV2 ? "" : "this file contains no V² points" }
                                }
                                Item {
                                    implicitWidth: t3ampBox.implicitWidth
                                    implicitHeight: t3ampBox.implicitHeight
                                    CheckBox { id: t3ampBox; text: "t3amp"; checked: root.haveT3amp; enabled: root.haveT3amp }
                                    ReasonTip { reason: root.haveT3amp ? "" : "this file contains no T3amp points" }
                                }
                                Item {
                                    implicitWidth: t3phiBox.implicitWidth
                                    implicitHeight: t3phiBox.implicitHeight
                                    CheckBox { id: t3phiBox; text: "t3phi"; checked: root.haveT3phi; enabled: root.haveT3phi }
                                    ReasonTip { reason: root.haveT3phi ? "" : "this file contains no T3phi points" }
                                }
                                Item {
                                    implicitWidth: diffvisBox.implicitWidth
                                    implicitHeight: diffvisBox.implicitHeight
                                    // Not a weight at all: differential phases enter reconstruct() through
                                    // use_diffphases = true, and only for a polychromatic dataset.
                                    CheckBox { id: diffvisBox; text: "diffvis"; checked: root.haveDiffvis; enabled: root.haveDiffvis }
                                    ReasonTip { reason: root.haveDiffvis ? "" : "this file has no differential phases (needs a polychromatic dataset)" }
                                }
                                Item {
                                    implicitWidth: fluxBox.implicitWidth
                                    implicitHeight: fluxBox.implicitHeight
                                    CheckBox { id: fluxBox; text: "flux †"; checked: false; enabled: false }
                                    ReasonTip { reason: root.engineExclusion }
                                }
                                Item { Layout.fillWidth: true }
                            }

                            Label {
                                Layout.fillWidth: true
                                wrapMode: Text.WordWrap
                                text: "† permanently unsupported by the imaging engines — nothing to do with this file.\n" +
                                      "weights here are 3 long, [V², T3amp, T3φ]; model fitting's are 7 long."
                                color: "#888"; font.pointSize: pt(baseFontPt - 2)
                            }
                        }
                    }

                    // ── engine (§5.2) ─────────────────────────────────────────
                    GroupBox {
                        title: "Engine"
                        Layout.fillWidth: true
                        Layout.margins: dp(6)

                        ColumnLayout {
                            anchors.fill: parent
                            spacing: dp(4)

                            ComboBox {
                                id: engineBox
                                Layout.fillWidth: true
                                model: engineModel
                                textRole: "name"
                                valueRole: "key"
                                onActivated: root.engine = engineModel.get(currentIndex).key
                                // wire: Julia probes for Pigeons and OIVI and sets
                                //       pigeonsAvailable / oiviAvailable, which greys these entries

                                delegate: ItemDelegate {
                                    id: engineDelegate
                                    required property int index
                                    required property var model
                                    width: engineBox.width
                                    enabled: root.engineAvailableFor(engineDelegate.model.key)
                                    highlighted: engineBox.highlightedIndex === engineDelegate.index
                                    contentItem: RowLayout {
                                        spacing: dp(8)
                                        Label {
                                            Layout.fillWidth: true
                                            text: engineDelegate.model.name
                                            elide: Text.ElideRight
                                            opacity: engineDelegate.enabled ? 1.0 : 0.5
                                        }
                                        Label {
                                            text: root.engineReasonFor(engineDelegate.model.key)
                                            visible: text.length > 0
                                            color: "#c62828"
                                            font.pointSize: pt(baseFontPt - 2)
                                        }
                                    }
                                }
                            }

                            Label {
                                Layout.fillWidth: true
                                text: engineBox.currentIndex >= 0
                                      ? engineModel.get(engineBox.currentIndex).entry + "()" : ""
                                color: "#666"; font.family: "monospace"
                                font.pointSize: pt(baseFontPt - 2)
                                elide: Text.ElideRight
                            }

                        }
                    }

                    // ── regularisation, one shape per engine (§5.3) ───────────
                    GroupBox {
                        title: "Regularisation and engine settings"
                        Layout.fillWidth: true
                        Layout.margins: dp(6)

                        ColumnLayout {
                            anchors.fill: parent
                            spacing: dp(6)

                            // The four engines take regularisation in four incompatible forms, so
                            // this is a stack of four panels rather than one list with fields that
                            // would be meaningless five sixths of the time.
                            StackLayout {
                                id: regStack
                                Layout.fillWidth: true
                                currentIndex: root.regPageFor(root.engine)

                                // 0 — VMLMB: ["name", μ, extra...] specs
                                ColumnLayout {
                                    spacing: dp(6)
                                    Loader {
                                        Layout.fillWidth: true
                                        sourceComponent: specRegPanel
                                    }
                                    RowLayout {
                                        spacing: dp(8)
                                        Label { text: "maxiter" }
                                        SpinBox {
                                            from: 1; to: 100000; stepSize: 50; value: 200; editable: true
                                            onValueModified: root.vmlmbMaxiter = value
                                        }
                                    }
                                }

                                // 1 — BSMEM: entropy mode, and one prior image
                                ColumnLayout {
                                    spacing: dp(6)
                                    Label {
                                        Layout.fillWidth: true
                                        wrapMode: Text.WordWrap
                                        text: "MaxEnt ignores named regularisers: the prior image is the regularisation, " +
                                              "and the stopping criterion is the 4-element method vector."
                                        color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                    }
                                    RowLayout {
                                        spacing: dp(8)
                                        Label { text: "method preset" }
                                        ComboBox {
                                            id: bsmemPreset
                                            Layout.preferredWidth: dp(200)
                                            model: ["1 — classic known noise", "4 — χ² = N", "custom"]
                                            currentIndex: 1
                                            onActivated: {
                                                if (currentIndex === 0) root.bsmemMethod = [1, 1, 1, 2]
                                                else if (currentIndex === 1) root.bsmemMethod = [4, 1, 1, 2]
                                            }
                                        }
                                    }
                                    RowLayout {
                                        spacing: dp(4)
                                        Label { text: "method" }
                                        Repeater {
                                            model: 4
                                            delegate: SpinBox {
                                                id: methodBox
                                                required property int index
                                                from: 0; to: 9
                                                value: root.bsmemMethod[methodBox.index]
                                                implicitWidth: dp(80)
                                                enabled: bsmemPreset.currentIndex === 2
                                                onValueModified: {
                                                    var m = root.bsmemMethod.slice()
                                                    m[methodBox.index] = value
                                                    root.bsmemMethod = m
                                                }
                                            }
                                        }
                                    }
                                    RowLayout {
                                        spacing: dp(8)
                                        CheckBox {
                                            text: "prior image"
                                            checked: root.bsmemUsePrior
                                            onToggled: root.bsmemUsePrior = checked
                                        }
                                        Button {
                                            text: "Choose…"
                                            enabled: root.bsmemUsePrior
                                            onClicked: bsmemPriorDialog.openAt("")
                                        }
                                        Label {
                                            Layout.fillWidth: true
                                            text: root.bsmemPriorPath === ""
                                                  ? "x_start is used as the entropy prior"
                                                  : root.bsmemPriorPath
                                            color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                            elide: Text.ElideLeft
                                        }
                                    }
                                    RowLayout {
                                        spacing: dp(8)
                                        Label { text: "maxiter" }
                                        SpinBox {
                                            from: 1; to: 100000; stepSize: 50; value: 200; editable: true
                                            onValueModified: root.bsmemMaxiter = value
                                        }
                                    }
                                }

                                // 2 — BSDMM: proximal weights and operator names, not specs
                                ColumnLayout {
                                    spacing: dp(6)
                                    Label {
                                        Layout.fillWidth: true
                                        wrapMode: Text.WordWrap
                                        text: "ADMM takes proximal blocks, not name/μ specs. A weight of 0 removes its block."
                                        color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                    }
                                    GridLayout {
                                        columns: 3
                                        columnSpacing: dp(8)
                                        rowSpacing: dp(4)

                                        Label { text: "mu_tv" }
                                        NumField {
                                            Layout.preferredWidth: dp(90); minimum: 0
                                            value: root.bsdmmMuTv
                                            onCommitted: (v) => root.bsdmmMuTv = v
                                        }
                                        ComboBox {
                                            Layout.preferredWidth: dp(130)
                                            model: ["tv", "l2smooth"]
                                            onActivated: root.bsdmmRegType = currentText
                                        }

                                        Label { text: "mu_group" }
                                        NumField {
                                            Layout.preferredWidth: dp(90); minimum: 0
                                            value: root.bsdmmMuGroup
                                            onCommitted: (v) => root.bsdmmMuGroup = v
                                        }
                                        ComboBox {
                                            Layout.preferredWidth: dp(130)
                                            model: ["sparsity", "tv"]
                                            onActivated: root.bsdmmGroupType = currentText
                                        }

                                        Label { text: "mu_cen" }
                                        NumField {
                                            Layout.preferredWidth: dp(90); minimum: 0
                                            value: root.bsdmmMuCen
                                            onCommitted: (v) => root.bsdmmMuCen = v
                                        }
                                        Label {
                                            text: "centering"
                                            color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                        }
                                    }
                                    Label {
                                        Layout.fillWidth: true
                                        wrapMode: Text.WordWrap
                                        text: "mu_group couples the channels, so it does nothing on a monochromatic dataset."
                                        color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                    }
                                    RowLayout {
                                        spacing: dp(8)
                                        Label { text: "maxiter" }
                                        SpinBox {
                                            from: 1; to: 100000; stepSize: 50; value: 200; editable: true
                                            onValueModified: root.bsdmmMaxiter = value
                                        }
                                    }
                                }

                                // 3 — SQUEEZE and tempering: six names of their own
                                ColumnLayout {
                                    spacing: dp(6)
                                    Label {
                                        Layout.fillWidth: true
                                        wrapMode: Text.WordWrap
                                        text: "Squeeze regularisers act on the integer flux histogram: tv, entropy and " +
                                              "compactness are spelled like the gradient engines' but are different maths, " +
                                              "and l0 exists here only because no gradient is taken."
                                        color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                    }
                                    Repeater {
                                        model: squeezeRegModel
                                        delegate: RowLayout {
                                            id: squeezeRow
                                            required property int index
                                            required property var model
                                            Layout.fillWidth: true
                                            spacing: dp(6)
                                            CheckBox {
                                                Layout.preferredWidth: dp(140)
                                                text: squeezeRow.model.name
                                                checked: squeezeRow.model.use
                                                onToggled: squeezeRegModel.setProperty(squeezeRow.index, "use", checked)
                                            }
                                            Label { text: "λ" }
                                            NumField {
                                                Layout.preferredWidth: dp(80); minimum: 0
                                                enabled: squeezeRow.model.use
                                                value: squeezeRow.model.lambda
                                                onCommitted: (v) => squeezeRegModel.setProperty(squeezeRow.index, "lambda", v)
                                            }
                                            Button {
                                                text: root.squeezePriorPath === "" ? "image…" : "image ✓"
                                                visible: squeezeRow.model.name === "priorimage"
                                                enabled: squeezeRow.model.use
                                                onClicked: squeezePriorDialog.openAt("")
                                            }
                                            Label {
                                                visible: squeezeRow.model.name === "priorimage"
                                                Layout.fillWidth: true
                                                text: "zero pixels are a hard mask; supplying one disables auto-centering"
                                                color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                                wrapMode: Text.WordWrap
                                            }
                                            Item { Layout.fillWidth: true; visible: squeezeRow.model.name !== "priorimage" }
                                        }
                                    }

                                    GridLayout {
                                        Layout.fillWidth: true
                                        columns: 4
                                        columnSpacing: dp(8)
                                        rowSpacing: dp(4)

                                        Label { text: "nelements" }
                                        SpinBox {
                                            from: 0; to: 1000000; stepSize: 100; value: 0; editable: true
                                            onValueModified: root.nelements = value
                                            // wire: 0 means default_nelements(data)
                                        }
                                        Label { text: "niter (sweeps)" }
                                        SpinBox {
                                            from: 1; to: 1000000; stepSize: 100; value: 200; editable: true
                                            onValueModified: root.niter = value
                                        }

                                        Label { text: "nchains" }
                                        SpinBox {
                                            from: 1; to: 256; value: 4; editable: true
                                            onValueModified: root.nchains = value
                                        }
                                        Label { text: "tmin" }
                                        NumField {
                                            Layout.preferredWidth: dp(90); minimum: 0
                                            value: root.tmin
                                            onCommitted: (v) => root.tmin = v
                                        }

                                        Label { text: "f_anywhere" }
                                        NumField {
                                            Layout.preferredWidth: dp(90); minimum: 0; maximum: 1
                                            value: root.fAnywhere
                                            onCommitted: (v) => root.fAnywhere = v
                                        }
                                        Label { text: "f_copycat" }
                                        NumField {
                                            Layout.preferredWidth: dp(90); minimum: 0; maximum: 1
                                            value: root.fCopycat
                                            onCommitted: (v) => root.fCopycat = v
                                        }

                                        Label { text: "cent_mult" }
                                        NumField {
                                            Layout.preferredWidth: dp(90); minimum: 0
                                            value: root.centMult
                                            onCommitted: (v) => root.centMult = v
                                        }
                                        CheckBox {
                                            Layout.columnSpan: 2
                                            text: "auto-centering"
                                            checked: root.autoCentering
                                            enabled: root.squeezePriorPath === ""
                                            onToggled: root.autoCentering = checked
                                        }
                                    }

                                    // The ladder is Pigeons' half of the run, so it appears only when
                                    // tempering is the engine — the same panel otherwise describes a
                                    // single annealing schedule.
                                    RowLayout {
                                        Layout.fillWidth: true
                                        visible: root.engine === "tempering"
                                        spacing: dp(8)
                                        Label { text: "replicas" }
                                        SpinBox {
                                            from: 2; to: 512; value: 10; editable: true
                                            onValueModified: root.nreplicas = value
                                        }
                                        Label { text: "rounds" }
                                        SpinBox {
                                            from: 1; to: 64; value: 10; editable: true
                                            onValueModified: root.nrounds = value
                                        }
                                        Item { Layout.fillWidth: true }
                                    }

                                    Label {
                                        visible: root.squeezeUseSparco
                                        text: "SPARCO model (SqueezeSparco)"
                                        font.bold: true
                                    }
                                    Loader {
                                        Layout.fillWidth: true
                                        active: root.squeezeUseSparco
                                        visible: root.squeezeUseSparco
                                        sourceComponent: sparcoPanel
                                    }
                                }

                                // 4 — VI: no regularisers exist; the prior is the model
                                ColumnLayout {
                                    spacing: dp(6)
                                    Label {
                                        Layout.fillWidth: true
                                        wrapMode: Text.WordWrap
                                        text: "Variational inference has no regularisers: the correlated-field prior is the " +
                                              "regularisation, and the result is a posterior rather than an image."
                                        color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                    }
                                    RowLayout {
                                        spacing: dp(8)
                                        Label { text: "engine" }
                                        ComboBox {
                                            Layout.preferredWidth: dp(140)
                                            model: ["map", "mgvi", "geovi", "hybrid"]
                                            currentIndex: 1
                                            onActivated: root.viEngine = currentText
                                        }
                                        Label { text: "samples" }
                                        SpinBox {
                                            from: 1; to: 4096; value: 8; editable: true
                                            onValueModified: root.viSamples = value
                                        }
                                    }
                                    GridLayout {
                                        Layout.fillWidth: true
                                        columns: 5
                                        columnSpacing: dp(8)
                                        rowSpacing: dp(4)

                                        Label { text: "" }
                                        Label { text: "mean"; color: "#666" }
                                        Label { text: "std";  color: "#666" }
                                        Label { text: ""; Layout.columnSpan: 2 }

                                        Label { text: "slope" }
                                        NumField { Layout.preferredWidth: dp(80); value: root.cfSlopeMean
                                                   onCommitted: (v) => root.cfSlopeMean = v }
                                        NumField { Layout.preferredWidth: dp(80); minimum: 0; value: root.cfSlopeStd
                                                   onCommitted: (v) => root.cfSlopeStd = v }
                                        Label { Layout.columnSpan: 2; text: "power-law slope, normal prior"
                                                color: "#888"; font.pointSize: pt(baseFontPt - 2) }

                                        Label { text: "fluct" }
                                        NumField { Layout.preferredWidth: dp(80); minimum: 0; value: root.cfFluctMean
                                                   onCommitted: (v) => root.cfFluctMean = v }
                                        NumField { Layout.preferredWidth: dp(80); minimum: 0; value: root.cfFluctStd
                                                   onCommitted: (v) => root.cfFluctStd = v }
                                        Label { Layout.columnSpan: 2; text: "lognormal"
                                                color: "#888"; font.pointSize: pt(baseFontPt - 2) }

                                        Label { text: "flex" }
                                        NumField { Layout.preferredWidth: dp(80); minimum: 0; enabled: root.cfUseIwp
                                                   value: root.cfFlexMean; onCommitted: (v) => root.cfFlexMean = v }
                                        NumField { Layout.preferredWidth: dp(80); minimum: 0; enabled: root.cfUseIwp
                                                   value: root.cfFlexStd; onCommitted: (v) => root.cfFlexStd = v }
                                        CheckBox { Layout.columnSpan: 2; text: "integrated Wiener layer"
                                                   checked: root.cfUseIwp; onToggled: root.cfUseIwp = checked }

                                        Label { text: "asp" }
                                        NumField { Layout.preferredWidth: dp(80); minimum: 0; enabled: root.cfUseIwp
                                                   value: root.cfAspMean; onCommitted: (v) => root.cfAspMean = v }
                                        NumField { Layout.preferredWidth: dp(80); minimum: 0; enabled: root.cfUseIwp
                                                   value: root.cfAspStd; onCommitted: (v) => root.cfAspStd = v }
                                        Label { Layout.columnSpan: 2; text: "asperity, with the Wiener layer"
                                                color: "#888"; font.pointSize: pt(baseFontPt - 2) }

                                        Label { text: "offset" }
                                        NumField { Layout.preferredWidth: dp(80); enabled: root.cfUseOffset
                                                   value: root.cfOffsetMean; onCommitted: (v) => root.cfOffsetMean = v }
                                        NumField { Layout.preferredWidth: dp(80); minimum: 0; enabled: root.cfUseOffset
                                                   value: root.cfOffsetStd; onCommitted: (v) => root.cfOffsetStd = v }
                                        CheckBox { Layout.columnSpan: 2; text: "explicit DC mode"
                                                   checked: root.cfUseOffset; onToggled: root.cfUseOffset = checked }
                                    }
                                }

                                // 5 — SPARCO hybrid: the same specs as VMLMB, plus the parametric star
                                ColumnLayout {
                                    spacing: dp(6)
                                    Loader {
                                        Layout.fillWidth: true
                                        sourceComponent: specRegPanel
                                    }
                                    Label {
                                        Layout.fillWidth: true
                                        text: "image regularisation above; the star and background below are fitted parameters."
                                        wrapMode: Text.WordWrap
                                        color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                    }
                                    Loader {
                                        Layout.fillWidth: true
                                        sourceComponent: sparcoPanel
                                    }
                                    RowLayout {
                                        spacing: dp(8)
                                        Label { text: "rounds" }
                                        SpinBox {
                                            from: 1; to: 200; value: 5; editable: true
                                            onValueModified: root.sparcoRounds = value
                                        }
                                        Label {
                                            text: "image / parameters alternate"
                                            color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                        }
                                    }
                                }
                            }

                            // Switching engines cannot carry everything across, and silently dropping
                            // what does not translate is how a run ends up regularised differently
                            // from what the panel shows.
                            Label {
                                Layout.fillWidth: true
                                visible: text.length > 0
                                wrapMode: Text.WordWrap
                                text: ""
                                color: "#ef6c00"; font.pointSize: pt(baseFontPt - 2)
                                // wire: Julia reports what a change of engine could not translate
                                //       (μ weights carry VMLMB ↔ BSDMM; named specs do not reach BSMEM)
                            }
                        }
                    }
                }
            }

            // ── run (§5.4) ────────────────────────────────────────────────────
            Rectangle {
                Layout.fillWidth: true
                Layout.preferredHeight: dp(70)
                color: "#f4f4f4"
                border.color: "#ddd"

                // Two rows: the controls, then what they had to say. On one row the feedback
                // shared its width with four buttons and a progress bar, so anything longer
                // than "ready" was elided down to nothing exactly when it mattered.
                ColumnLayout {
                    anchors.fill: parent
                    anchors.margins: dp(6)
                    spacing: dp(4)

                    RowLayout {
                        Layout.fillWidth: true
                        spacing: dp(8)

                        Item {
                            implicitWidth: runButton.implicitWidth
                            implicitHeight: runButton.implicitHeight
                            Button {
                                id: runButton
                                text: "Run"
                                enabled: root.canRun
                                onClicked: root.startRun(false)
                            }
                            ReasonTip { reason: root.blockedReason }
                        }
                        // Continue rather than restart: x_start is the image already on screen, so
                        // more iterations can be added, or a regulariser introduced part-way,
                        // without going back to a generic start.
                        Button {
                            id: continueButton
                            text: "Run from previous"
                            enabled: root.canRun && root.hasResult
                            ToolTip.visible: hovered
                            ToolTip.text: root.hasResult
                                          ? "continue from the current image"
                                          : "nothing reconstructed yet"
                            onClicked: root.startRun(true)
                        }
                        Button {
                            id: stopButton
                            text: "Stop"
                            enabled: root.running
                            // wire: ask the running task to stop; SQUEEZE can return its best chain so far
                            onClicked: { }
                        }
                        // Nothing in the reconstruction holds the image still -- V² and closure
                        // phase are both translation-invariant -- so a result routinely sits off
                        // to one side of the field. This puts it back in the middle.
                        Button {
                            id: recenterButton
                            text: "Recenter"
                            enabled: root.hasResult && !root.running
                            ToolTip.visible: hovered
                            ToolTip.text: root.hasResult
                                          ? "shift the centroid to the centre of the field"
                                          : "nothing reconstructed yet"
                            onClicked: {
                                root.statusText = Julia.shell_recenter_image()
                                imageArea.update()
                                root.consoleChanged()
                            }
                        }
                        // Only while a run is going. Afterwards the space belongs to the result
                        // line, which otherwise overruns the bar and prints on top of it.
                        ProgressBar {
                            Layout.fillWidth: true
                            visible: root.running
                            from: 0; to: 1
                            value: Math.max(0, root.progress)
                            indeterminate: root.running && root.progress < 0
                        }
                    }

                    Label {
                        // Once a run has finished, show what it produced. "ready" is what the
                        // control said before the run and says nothing a user wants after it.
                        text: root.running
                              ? (root.roundText.length > 0 ? root.roundText : root.statusText)
                              : root.hasResult ? root.statusText
                              : (root.blockedReason.length > 0 ? root.blockedReason : "ready")
                        color: (!root.running && root.blockedReason.length > 0) ? "#c62828" : "#666"
                        elide: Text.ElideRight
                        Layout.fillWidth: true
                        ToolTip.visible: hovered && root.hasResult
                        ToolTip.text: root.statusText
                        HoverHandler { id: statusHover }
                        property bool hovered: statusHover.hovered
                    }
                }
            }
        }

        Rectangle { Layout.preferredWidth: 1; Layout.fillHeight: true; color: "#ddd" }

        // ── right: what came out ──────────────────────────────────────────────
        ColumnLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            Layout.margins: dp(6)
            spacing: dp(6)

            // ── result mode (§5.5) ────────────────────────────────────────────
            //
            // The switch is always present, even for the deterministic engines, because it is
            // the shape of the result that differs: annealing, tempering and VI hand back a
            // distribution, and a per-pixel σ map is the thing the other four cannot give.
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(8)

                Label { text: "Result:" }
                Item {
                    implicitWidth: modeRow.implicitWidth
                    implicitHeight: modeRow.implicitHeight

                    RowLayout {
                        id: modeRow
                        spacing: dp(8)
                        enabled: root.stochastic

                        RadioButton {
                            text: "mean image"
                            checked: root.resultMode === "mean"
                            // wire: push the posterior mean into the canvas
                            onClicked: root.resultMode = "mean"
                        }
                        RadioButton {
                            text: "uncertainty map"
                            checked: root.resultMode === "sigma"
                            // wire: push the per-pixel σ map into the canvas
                            onClicked: root.resultMode = "sigma"
                        }
                        RadioButton {
                            text: "sample"
                            checked: root.resultMode === "sample"
                            // wire: push sample sampleIndex into the canvas
                            onClicked: root.resultMode = "sample"
                        }
                    }
                    ReasonTip {
                        reason: root.stochastic ? ""
                              : "this engine returns one image, not a distribution — there is no mean, σ or sample to show"
                    }
                }

                Slider {
                    id: sampleSlider
                    Layout.preferredWidth: dp(180)
                    from: 0
                    to: Math.max(1, root.sampleCount - 1)
                    stepSize: 1
                    snapMode: Slider.SnapAlways
                    value: root.sampleIndex
                    enabled: root.stochastic && root.resultMode === "sample" && root.sampleCount > 1
                    // wire: show that sample of the ensemble
                    onMoved: root.sampleIndex = Math.round(value)
                }
                Label {
                    text: root.sampleCount > 0
                          ? (root.sampleIndex + 1) + " / " + root.sampleCount
                          : "no samples"
                    color: "#888"; font.pointSize: pt(baseFontPt - 2)
                }
                Item { Layout.fillWidth: true }
                Button {
                    text: "Save FITS…"
                    enabled: root.hasResult
                    // wire: write the current image (or the whole ensemble) to disk
                    onClicked: { }
                }
            }

            // ── canvas (§5.7) ─────────────────────────────────────────────────
            //
            // This perspective has its OWN figure, `imagePlot`, built beside the Exploring
            // one before the window exists. Sharing a single canvas would mean a
            // reconstruction wiping the plot being read on the other tab, and this tab's own
            // result appearing somewhere else. East is left and North is up, matching imdisp.
            Rectangle {
                id: imageCanvasMount
                Layout.fillWidth: true
                Layout.fillHeight: true
                Layout.minimumHeight: dp(240)
                color: "#f4f4f4"
                border.color: "#ddd"

                MakieArea {
                    id: imageArea
                    anchors.fill: parent
                    scene: imagePlot
                }

                // Over the canvas until there is something on it: an empty axis reads as a
                // broken reconstruction rather than as one that has not been run.
                Rectangle {
                    anchors.fill: parent
                    visible: !root.hasResult
                    color: "#f4f4f4"
                    Label {
                        anchors.centerIn: parent
                        horizontalAlignment: Text.AlignHCenter
                        color: "#888"
                        text: "no reconstruction yet"
                    }
                }
            }

            // ── colormap ──────────────────────────────────────────────────────
            //
            // Under the image because it describes what is on it. "hot" is `gist_heat`, which
            // is what `imdisp` uses, so a figure here and a figure from a script are the same
            // picture rather than merely similar ones. The two greys are for faint structure:
            // which of them reads better depends on whether it will end up on a white page.
            // A Flow rather than a row, so the names wrap onto a second line on a narrow
            // image panel instead of being pushed off the edge.
            Flow {
                Layout.fillWidth: true
                // A Flow's implicit height covers one row: it is computed before the layout
                // hands it a width, so it cannot know how many rows the buttons wrap onto.
                // childrenRect is measured after the arrangement, when the width is known.
                Layout.preferredHeight: childrenRect.height
                spacing: dp(4)
                enabled: root.hasResult

                Label {
                    text: "Colour map"
                    color: "#666"
                    font.pointSize: pt(baseFontPt - 1)
                    anchors.verticalCenter: undefined
                }
                Repeater {
                    model: root.colormapNames
                    Button {
                        required property string modelData
                        text: modelData
                        checkable: true
                        checked: root.colormap === modelData
                        onClicked: {
                            root.colormap = modelData
                            Julia.shell_image_colormap(modelData)
                            imageArea.update()
                        }
                    }
                }
            }

            // ── diagnostics (§5.5) ────────────────────────────────────────────
            //
            // Only the stochastic engines have these: "χ² stopped moving" is not convergence
            // for a sampler, and a tempering run fails in ways no χ² trace can show.
            GroupBox {
                title: root.engine === "tempering" ? "Tempering diagnostics"
                     : root.isSqueeze              ? "Annealing diagnostics"
                                                   : "Inference diagnostics"
                Layout.fillWidth: true
                visible: root.stochastic

                ColumnLayout {
                    anchors.fill: parent
                    spacing: dp(4)

                    Repeater {
                        model: root.diagRows

                        delegate: RowLayout {
                            id: diagRow
                            required property var modelData
                            Layout.fillWidth: true
                            spacing: dp(8)

                            Rectangle {
                                implicitWidth: dp(12); implicitHeight: dp(12)
                                radius: width / 2
                                color: root.lightColor(diagRow.modelData.state)
                                border.color: "#999"
                            }
                            Label {
                                Layout.preferredWidth: dp(120)
                                text: diagRow.modelData.label
                            }
                            Label {
                                Layout.fillWidth: true
                                text: diagRow.modelData.value
                                color: "#666"
                                font.family: "monospace"
                                font.pointSize: pt(baseFontPt - 2)
                                elide: Text.ElideRight
                            }
                            ReasonTip { reason: diagRow.modelData.tip }
                        }
                    }

                    // A single min(α) says a rung is blocked; it does not say which. The ladder
                    // does, which is the difference between a diagnosis and a complaint.
                    RowLayout {
                        Layout.fillWidth: true
                        visible: root.engine === "tempering"
                        spacing: dp(2)
                        Label {
                            Layout.preferredWidth: dp(120)
                            text: "ladder"
                        }
                        Repeater {
                            model: root.diagSwapPerRung
                            delegate: Rectangle {
                                id: rung
                                required property var modelData
                                Layout.fillWidth: true
                                Layout.preferredHeight: dp(14)
                                color: root.lightColor(rung.modelData <= 0 ? "bad"
                                                     : (rung.modelData < 0.3 ? "warn" : "ok"))
                                border.color: "#fff"
                                ReasonTip { reason: "rung swap acceptance " + root.fmt(rung.modelData, 2) }
                            }
                        }
                        Label {
                            visible: root.diagSwapPerRung.length === 0
                            Layout.fillWidth: true
                            text: "no rounds completed yet"
                            color: "#888"; font.pointSize: pt(baseFontPt - 2)
                        }
                        // wire: Julia pushes the Pigeons report — restarts, Λ, swap and explorer
                        //       acceptance, autocorrelation, log(Z₁/Z₀) — or SQUEEZE's
                        //       diagnostics NamedTuple, once per round
                    }
                }
            }

            // ── how well each observable is fitted ────────────────────────────
            //
            // This is where the L-curve used to be. A μ sweep reruns the whole reconstruction
            // per point and answers a question nobody has yet on a first pass; what a run
            // actually raises is "which observable is being fitted, and which is being given
            // up", and one aggregate chi2 cannot say.
            GroupBox {
                title: "Reduced χ² per observable"
                Layout.fillWidth: true

                ColumnLayout {
                    spacing: dp(4)

                    Repeater {
                        model: breakdownModel
                        delegate: RowLayout {
                            id: bdRow
                            required property var model
                            Layout.fillWidth: true
                            spacing: dp(8)

                            Label {
                                text: bdRow.model.obs
                                Layout.preferredWidth: dp(70)
                                font.bold: true
                                color: bdRow.model.used ? "#222" : "#999"
                            }
                            Label {
                                text: bdRow.model.chi2r
                                Layout.preferredWidth: dp(90)
                                horizontalAlignment: Text.AlignRight
                                // A reduced chi2 near 1 is the fit the error bars describe;
                                // far above means unfitted, far below means overfitted.
                                color: !bdRow.model.used ? "#999"
                                     : Math.abs(parseFloat(bdRow.model.chi2r) - 1) < 0.5 ? "#2e7d32"
                                     : parseFloat(bdRow.model.chi2r) > 5 ? "#c62828" : "#b06000"
                            }
                            Label {
                                text: bdRow.model.npts + " pts"
                                color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                Layout.preferredWidth: dp(80)
                            }
                            Label {
                                // An unticked observable was not fitted, so its chi2 is a
                                // prediction rather than a fit — worth saying, not hiding.
                                visible: !bdRow.model.used
                                text: "not fitted — predicted"
                                color: "#888"; font.pointSize: pt(baseFontPt - 2)
                            }
                            Item { Layout.fillWidth: true }
                        }
                    }

                    Label {
                        visible: breakdownModel.count === 0
                        text: "run a reconstruction to see this"
                        color: "#888"; font.pointSize: pt(baseFontPt - 2)
                    }
                }
            }
        }
    }
}
