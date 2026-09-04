# Variational inference for optical interferometry.
#
# This was the OIVI.jl package. It is here rather than beside OITOOLS for one reason: it is
# built ON OIdata, setup_ft and the chi2 kernels, so as a separate package it could only ever
# be DOWNSTREAM of OITOOLS — which is how it drifted. Its manifest pinned OITOOLS 0.9.5 against
# 0.13.x, and the breakage that caused (Float32 visibilities against Float64 signatures) sat
# unnoticed until someone ran it. In one repository the operators and the data structures they
# read change in the same commit and are tested by the same run.
#
# VarInf stays a separate package on purpose: it holds the generic VI core — the
# AbstractInferenceProblem protocol, NIFTy-style CG / Newton-CG, correlated fields, and the
# generic reconstruct_* — and nothing interferometric. It is reusable for other domains, so
# absorbing it would be the wrong direction.
module OITOOLSVarInfExt

using FFTW
using LinearAlgebra
using Printf
using Statistics

using OITOOLS
using OITOOLS: OIdata, FlatModel, NFFTCell, DFTCell, mod360, ft_eltype,
               scatter_obs_cotangent!, _pad_weights, _chi2_terms

using VarInf
using VarInf: _nifty_cg, _nifty_newton_cg, _spectral_latent_size
import VarInf: AbstractInferenceProblem,
               energy_and_gradient, latent_size, transformation,
               right_sqrt_metric, left_sqrt_metric, data_size,
               reconstruct_map, reconstruct_mgvi,
               reconstruct_geovi, reconstruct_hybrid

const VIDIR = joinpath(dirname(@__DIR__), "src", "vi")

include(joinpath(VIDIR, "sky.jl"))          # correlated-field sky + limb-darkened disc weight
include(joinpath(VIDIR, "diffphase.jl"))    # differential phase, cross-channel
include(joinpath(VIDIR, "observe.jl"))      # image -> visibilities -> observables, + JVPs
# `pointsource.jl` before `problem.jl`: `PointSourceProblem` takes a `PointSourceParams` as a
# FIELD, so that type must exist when the struct is defined. The reverse reference — point
# sources building a `PointSourceProblem` — happens inside function bodies, which resolve late.
include(joinpath(VIDIR, "pointsource.jl"))  # analytic point sources and their geoVI problem
include(joinpath(VIDIR, "problem.jl"))      # likelihood energy + the VarInf protocol wrappers
include(joinpath(VIDIR, "reconstruct.jl"))  # MAP / MGVI / geoVI / hybrid entry points

export harmonic_smooth, make_smoothing_kernel
export CorrFieldConfig, FourierGridInfo, fourier_mode_distributor, amplitude_spectrum
export limb_weight, SkyModelParams, sky_forward, sky_adjoint, frozen_spectral_range,
       frozen_chromatic_range
export report_chi2, minisanity
export DiffPhaseConfig, ObservationConfig, ObsContext, observe, observe_adjoint, observe_data
export build_diffphase_config, diffphi_forward, diffphi_adjoint, diffphi_jvp
export sky_forward_jvp, observe_jvp
export energy_fg, energy_fg!, total_latent_size, n_model_params
export reconstruct_map, reconstruct_mgvi, reconstruct_geovi, reconstruct_hybrid
export PointSourceParams, ps_latent_size, ps_unpack, ps_energy_fg
export PointSourcePosterior, reconstruct_pointsource, ps_report_chi2
export AbstractInferenceProblem, InterferometricProblem, PointSourceProblem
export energy_and_gradient, latent_size, transformation
export right_sqrt_metric, left_sqrt_metric, data_size

end # module
