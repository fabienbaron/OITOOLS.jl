# Demo data

The OIFITS files that ship with OITOOLS, and what each one is good for. The numbers are what
`readoifits` reports for the first epoch, so they say what a file can actually exercise — a
dataset with `vis=0` cannot show you a differential-phase plot however the panel is set.

## Single files

| file | target | V² | T3φ | VIS | λ channels | band | tel |
|---|---|---|---|---|---|---|---|
| `AX_Cir.oifits` | AX Cir | 900 | 600 | 0 | 3 | H | 4 |
| `AZ_Cyg_2011.oifits` | AZ Cyg | 916 | 1199 | 0 | 8 | H | 6 |
| `AZ_Cyg_2014.oifits` | AZ Cyg | 600 | 586 | 0 | 16 | H | 6 |
| `AlphaCenA.oifits` | α Cen A | 324 | 156 | 0 | 12 | H | 7 |
| `AlphaCenB.oifits` | α Cen B | 432 | 226 | 0 | 12 | H | 7 |
| `Bet_Lyr6T.oifits` | β Lyr | 3836 | 4636 | 0 | 8 | H | 6 |
| `HD140573.oifits` | BSC3482 | 1339 | 609 | 0 | 28 | 0.56–0.86 µm | 3 |
| `Iota_Peg4T.oifits` | ι Peg | 96 | 64 | 0 | 8 | H | 4 |
| `Iota_Peg6T.oifits` | ι Peg | 636 | 800 | 0 | 8 | H | 6 |
| `MWC275_4T.oifits` | MWC 275 | 522 | 402 | 4024 | 18 | H | 10 |
| `MWC480.oifits` | HD 31648 | 1233 | 1511 | 0 | 8 | H | 6 |
| `Pi1_Gru.oifits` | π¹ Gru | 909 | 603 | 0 | 6 | H | 7 |
| `Rho_Cas.oifits` | ρ Cas | 1843 | 2460 | 26445 | 41 | K | 6 |
| `V1295_Aql.oifits` | HD 190073 | 955 | 840 | 2225 | 5 | H | 15 |
| `polaris.oifits` | HD 8890 | 3194 | 1928 | 50524 | 98 | H | 6 |

Picking one:

- **`V1295_Aql.oifits` is the workhorse.** Polychromatic without being huge, and what most of the
  imaging and model-fitting examples and tests are written against. Start here.
- **`AlphaCenA.oifits` for a limb-darkened disc** — a resolved single star, which is what the
  limb-darkening laws in the Model perspective are for.
- **`polaris.oifits` and `Rho_Cas.oifits` are the big ones**, 98 and 41 channels. Use them to see
  how a view behaves with real spectral coverage, and expect reconstructions to take a while.
- **`HD140573.oifits` is the only visible-band file** (0.56–0.86 µm); everything else is H or K.
- **`Beauty_Contest_2004/2004-data1.oifits` is monochromatic**, which makes it the degenerate
  case for anything that colours by wavelength.

## Directories

| directory | what is in it |
|---|---|
| `Beauty_Contest_2004/` | The 2004 imaging beauty contest: `2004-data1.oifits`, small, monochromatic and well understood. Its truth image is `images/`. |
| `Beauty_Contest_2026/` | The 2026 contest: `OBJECT1_LM`, `OBJECT1_N`, `OBJECT2_K`. The only shipped files carrying **differential phase and OI_FLUX**, so the only ones on which the diffphi and flux views show anything. |
| `images/` | FITS images rather than interferometric data: truth images for the beauty contest. Kept apart so the file picker's OIFITS listing is data only. |

## A note on the names

Several files were renamed in 0.15: `2019_v1295Aql.WL_SMOOTH.A.oifits` → `V1295_Aql.oifits`,
`AZCYG2011_FINAL2018` → `AZ_Cyg_2011`, `AZCYG2014NEW` → `AZ_Cyg_2014`, `betlyr6t` → `Bet_Lyr6T`,
`MWC275_T4a` → `MWC275_4T`, `iota_peg4t`/`iota_peg6t` → `Iota_Peg4T`/`Iota_Peg6T`,
`rho_Cas_example` → `Rho_Cas`, `pigru` → `Pi1_Gru`, `AXCir` → `AX_Cir`. A script of your own that
refers to the old names needs updating; everything in this repository was swept.

### `o_Leo/` and `HD210763/` — files that break the standard

One file from each of PMOIRED's o Leo (EX5) and HD 210763 (EX6) examples, kept because they are
the only examples of these defects, which `load_oifits` reads around (`PLAN_C14.md`):

- **`o_Leo/GRAVI.2022-02-28…_singlesciviscalibrated.fits`** — GRAVITY, SC and FT. Its OI_ARRAY
  declares revision 2 without the FOV column; its OI_FLUX tables name their data `FLUX`, not
  `FLUXDATA`; a `TELLURICS` HDU. ~9750 V², ~6530 T3, ~9790 VIS, 6536 flux.
- **`HD210763/2022-08-03-calibrated_cleaned.fits`** — a **ZIP archive** holding one FITS file
  (written by OICLEAN). Revisions 5 and 3, which do not exist; OI_VIS without VISPHI; half the rows
  with TARGET_ID 0 (not in OI_TARGET) and OI_FLUX with meaningless TARGET_IDs, so those rows and
  that table are left out. 3588 V², 1651 T3, 3771 VIS amplitudes.
