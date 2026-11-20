# PLAN amendment 0.5.11: tbc-fill-v1 graded: TiN's 5-6x deficit below its own DMM was mostly a binning artefact (F2 passed), but the unsampled bins reach 3.2 THz and still inflate TiN where merging stops early (F1, F4 failed); next, a conserving "cover" band rule (tests V0-V5)

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-20. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

Sections 1-3 grade job `tbc-fill-v1` against amendment 0.5.8 section 4 (written before it ran;
thresholds unchanged). Section 4 is written **before** any code for the next stage exists.
Raw output, unedited: `results/tin/fill/` (`run.txt`, `fill-summary.json`, `stdout-fill.txt`,
`tests.txt`, `versions.txt`, `sha256.txt`, `compile.txt`, `fetch.txt`, `job.json`) and our grades
`results/tin/fill/grades-fill-v1.json`. `tbc-anh-v2` is still queued; no module was changed.

## 1. The run

python:3.12, `tbc_fill.py` 0.1.0 with `tbc_bte`, `tbc_anh`, `tbc_transport` 0.1.0 fetched from
main; the four orb-v3 npz files of `tbc-tin-v1` reloaded. All steps exit 0, 263 s; unit tests
pass. Recovered cell volumes equal a^3/4 to 1e-4. TiN's lowest positive frequency on the 16^3
lifetime mesh is 0.568 THz (diamond 3.705, Si 0.848, Ge 0.440).

## 2. Grades (MW m^-2 K^-1)

| pair | T (K) | G base | G merge | merge / base | DMM (merged bands) | merge / DMM | nofill (info) |
|---|---|---|---|---|---|---|---|
| TiN/C | 300 | 43.3 | 192.9 | 4.46 | 266.4 | 0.72 | 56.6 |
| TiN/C | 900 | 52.8 | 235.8 | 4.46 | 363.9 | 0.65 | 80.8 |
| TiN/Si | 300 | 52.3 | 317.0 | 6.06 | 315.5 | 1.00 | 85.5 |
| TiN/Si | 700 | 61.5 | 357.4 | 5.81 | 340.3 | 1.05 | 113.4 |
| TiN/Ge | 300 | 51.5 | 169.6 | 3.29 | 272.6 | 0.62 | 70.2 |
| TiN/Ge | 700 | 59.3 | 188.0 | 3.17 | 288.3 | 0.65 | 89.4 |

| | test | result | prior | grade |
|---|---|---|---|---|
| F0 | closure: base equals tbc-diag-v1's, all rows | all 6 rows | 95 % | **PASS** |
| F1 | TiN 300 K: filled share >= 0.60 AND every filled bin below 3.0 THz | share 0.663; highest filled bin 3.2 THz | 80 % | **FAIL** |
| F2 | **main**: G_merge >= 2 x G_base, all rows | 3.17-6.06 | 65 % | **PASS** |
| F3 | G_merge x 1.6 >= DMM_merge, all rows | 5 of 6; TiN/Ge 300 K 271.4 against 272.6 | 45 % | **FAIL** |
| F4 | TiN model conductivity after merging within 5 % of phono3py, all rows | TiN/C exact; TiN/Si +13 % / +11 %; TiN/Ge +22 % / +18 % | 85 % | **FAIL** |

Expectations met: F0 and F2 (2 of the 3 written with a prior above 50 %); F1 and F4 were missed.
F3's 45 % was a minority view; it failed by 0.4 % in one row.

## 3. Reading (0.5.8 sec. 4 rules, written before the job), applied

- **F2 PASS, F1 FAIL**: *merging helps but the filled bins are not only the low-frequency gap; the
  printed bin table decides where they are before any rule is frozen.* The table (run.txt): TiN's
  filled bins at 300 K start at 0.0, 0.2, 0.6, 0.8, 1.0, 1.2, 1.4, 2.0, 2.2, 2.4 and 3.2 THz,
  interleaved with sampled bins at 0.4, 1.6, 1.8, 2.6, 2.8 and 3.0 THz. So they are the whole
  acoustic range below 3.4 THz that a 16^3 mesh samples only sparsely, not only the gap below its
  lowest phonon. Diamond has 28 filled bins up to 10.4 THz, but they carry only 12 % of its model
  conductivity; Si 6 bins (2 %) and Ge 2 (0.2 %).
- **F4 FAIL**: *interior filled bins still inflate the model; fix those before using merge as the
  default.* Merging stops where both sides are first sampled (4.8 THz with diamond, 1.0 with Si,
  0.6 with Ge), so with Si and Ge, TiN's filled bins at 1.2-3.2 THz keep invented lengths.
- **What is retired and what is not.** The reading of amendment 0.5.7 sec. 3 ("TiN interfaces sit
  5-6x below their own DMM because of multiband non-equilibrium; the electron bath is the first
  suspect") is **retired**: most of that deficit came from the filling rule. G rises 3.2-6.1x when
  low-frequency lengths are made from real conductivities. What is left is not yet a physical
  number: TiN/Si sits at its DMM and TiN/C and TiN/Ge at 0.62-0.72 of it, but in the Si and Ge rows
  part of TiN is still inflated (F4). Stage 3c (electron bath) waits for section 4.
- Information, not graded: merge plus the printed inelastic term against Khan et al. (2024):
  TiN/C 199 / 449 (300 K) and 280 / 910 (900 K); TiN/Si 373 / 507 and 530 / 916; TiN/Ge 231 / 190
  and 376 / 281. The no-fit model is now below the diamond and Si measurements and above Ge.
- Every earlier pair (Al/Si, Al/GaN, ZnO/GaN, AlN/GaN) used the same filling rule and is re-solved
  under section 4's rule before any number from it is used again.

## 4. Next stage, pre-registered now: the "cover" band rule (job `tbc-cover-v1`, tests V0-V5)

**Rule.** A new module `tbc_cover.py` 0.1.0 (tbc_bte, tbc_fill and every published module stay as
they are). Starting from 0 THz, consecutive fine bins (0.2 THz) are grouped into one band until the
group contains at least one bin sampled on the lifetime mesh **on both sides** of the interface;
the band then closes and the next starts. A tail left without a sampled bin on both sides joins the
last band. In each band, K (radiation conductance) and kappa are summed per side, the length is
kappa_band / (2 K_band), and the band-summed DMM is recomputed on the same bands. No length is
interpolated, so each side's model conductivity equals its binned phono3py value by construction.
The first band is the `merge` band of tbc-fill-v1, so the rule extends it to every band. It is
applied to the same six TiN rows and, where the npz files of earlier stages are published, to
Al/Si, Al/GaN, ZnO/GaN and AlN/GaN at 300 K (orb-v3 where available, otherwise the potential the
earlier stage used; reported as such).

| | test | PASS if | prior |
|---|---|---|---|
| V0 | gate: unit tests pass; the first cover band equals tbc-fill-v1's merged range in all 6 TiN rows; G_base and G_merge reproduce tbc-fill-v1 within max(0.5, 1 %) | all hold | 85 % |
| V1 | conservation: each side's model conductivity within 1 % of its binned phono3py value, every row | all rows | 90 % |
| V2 | **main**: G_cover within 25 % of G_merge in the TiN/C rows (where merge already covers all of TiN's filled bins) and G_cover >= G_merge in the TiN/Si and TiN/Ge rows (the invented long lengths that remain there depress G, as they did in base) | all 6 rows | 55 % |
| V3 | G_cover / DMM_cover in [0.4, 1.4], all TiN rows (below 1: multiband non-equilibrium as reported by Han and Lee, arXiv:2309.13187; up to 1.37: our own gray-length value of 0.5.8) | all 6 | 70 % |
| V4 | coarse-banding check: G_cover with bands at least 0.4 THz wide (two fine bins minimum) within 10 % of G_cover, all TiN rows | all 6 | 60 % |
| V5 | earlier pairs: G_cover / G_base >= 1.2 for ZnO/GaN (its deficit below DMM was the largest of them, 190 against 314) | holds | 50 % |

Reported, never graded: band count and widths per row; the cover result plus the printed
inelastic term against the measured values; the four earlier pairs' G_base, G_cover and DMM.

**Reading rules, written now.**
- V0 FAIL: the reconstruction differs; V1-V5 indicative only.
- V1 PASS and V2 PASS: the cover rule becomes the default length rule for every new stage, and
  stage 3c (electron bath in TiN, with orb-v3) is designed on it.
- V2 FAIL: the band boundaries, not only the filling, move G (or the interior lengths were not
  depressing G as we think); we then need a denser lifetime mesh
  (TiN at 24^3 or more) before trusting any band rule, and say so.
- V3 FAIL: what remains below (or above) DMM after the artefact is removed is larger than the
  range other Boltzmann studies report and must be explained before stage 3c.
- V5 informs whether the earlier pairs' milder deficits share the artefact.

## 5. Limits

The universal potential under-predicts the bulk conductivity of Si, Ge and diamond 2-3x (0.5.7
sec. 6). The measured values are blind figure readings (+-5-30). The cover rule trades the
invented lengths for coarser bands; V4 measures how much the coarseness matters, and a denser
lifetime mesh is the physical fix. The no-fit model is still not a candidate Grand Challenge entry.
