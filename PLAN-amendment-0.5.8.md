# PLAN amendment 0.5.8: tbc-diag-v1 graded (D0-D3 pass), and a defect the tests did not look for

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-20. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

Sections 1-2 grade job `tbc-diag-v1` against amendment 0.5.7 section 4 (written before it ran).
Section 3 is a post-hoc observation, kept apart from the tests. Section 4 is new and is written
**before** job `tbc-fill-v1` produces any number. Nothing earlier changes.
New files: `results/tin/diag/` (the job's raw outputs), `tbc_fill.py` 0.1.0, `test_tbc_fill.py`.

## 1. What ran

`tbc-diag-v1` (tbc_diag 0.1.0, tbc_bte 0.1.0, tbc_anh 0.1.0; 306 s, all steps exit 0, unit tests
pass) reloaded the four orb-v3 npz files of `tbc-tin-v1` and re-solved TiN on diamond, Si and Ge
(001) at 300 K and at 900 / 700 / 700 K. The recovered cell volumes match a^3/4 of the relaxed
cells to 1e-4.

## 2. Grades (MW m^-2 K^-1)

| pair | T | G base | own DMM | G uni | uni1 (TiN gray) | uni2 (partner gray) | uni / DMM |
|---|---|---|---|---|---|---|---|
| TiN/C | 300 | 43.3 | 264.3 | 307.7 | 331.0 | 42.8 | 1.16 |
| TiN/C | 900 | 52.8 | 361.8 | 419.1 | 431.8 | 52.6 | 1.16 |
| TiN/Si | 300 | 52.3 | 315.3 | 431.8 | 398.5 | 52.1 | 1.37 |
| TiN/Si | 700 | 61.5 | 340.2 | 454.5 | 424.8 | 61.2 | 1.34 |
| TiN/Ge | 300 | 51.5 | 272.4 | 427.9 | 194.0 | 55.0 | 1.57 |
| TiN/Ge | 700 | 59.3 | 288.1 | 447.9 | 209.3 | 64.2 | 1.55 |

| | test | result | prior |
|---|---|---|---|
| D0 | closure (G and DMM within max(1, 2 %) of the printed values) | **PASS**, all 6 rows | 75 % |
| D1 | lengths x0.5 and x2 within 2 % | **PASS** (identical to 0.1) | 85 % |
| D2 | 800 nodes within 1 % | **PASS** (identical to 0.1) | 90 % |
| D3 | gray lengths restore DMM (main): G_uni >= G_dmm | **PASS**, all rows (ratios 1.16-1.57) | 70 % |
| D4 | attribution (information) | the TiN side carries the gap in every row: share of ln(G_uni/G_base) recovered by uni1 0.96-1.04 (C, Si) and 0.62-0.63 (Ge); by uni2 -0.01 to 0.04 | none |

`s1x0.01` (information only) equals base to 0.2, as amendment 0.5.9 section 4 predicted: scaling one
side's lengths cannot change G, so it cannot stand in for an electron bath.

**Reading by 0.5.7's rules.** D0-D2 and D3 pass, so the deficit below DMM is band non-equilibrium
inside the model, not a defect of the solver: with one length per side the same code gives the
dressed value above DMM. D4 names the TiN side in every row. 0.5.7's rule then says stage 3 adds an
electron bath in metals first. Section 3 shows why we do not apply that step yet.

## 3. Post-hoc observation (not a test): the TiN side's model conductivity is three times its input

The job printed each side's two-stream bulk conductivity, sum_b 2 K_b l_b. By construction
(`tbc_bte` docstring) it should equal the binned phono3py conductivity the lengths come from. For
TiN it does not:

| T | TiN model kappa | TiN binned phono3py kappa | partner model kappa (C / Si / Ge) | partner phono3py |
|---|---|---|---|---|
| 300 K | **132.3** | 44.6 | 881.5 / 62.9 / 20.6 | 780.8 / 61.5 / 20.5 |

The cause is in `tbc_bte.lengths_from_kappa`. A bin with a radiation conductance K_b > 0 on the dense
harmonic mesh but no mode on the 16^3 lifetime mesh is "filled" with a length log-interpolated from
its neighbours, and below the first sampled bin `np.interp` copies the first sampled length. TiN has
11 such bins (diamond 28, Si 6, Ge 2). Since the sampled bins hold at most 44.6 W/m/K, the filled
TiN bins carry at least (132.3 - 44.6) / 132.3 = **66 %** of the model's TiN conductivity, with
lengths no phono3py mode produced. The diagnostic's own numbers point the same way: TiN's
conductivity-weighted mean length is 1,390 nm over all bands but 106 nm over the bands diamond can
receive, and only 26 % of TiN's model conductivity sits in bands diamond shares (Si 78 %, Ge 99 %).
We expect, but have not shown, that the filled bins are TiN's lowest-frequency bins (below about
2.2 THz, under the lowest phonon of the 16^3 mesh), where diamond's dense mesh also has no modes.

A rough estimate shows this can produce the whole deficit. If two thirds of the heat in TiN rides in
bands that cannot cross and must first hand their heat to other bands over about 1.4 um, the extra
resistance is of order (2/3)^2 x 1.4 um / (44.6 W/m/K) = 1.4e-8 m^2 K/W, i.e. a ceiling near
70 MW m^-2 K^-1, close to the 43-61 the model gives. **So the 5-6x deficit of 0.5.7 section 3 may be
an artefact of the frequency discretisation, not physics.** D3 could not see this: gray lengths
discard exactly the band structure the filled bins distort. We missed it when 0.5.7 was written
because we checked the binned sums (44.6 = 44.6) but not the model's own sum.

Consequences, before any new number:
- 0.5.7's T-grades stand as graded: they grade the model as it was run.
- 0.5.7 section 3's reading ("the multiband structure" causes the deficit; "TiN's missing electron
  bath" as the first suspect) and this amendment's D-reading are **suspended** until `tbc-fill-v1`
  is graded. Stage 3c (electron bath) is not designed until then.
- Every earlier pair (Al/Si, Al/GaN, ZnO/GaN, AlN/GaN) used the same filling rule. Their milder
  deficits below DMM (ZnO/GaN 190 against 314) may share the cause. They are re-checked after
  section 4 is graded.

## 4. Job `tbc-fill-v1` and its pre-registered grades (written before it ran)

`tbc_fill.py` 0.1.0 reloads the same npz files and re-solves each pair with (a) `base`, the job's
lengths; (b) `merge`: every bin from 0 up to the first bin sampled on **both** sides' lifetime meshes
is merged into one band per side (K and kappa summed), so every low-frequency length is a real
kappa / (2 K) and each side's bulk conductivity is conserved, with the band-summed DMM recomputed on
the merged bands (`dmm_m`); (c) `nofill` (information): filled bins get a length of 1e-12 m, so they
carry no bulk heat and stay at the bath temperature; (d) `uni_m` (information): merged bands with
one length per side; (e) merge at 800 nodes (information). It prints the lowest 30 bins of each
material at 300 K (K, kappa, length, filled flag).

| | test | PASS if | prior |
|---|---|---|---|
| F0 | closure | base within max(0.5, 1 %) of tbc-diag-v1's base, all 6 rows | 95 % |
| F1 | the filled bins are the low-frequency gap | TiN at 300 K: filled share of model conductivity >= 0.60 AND every filled bin starts below 3.0 THz | 80 % |
| F2 | **merging removes most of the deficit (main)** | G_merge >= 2 x G_base, all 6 rows | 65 % |
| F3 | the rest is in the literature's range | G_merge x 1.6 >= DMM_merge, all 6 rows (Han and Lee, arXiv:2309.13187: non-equilibrium raises resistance 2-3x over Landauer; 1.6 is the factor we saw for ZnO/GaN and AlN/GaN) | 45 % |
| F4 | consistency | TiN model conductivity after merging within 5 % of the binned phono3py value, all 6 rows | 85 % |

Reported, never graded: `nofill`, `uni_m`, merge at 800 nodes, the merged range, and merge + the
printed G_inel against Khan et al.'s measured values.

**Reading rules.**
- F0 FAIL: the reconstruction differs from tbc-diag-v1; F1-F4 are indicative only.
- F1 PASS and F2 PASS: the 5-6x deficit of 0.5.7 section 3 is mostly a discretisation artefact.
  We say so in the README and in any article on TiN, retire the "below its own DMM" reading, make a
  conserving length rule the default in a NEW module (tbc_bte stays as published), and re-solve
  every earlier pair with it before stage 3c is designed.
- F2 PASS, F1 FAIL: merging helps but the filled bins are not only the low-frequency gap; the
  printed bin table decides where they are before any rule is frozen.
- F2 FAIL: filling is not the main cause; 0.5.7's reading and the electron-bath step (stage 3c)
  come back, with the filled-bin inflation listed as a known secondary defect.
- F3 tells whether what remains after merging is in the range other Boltzmann studies report; F3
  FAIL with F2 PASS means a residual non-equilibrium larger than published, to be explained.
- F4 FAIL: interior filled bins (not the low-frequency gap) still inflate the model; fix those
  before using merge as the default.

## 5. Limits

All of section 3 is post-hoc and rests on printed sums; only `tbc-fill-v1` can confirm where the
filled bins are. The rough estimate in section 3 is an order of magnitude, not a fit. The measured
values are blind figure readings (+-5-30). The universal potential under-predicts the bulk
conductivity of Si, Ge and diamond 2-3x (0.5.7 section 6). The no-fit model is still not a
candidate Grand Challenge entry.
