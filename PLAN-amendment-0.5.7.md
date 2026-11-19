# PLAN amendment 0.5.7: TiN interfaces to 900 K graded (stage 3b fails), and why the model falls below its own DMM

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-19. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

Sections 1-3 grade job `tbc-tin-v1` against amendment 0.5.2 (written before it ran). Section 4 is
new and is written **before** job `tbc-diag-v1` produces any number. Nothing earlier changes.
New files: `results/tin/orb-v3/` (the job's raw outputs, including the four per-material npz
files), `results/tin/grades-orb-v3.json`, `tbc_diag.py` 0.1.0, `test_tbc_diag.py`, `fetch_npz.py`.

## 1. What ran

`tbc-tin-v1` computed the full orb-v3 (primary) half: three-phonon lifetimes of TiN, diamond, Si
and Ge (300 K conductivities 44.6, 780.8, 61.5 and 20.5 W m^-1 K^-1), the dressed elastic G_BTE,
the no-fit three-phonon channel G_inel and G_total for TiN/C, TiN/Si and TiN/Ge (001), 100-1000 K.
The MACE-MPA-0 cross-check was killed by the 3 GB memory cap while relaxing TiN, so the summary
file and the automatic grades were never written. The orb-v3 numbers were printed line by line
(`results/tin/orb-v3/stdout-tin.txt`, lines `G orb-v3 ...`), and we graded them by hand with the
rules of `tbc_tin.grade_t`. Every grade below holds for any value within +-0.5 of the printed
integers.

## 2. Grades (orb-v3; MW m^-2 K^-1)

| pair | T | G_BTE | G_inel | G_total | own DMM | measured |
|---|---|---|---|---|---|---|
| TiN/C | 300 | 43 | 5.7 | 49 | 264 | 449 |
| TiN/C | 900 | 53 | 44.1 | 97 | 362 | 910 |
| TiN/Si | 300 | 52 | 56.4 | 109 | 315 | 507 |
| TiN/Si | 700 | 61 | 172.4 | 234 | 340 | 916 |
| TiN/Ge | 300 | 51 | 61.8 | 113 | 272 | 190 |
| TiN/Ge | 700 | 59 | 188.2 | 248 | 288 | 281 |

| | test | result | prior |
|---|---|---|---|
| T0 | sanity | **INCOMPLETE**: all pairs computed; the frequency and imaginary-mode checks were only in the summary the killed job never wrote | 80 % |
| T1 | elastic part rises less than data | **PASS** (1.233 vs 2.027, 1.173 vs 1.807, 1.157 vs 1.479) | 85 % |
| T2 | level within 25 % | **FAIL** 0/3 (C -89 %/-89 %, Si -79 %/-74 %, Ge -41 %/-12 %) | 30 % |
| T3 | rise within 15 % (main) | **FAIL** 1/3 (C 1.980 vs 2.027 passes; Si 2.147 vs 1.807; Ge 2.195 vs 1.479) | 30 % |
| T4 | closer to data than own DMM | **FAIL** (813 vs 548, 682 vs 576, 33 vs 7) | 70 % |
| T5 | ordering Si > C > Ge | **FAIL** (109, 49, 113) | 50 % |
| T6 | Ge/Si at 300 K in [0.25, 0.55] | **FAIL** (1.037) | 55 % |

**Reading (by amendment 0.5.2's rules).** T1 passes: the model's elastic part cannot produce the
measured rise, as Khan et al. concluded. T3 fails, but not in the direction 0.5.2's rule
anticipated ("rising too little"). For TiN/diamond the rise is right (1.98 against 2.03), yet at a
level nine times too low. For TiN/Si and TiN/Ge the model rises too fast, because its three-phonon
channel grows about threefold from 300 to 700 K and dominates the total. As 0.5.2 asked, we report
"right trend, wrong level" (TiN/C) and "wrong trend" (Si, Ge) as they are. **The no-fit model as it
stands is not a candidate entry for TiN interfaces.** The priors for T4-T6 (50-70 %) were
over-optimistic, and we say so.

## 3. The defect that matters most

The model's elastic part lies **far below its own DMM**: 43 against 264 (TiN/C), 52 against 315
(TiN/Si), 51 against 272 (TiN/Ge), factors of 5.3-6.1. The same happened, more mildly, for ZnO/GaN
(190 against 314, amendment 0.5.6) and AlN/GaN. In the one-band limit the solver reproduces the
Simons / Zeng-Chen dressed conductance, which always lies **above** DMM. So the deficit comes from
the multiband structure: different bands relax to the local bath over very different lengths, and
heat entering through a few bands has to be handed to the rest.

This direction is known. Han and Lee ([arXiv:2309.13187](https://arxiv.org/abs/2309.13187)) solved
the Peierls-Boltzmann equation with DMM transmission and ab initio three-phonon scattering for 36
III-V interfaces. Non-equilibrium phonons made the total resistance two to three times the Landauer
value, and more so for larger Debye-temperature mismatch. Feng, Zhong, Shi and Ruan (PRB 99,
045301, 2019) also find strong modal non-equilibrium near interfaces. Our factors of about 1.6
(ZnO/GaN, AlN/GaN) fall in Han and Lee's range. The factors of 5-6 for TiN on much stiffer or
softer partners are larger, and they could be either real or a defect.

**Physics the model leaves out.** TiN is a metal. Its phonons scatter on electrons
(g ~ 1e18 W m^-3 K^-1 per Khan et al.), which hold every TiN phonon band close to the electron
temperature within nanometres. Our TiN side has only phonon-phonon lifetimes, so it can sit far
from equilibrium. If the TiN side causes most of the deficit, the missing electron bath is the
first thing stage 3 must add. If the partner side causes it, the deficit is Han and Lee's effect
in diamond, Si or Ge, and the measured values (which exceed even our DMM for TiN/C and TiN/Si) then
need a channel we do not have.

## 4. Job `tbc-diag-v1` and its pre-registered grades (written before it ran)

`tbc_diag.py` reloads the four npz files, recovers the cell volume from the stored phono3py
conductivity (checked against a^3/4), rebuilds the exact band conductances and two-stream lengths
of the job (oriented (001), isotopes included), and re-solves each pair at 300 K and T_high with
the unchanged `tbc_bte.solve_two_stream`. Variants: `uni` (every band on both sides gets the same
length, keeping each side's bulk conductivity and the per-band DMM transmission), `uni1` (TiN side
only), `uni2` (partner side only), lengths x0.5 and x2, 800 nodes, and, for information only,
`s1x0.01` (TiN lengths divided by 100, a crude stand-in for electron-pinned TiN phonons).

| | test | PASS if | prior |
|---|---|---|---|
| D0 | closure | recomputed G_BTE and G_dmm within max(1, 2 %) of the printed values, all 6 rows | 75 % |
| D1 | scale invariance | x0.5 and x2 within 2 % of base, all rows | 85 % |
| D2 | grid | 800 nodes within 1 % of base, all rows | 90 % |
| D3 | **gray lengths restore DMM (main)** | G_uni >= G_dmm, all rows | 70 % |
| D4 | attribution (information) | for each row, the share of the log-gap ln(G_uni/G_base) recovered by uni1 and by uni2; the larger names the side | none |

**Reading rules.**
- D0-D2 PASS and D3 PASS: the deficit below DMM is band non-equilibrium inside the model, not a
  solver defect. Stage 3 must then decide, with the literature above, whether that much
  non-equilibrium is physical. If D4 names the TiN side, stage 3 adds an electron bath in metals
  (electron-phonon lifetimes in the TiN lengths) before anything else.
- D3 FAIL with D0-D2 PASS: even with equal lengths the solver gives less than DMM, so the per-band
  coupling or the solver has a defect. Stage 3 does not freeze until it is found and fixed.
- D1 or D2 FAIL: numerical defect; fix first, then re-run.
- D0 FAIL: the reconstruction does not reproduce the job. D1-D4 are then read as indicative only.
- `s1x0.01` is reported, never graded.

## 5. Other consequences

- MACE-MPA-0 for TiN is not re-run now. The defect in section 3 sits in the transport model, not
  in the potential, so a second potential would not change what we do next.
- T0 stays incomplete. `tbc-diag-v1` prints the top frequency and imaginary share on the
  lifetime mesh (16^3), which is coarser than the mesh T0 named, as information only.
- The Grand Challenge entry stays far off: only a frozen model that meets every benchmark goes to
  Sheehan for approval.

## 6. Limits

The measured values are blind figure readings (+-5-30). The model numbers are printed integers.
The (001) orientation and abrupt epitaxial interfaces are assumptions. The universal potential
under-predicts the bulk conductivity of Si (61.5 against about 145 W m^-1 K^-1 measured), Ge
(20.5 against about 60) and diamond, so every length is short by a factor of 2-3. The model is
scale invariant, so this affects G only through the length ratios between bands and between sides.
