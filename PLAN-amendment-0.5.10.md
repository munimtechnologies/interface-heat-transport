# PLAN amendment 0.5.10: stage 1 graded (tbc-base-v2): the elastic models saturate and bracket Al/Si, but do not get Al/GaN's ordering wrong; the two potentials disagree by 16 %

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-20. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

Grades job `tbc-base-v2` against PLAN.md (B1-B3, E1-E4) and amendment 0.5.4 sec. 4 (M0, N1), all
written before it ran. Thresholds unchanged. Raw output, unedited: `results/base-v2/`
(`tbc-summary.json`, `run.txt`, `tests.txt`, `versions.txt`, `sha256.txt`, `progress.txt`,
`job.json`). `tbc-fill-v1` and `tbc-anh-v2` are still queued; nothing here changes their rules, and
no module was changed.

## 1. The run

python:3.12, `run_tbc_safe.py` 0.1.0 with `run_tbc.py`, `tbc_transport.py`, `tbc_phonons.py` 0.1.0
fetched from main, MACE-MPA-0 (medium) and orb-v3 (conservative, inf, omat). Run step exit 0 after
682 s; peak child memory 3.18 GB (MACE) and 2.41 GB (orb-v3). One MACE child (diamond, float64)
was killed at the memory cap and, as the driver was built to do, retried once in float32, which
finished. The job runner records the job as "ran out of memory" because of that kill; every
material finished. 17 of 17 unit tests passed.

## 2. Grades

| | statement | prior | result | grade |
|---|---|---|---|---|
| M0 | tests pass and >= 9 of 10 materials per potential | 80 % | 17/17 tests; 10/10 and 10/10 (MACE diamond in float32) | **PASS** |
| B1 | all unit tests pass | - | 17/17 | **PASS** |
| B2 MACE-MPA-0 | every material within 15 % of its reference top frequency | FAIL expected (0.5.4) | Al 0.82, Si 0.81, GaN 1.26, AlN 1.32 | **FAIL** (as expected) |
| B2 orb-v3 | same | none | GaN 1.28, AlN 1.27; all others within 15 % (Al 9.88 against 9.7 THz, Si 15.37 against 15.6) | **FAIL** |
| B3 | MACE and orb-v3 isotropic full-band DMM for Al/Si at 300 K within 15 % | - | 369 against 441 MW m^-2 K^-1 (16.3 %) | **FAIL** |
| E1 | oriented full-band DMM G(700 K)/G(400 K) < 1.10, Al/Si and Al/GaN | 95 % | MACE 1.022 / 1.022; orb 1.027 / 1.029 | **PASS** (both) |
| E2 | Al/Si 300 K with Al's electron-phonon term in series: conventional DMM < 379 < Simons-corrected DMM | 55 % | MACE 281 < 379 < 391; orb 320 < 379 < 426 | **PASS** (both) |
| E3 | oriented DMM G(Al/GaN)/G(Al/Si) < 1.12 at 300 K | 70 % | MACE 1.135 (Simons 1.128); orb 1.275 (1.370); measured 1.116 | **FAIL** (both) |
| E4 | Cu/TiC/diamond over Cu/diamond > 1.48, conventional and corrected | 65 % | MACE 6.20 / 7.31; orb 7.64 / 8.91 | **PASS** (both) |
| N1 | MACE's oriented full-band DMM for Al/Si at 300 K below orb-v3's | 65 % | 366 < 434 | **PASS** |

Expectations met: M0, the certain B2 failure for MACE, E1, E2, E4 and N1 (6 of 7 written with a
prior above 50 %); E3 was missed. B2's failures for GaN and AlN are against PLAN.md's references,
which amendment 0.5.4 sec. 3(b) showed are the wrong mode (E1(TO) instead of the silent B1(high)
mode, 20.67 THz measured for GaN). Against B1(high), GaN reads 1.021 (MACE) and 1.039 (orb-v3)
(information, printed by the job). The grade stands as written.

## 3. Meaning (PLAN.md and amendment 0.5.4, written before the job), applied

- **E1 PASS**: elastic models cannot give the measured rise of G with temperature (Li 2022, Al/Si
  379 -> 530 MW m^-2 K^-1 from 300 to 700 K; both potentials give +2-3 % from 400 to 700 K). An
  inelastic channel is needed. This is what stages 2-3 (`tbc-anh`) already pursue.
- **E2 PASS**: at 300 K the non-equilibrium (Simons) correction alone brackets the Al/Si
  measurement for both potentials. The frozen model must choose between the two readings from
  physics, not by picking the closer one.
- **E3 FAIL**: the plan expected the elastic model to get the Al/GaN vs Al/Si ordering wrong
  because of GaN's acoustic-optical gap. It does not: MACE gives 1.135 against the measured 1.116,
  and orb-v3 overshoots (1.275). So the plan's "clear target for the inelastic channel (it must
  lift Al/GaN more than Al/Si)" is **withdrawn**; with orb-v3 the inelastic channel must instead
  lift Al/Si more.
- **E4 PASS**: the elastic DMM overstates the TiC interlayer gain 4-6 times (6.2-8.9 against the
  measured 1.48). The interlayer benchmark needs physics inside TiC (its own electron-phonon
  coupling) or at the TiC/diamond bond, as the plan said.
- **N1 PASS with B3 FAIL**: amendment 0.5.4's rule applies. From now on orb-v3, which reproduces
  every other measured top frequency within 6.1 % (Al +1.9 %, Si -1.5 %, Cu +2.1 %, C +2.4 %,
  Ni -0.2 %, Au -6.1 %) and GaN's B1(high) within 3.9 %, is the potential for
  the phonon side of any new stage; MACE-MPA-0 is reported beside it. Jobs already queued
  (`tbc-anh-v2`, `tbc-fill-v1`) keep the choices they were pre-registered with.

## 4. Observations (after the output; not tests)

- MACE-MPA-0's 18 % softening of Al and Si lowers its DMM values for Al/Si by about 16 % against
  orb-v3 (369 against 441 isotropic), which is roughly what a sound-speed scaling predicts.
- At 300 K orb-v3's conventional oriented DMM is above the measured values for Al/Si (434 against
  379) and Al/GaN (554 against 423), before any electron-phonon resistance; with Al's
  electron-phonon term in series it is 320 for Al/Si (below 379). The electron side therefore
  matters for any no-fit prediction with orb-v3.

## 5. Next

`tbc-fill-v1` (F0-F4) and `tbc-anh-v2` are graded when they return (amendment 0.5.11 onwards). Any
new stage is pre-registered with orb-v3 as the primary phonon potential.
