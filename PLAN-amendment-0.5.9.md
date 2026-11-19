# PLAN amendment 0.5.9: stage 2a graded on orb-v3 (tbc-anhpairs-v1), and what Q1 means for the TiN diagnosis

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-19. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text and data: CC BY 4.0.

Grades job `tbc-anhpairs-v1` against amendment 0.2 (S2, S3, S5, written before stage 2a ran) and
amendment 0.5.5 (Q1, A0 recheck). Amendment 0.5.8 is reserved for the grades of `tbc-diag-v1`
(amendment 0.5.7 sec. 4), which is still queued. Raw outputs, unedited:
[results/anh/pairs/](results/anh/pairs/) (`run.txt`, `tba-summary.json`, `q1.txt`, `tests.txt`,
`versions.txt`, `sha256.txt`, `pairs-progress.txt`). Hand grades:
[results/anh/pairs/grades-0.5.9.json](results/anh/pairs/grades-0.5.9.json).

## 1. Grades

| | statement (amendment) | prior | result | grade |
|---|---|---|---|---|
| A0 | per-bin sum = phono3py conductivity (0.5.5 recheck) | - | ratio 1.000 | **PASS** |
| A1 | bulk kappa300 in experimental ranges (0.2) | 60 % | Si 62, C 888 low; GaN 141, AlN 169 in range | **FAIL** (as in 0.5.5) |
| Q1 | G unchanged (< 1 %) when one side's lengths are scaled; changes when relative lengths change (0.5.5) | 90 % | 3 of 3 checks pass | **PASS** |
| S2 | G_BTE(700 K) / G_BTE(400 K) < 1.10 for Al/Si and Al/GaN, oriented (0.2, main) | 85 % | Al/Si 1.016, **Al/GaN 1.143** | **FAIL** |
| S3 | oriented Al/Si, 300 K: conventional DMM < G_BTE < gray-Simons DMM (0.2) | 65 % | 456 < 648 < 710 | **PASS** |
| S5 | removing isotopes changes G_BTE at 300 K by < 5 % (0.2) | 70 % | Al/Si 0.997, **Al/GaN 1.218** | **FAIL** |
| S6 | orb-v3 and MACE-MPA-0 agree within 15 % (0.2) | 70 % | one potential so far | **UNDETERMINED** (tbc-anh-v2 queued) |

Priors met: Q1, S3 (and A0); missed: S2, S5 (and A1). We said S2 and S5 would pass at 85 % and
70 %; both failed on Al/GaN.

## 2. Reading S2 by amendment 0.2's rule

0.2: "FAIL means the non-scaling ingredients (quantum heat capacities, isotopes) already give a
rise", and we then measure how much of the data they explain. Measured (Li et al. 2022, MW m^-2 K^-1,
400 -> 700 K): Al/Si 400 -> 530 (x1.325), Al/GaN 432 -> 546 (x1.264). The elastic model with bulk
lifetimes gives x1.016 and x1.143. In log terms it explains **about 6 % of the Al/Si rise and about
57 % of the Al/GaN rise**. So S2 failing does not mean the elastic model accounts for the data: for
Al/Si it accounts for almost none of the rise. The data still need an inelastic interface channel
(stage 2b, `tbc_inel.py`) for Al/Si, as 0.2 said for the PASS branch. For Al/GaN about half of the
rise comes from the GaN side's high-frequency modes filling up with temperature (G_BTE climbs from
145 at 100 K to 692 at 700 K).

S5 fails for the same side: isotope scattering in GaN (strong, rising as omega^4) changes the
spectral shape of its lifetimes. By Q1, only that shape matters to G, so a 22 % effect is
consistent with the model's structure.

## 3. Information (not graded)

Oriented interfaces, 300 K, MW m^-2 K^-1:

| interface | DMM | G_BTE | + Al electron-phonon in series | measured |
|---|---|---|---|---|
| Al/Si | 456 | 648 | 422 | 379 (Li 2022) |
| Al/GaN | 578 | 522 | 365 | 423 (as printed by the job) |
| AlN/GaN | 482 | 304 | - | 324 (UCLA); 483 (Walwil) |

AlN/GaN at 100 K: G_BTE 103 against 117 measured (UCLA). With no fitted parameter, the two metal
interfaces land at +11 % and -14 % and AlN/GaN at -6 % (UCLA) at 300 K. Three interfaces at one
temperature are not a test and nothing here is graded.

## 4. What Q1 means for the TiN diagnosis (amendment 0.5.7)

Q1 passed: in this model, multiplying every phonon length on one side by a constant leaves G
unchanged. Two consequences, written before `tbc-diag-v1` returns:

1. **`s1x0.01` cannot test the electron-bath idea.** 0.5.7 added it (TiN lengths divided by 100) as a
   crude stand-in for TiN phonons pinned by electrons. By Q1 it must return the base value within
   about 1 %. If it does, that says nothing about electrons. If it does not, Q1 and the diagnostic
   disagree, which is a defect to chase. (It stays information only, as 0.5.7 said.)
2. **Stage 3c must model the electron bath as an energy exchange, not as shorter lengths.** In a
   metal, electron-phonon coupling g ties every phonon band to the electron temperature: each band
   exchanges heat with a third reservoir. In the two-stream equations that is a band-to-electron
   term, G_ep,j (T_j - T_e), with the electrons' own conductivity carrying heat away; shortening the
   phonon-phonon lengths is not equivalent. The design of stage 3c, pre-registered before any run,
   will add that term with g from the literature (TiN: about 1e18 W m^-3 K^-1, Khan et al.; Al: 2.45e17,
   already used in stage 1), split across bands by phonon density of states, and no fitted
   parameter.

Al/Si sits **above** its DMM (648 against 456) while every TiN interface sits 5-6 times **below**
(0.5.7). The two runs use the same solver, so D3 of `tbc-diag-v1` (uniform lengths, G_uni >= G_dmm)
remains the decisive check on whether the TiN deficit is a defect or the physics of band
non-equilibrium.

## 5. Next

- `tbc-diag-v1`: grade D0-D4 in amendment 0.5.8 by 0.5.7 sec. 4's rules.
- `tbc-anh-v2` (MACE-MPA-0): S6 and the MACE grades.
- Then stage 3c (electron-phonon exchange in metals) and stage 2b (inelastic channel), each
  pre-registered on TiN/X, Al/Si and Al/GaN before it runs.

## 6. Limits

The lifetimes come from one universal potential (orb-v3) that under-predicts silicon's and
diamond's bulk conductivity by about 2.5 times (A1). By Q1 that overall scale does not enter G, but
spectral-shape errors do, and S6 (the second potential) is still to come. The measured values are
figure readings from the cited papers.
