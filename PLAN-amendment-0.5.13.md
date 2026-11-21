# PLAN amendment 0.5.13: tbc-cover-v1 graded: the gate failed on diamond's first band (so every other test is indicative only), the main test V2 failed by 2-5 %, and TiN on diamond and on Ge stays at about 0.6 of its own DMM once binning is removed as a cause

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-21. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

Grades of job `tbc-cover-v1` (tbc_cover.py 0.1.0; amendment 0.5.11 sec. 4, specification of 0.5.12
sec. 2), applied with the reading rules of 0.5.11 as written. Raw output: `results/tin/cover/`
(run.txt, cover-summary.json, tests, sha256, versions, job.json). Units: MW m^-2 K^-1.

## 1. Grades

| | test | prior | result | grade |
|---|---|---|---|---|
| V0 | gate: tests pass; first cover band = tbc-fill-v1's merged range in all 6 TiN rows; G_base, G_merge reproduced within max(0.5, 1 %) | 85 % | tests pass; G_base and G_merge reproduced in all rows; **first band 0-0.6 THz in both TiN/C rows against the merged 0-4.8 THz** (TiN/Si 0-1.0 and TiN/Ge 0-0.6 match) | **FAIL** |
| V1 | each side's model conductivity within 1 % of its binned value | 90 % | exact in every row (no conductivity lost, nothing filled) | PASS (indicative) |
| V2 | **main**: G_cover within 25 % of G_merge (TiN/C) and G_cover >= G_merge (TiN/Si, TiN/Ge) | 55 % | TiN/C 0.878, 0.911; TiN/Si 0.971, 0.977; TiN/Ge 0.946, 0.962 | **FAIL** (indicative) |
| V3 | G_cover / DMM_cover in [0.4, 1.4], all TiN rows | 70 % | TiN/C 0.601, 0.563; TiN/Si 0.974, 1.025; TiN/Ge 0.585, 0.624 | PASS (indicative) |
| V4 | bands >= 0.4 THz change G_cover by < 10 % | 60 % | TiN/C 1.011, 1.013; TiN/Ge 1.083, 1.082; **TiN/Si 1.121, 1.127** | FAIL (indicative) |
| V5 | ZnO/GaN G_cover / G_base >= 1.2 | 50 % | 255.4 / 189.9 = 1.345 | PASS (indicative) |

**Why V0 failed.** tbc-fill-v1's merge band in TiN/C runs to 4.8 THz because it closes only once
every filled TiN bin is covered. The cover rule closes a band as soon as each side holds one sampled
bin, which in TiN/C already happens at 0.6 THz; the next band is 0.6-4.8 THz. Amendment 0.5.11's
sentence "the first band is the merge band of tbc-fill-v1" was therefore a wrong expectation for
diamond, not a reconstruction error (G_base and G_merge reproduce). The rule says V0 FAIL makes V1-V5
indicative only, and we keep to it.

## 2. Numbers (300 K / high T)

| pair | base | merge | cover | DMM (cover bands) | measured (Khan 2024 / Li 2022, as the job compares them) |
|---|---|---|---|---|---|
| TiN/C 300 / 900 K | 43.3 / 52.8 | 192.9 / 235.8 | **169.3 / 214.8** | 281.7 / 381.3 | 449 / 910 |
| TiN/Si 300 / 700 K | 52.3 / 61.5 | 317.0 / 357.4 | **307.7 / 349.2** | 315.9 / 340.8 | 507 / 916 |
| TiN/Ge 300 / 700 K | 51.5 / 59.3 | 169.6 / 188.0 | **160.4 / 180.9** | 274.2 / 289.9 | 190 / 281 |
| Al/Si 300 K | 648.3 | 651.3 | 662.5 | 456.5 | 379 (Li 2022) |
| Al/GaN 300 K | 522.0 | 570.6 | 638.9 | 578.6 | 422 (Li 2022) |
| AlN/GaN 300 K | 303.6 | 416.6 | 443.2 | 484.5 | 324 |
| ZnO/GaN 300 K | 189.9 | 193.3 | 255.4 | 314.2 | 201 |

Cover plus the tbc-tin-v1 inelastic term (information): TiN/C 175 / 259, TiN/Si 364 / 522, TiN/Ge
222 / 369 against 449 / 910, 507 / 916, 190 / 281.

## 3. Reading (indicative, because V0 failed)

- **Binning is no longer the explanation.** Once each side's conductivity is conserved exactly
  (V1), the band rule moves TiN's G by only 2-12 % from the merge result (V2's ratios). The
  5-6x deficit of the base run was the binning artefact amendment 0.5.11 found; what is left is
  not.
- **What is left:** on diamond and Ge, TiN's elastic-phonon G sits at 0.56-0.62 of its own DMM; on
  Si it equals DMM. This is inside the range other Boltzmann studies report (V3), so it is not
  yet a contradiction.
- **Against measurement** (cover alone; cover plus the inelastic term in brackets), model /
  measured is: TiN/C 0.38 and 0.24 (0.39, 0.28); TiN/Si 0.61 and 0.38 (0.72, 0.57); TiN/Ge 0.84
  and 0.64 (1.17, 1.31). For the earlier pairs at 300 K the cover value is too HIGH: Al/Si 1.75,
  Al/GaN 1.51, AlN/GaN 1.37, ZnO/GaN 1.27. It does not beat DMM on these benchmarks, so **it is
  not a candidate Grand Challenge entry**, and we say so.
- **V2 FAIL rule (0.5.11):** a denser lifetime mesh (TiN at 24^3 or more) is needed before any band
  rule is trusted. **V4 FAIL** on TiN/Si says the same: G there still depends on band width by 12 %.

## 4. What happens to this lane

The measured TiN/C and TiN/Si values are far above every elastic-phonon number we can make, and
above DMM itself. That points at channels the model does not have (electron-phonon coupling across
the metal interface, or inelastic phonon processes larger than our estimate), not at
bookkeeping. Lane B is **parked**: no new job until there is a concrete physical model for the
missing channel that can be pre-registered with no fitted parameter. The denser-mesh job is noted,
not filed; it would sharpen a number that is already known to be too low.

## 5. Limits

As in 0.5.11 sec. 5: the universal potential under-predicts bulk conductivity of Si, Ge and diamond
2-3x; the measured values are blind figure readings (+-5-30); lifetime meshes are coarse. All
work was carried out with AI agents and has not been checked by a specialist.
