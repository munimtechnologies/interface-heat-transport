# PLAN amendment 0.5.1: the ZnO/GaN bound of amendment 0.5 does not hold up, and three supplementary grades for stage 3a

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-19. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

Written before any output of `tbc-zno-v1` (or of stages 1, 2a, 2b) exists; all those jobs are still
queued. Nothing earlier is changed: the thresholds Z0-Z5 of amendment 0.5 stay as filed. New files:
`data/read-gaskins2018.json`, `tbc_dressed.py` 0.1.0, `test_tbc_dressed.py`,
`results/nonmetal/bounds-0.5.1.json`.

## 1. Summary

Amendment 0.5 said that at central readings the flat ZnO/GaN conductance (300-500 K) leaves room
for a three-phonon channel of at most 16.5 MW m^-2 K^-1 (3.4 %) at 300 K. That number is not
robust, for two reasons we found this run:

1. **The published error bars.** A blind reading of the published paper shows that the error bars
   are 95 % confidence contours of a separate fit at each temperature. If those errors are
   independent, the same data allow a channel of up to **229 MW m^-2 K^-1 (47 % of G)** at 300 K.
2. **The premise.** The bound assumed the elastic part cannot fall with temperature. That is true
   for the conventional Landauer conductance, but a TDTR fit reports the "dressed" conductance,
   referred to the local temperatures next to the interface. Shi, Yang, Fisher and Ruan
   (arXiv:1812.07910) say their dressed ZnO/GaN conductance falls at high temperature. We show
   below that this can happen exactly, with no inelastic physics at all.

So the ZnO/GaN measurement alone does **not** rule out a sizeable inelastic channel. Our earlier
reading of the gap above AGF as "missing elastic physics" was also anticipated by Shi et al.: their
non-equilibrium (dressed) treatment "nearly doubles" the DMM conductance at Al/Si and ZnO/GaN. Our
stage 2a/3a model `tbc_bte` already computes a dressed conductance (its one-band limit is the
Simons/Zeng-Chen corrected DMM), so the queued jobs are not affected. Only the data-only bound and
the reading rules built on it change.

## 2. What the published paper says (blind reading, `data/read-gaskins2018.json`)

MV2's hands read Gaskins et al., Nano Lett. 18, 7469 (2018) on 2026-11-19, without sight of any of
our files. One ZnO/GaN series (8 points, 79-500 K; the 180 and 95 nm films are averaged). Both axes
are logarithmic. Reading uncertainty is +-10 on points.

| T (K) | G | 95 % bar | arXiv v1 reading (amendment 0.5) |
|---|---|---|---|
| 300 | 490 | +107 / -111 (text: +150 / -110) | 491 |
| 385 (nearest point to 400) | 524 | +121 / -124 | 520 at 400 K (interpolated: 521) |
| 500 | 503 | +127 / -131 | 502 |

Model curves: AGF 300 / 305, DMM 250 / 280 (arXiv reading 240 / 278), phonon radiation limit
(zeta = 1, isotropic) 680 / 780, and the same with AGF dispersion 1190 / 1330, at 300 / 500 K. The
error bars come from "95 % confidence contour" fits over h(ZnO/GaN) and h(Al/ZnO), with the ZnO and
GaN conductivities fixed from control samples. Whether the errors are correlated between
temperatures is not stated. The published correction (Nano Lett. 19, 1408, 2019) only adds a
compiled Supporting Information PDF; no value changes. So the published data agree with our arXiv
reading to within 1 MW m^-2 K^-1.

**The bound under each error model** (`tbc_dressed.numbers()`, hand-checked):

| error model | largest I(300 K) | share of G |
|---|---|---|
| fully shared between temperatures (central published values) | 19.5 | 4.0 % |
| fully shared, plus the +-10 reading error | 49.5 | 10 % |
| independent, 300 and 500 K only (one-sided 95 %) | 232.1 | 47 % |
| independent, all three points, any non-decreasing elastic part (profile, Delta chi^2 = 2.71) | **229.05** | **47 %** |

Fitting at each temperature separately points to mostly independent errors. The shared inputs
(transducer thickness, heat capacities) that would move the whole curve together are not part of
the quoted bars. Until the authors say otherwise, the honest number is about 229.

## 3. The dressed conductance can fall (`tbc_dressed.dressed_dmm`)

Shi et al. define the measured temperature on each side as a weighted mean of modal
equivalent-equilibrium temperatures, T_lambda,1 = T_e,1 - t12 dT_e / 2 and
T_lambda,2 = T_e,2 + t21 dT_e / 2. For a DMM interface (t12 + t21 = 1 where both sides have modes,
0 elsewhere), with weights C M_i (heat capacity times number of modes), this reduces exactly to

    1/G_d = 1/G_0 - 1/(2 R_1) - 1/(2 R_2),

with G_0 the conventional DMM conductance and R_1, R_2 the radiation limits of the two sides. The
unit test checks the identity on 50 random spectra. It is the gray Simons formula, but with full
spectral sums, so it holds for any DMM spectrum. When the softer side (ZnO) has saturated and the
stiffer side's high modes (GaN up to about 22 THz) are still being excited, R_2 grows, the last term
shrinks and G_d falls. In the unit test's two-band example, G_0 and every heat capacity grow from
300 to 500 K, yet G_d falls from 0.8750 to 0.8287 (-5.3 %).

**Corrected bound.** If the elastic part changes by kappa_E = E(T2)/E(T1) and the channel grows by
at least rho = T2/T1,

    I(T1) <= max(0, G2 - kappa_E G1) / (rho - kappa_E).

kappa_E = 1 gives amendment 0.5. The worst a DMM interface can do without a model is crude: hold
G_0 fixed and grow both radiation limits by the 22 THz heat-capacity ratio (1.793 from 300 to 500
K), with G_d/G_0 = 2. That gives kappa_E >= 0.693. Even with fully shared errors, the bound then
rises to 167.7 (34 %). So the size of the channel the data allow depends on what the elastic part
does, and only a model can say that. Stage 3a computes it.

## 4. Supplementary grades of `tbc-zno-v1` (pre-registered now, before any output)

These are graded with `tbc_dressed.grade_supp` on `tbz-summary.json["results"][primary]`, the same
primary model as Z0-Z5. They add to Z0-Z5 and change none of them.

- **Z1g** (the corrected bound with the model's own elastic trend): with
  kappa_E = G_BTE(500)/G_BTE(300) for ZnO/GaN, the worst-case reading of amendment 0.5 (471.36 and
  522.08) gives B_g = max(0, 522.08 - 471.36 kappa_E) / (5/3 - kappa_E). PASS if
  G_inel(300) <= B_g. Expected PASS (60 %).
- **Z6** (Shi et al.'s prediction in our model): G_BTE(ZnO/GaN, 500 K) < G_BTE(400 K). Expected
  FAIL (65 %): isotope and three-phonon scattering in our two-stream model spread the non-equilibrium,
  and the published fall may need their specific weights.
- **Z7** (the non-equilibrium enhancement): G_BTE/G_DMM for ZnO/GaN at 300 K in [1.5, 2.5]. Expected
  PASS (65 %).

**Revised expectations for the filed grades (thresholds unchanged).** Amendment 0.5 expected Z3
(G_total(300) in [380, 640]) to FAIL low because AGF gives 300. That reasoning ignored that our G_BTE
is a dressed conductance, about twice DMM (2 x 250 = 500 for the published DMM). New expectation:
Z3 PASS (55 %). Z5 (AlN/GaN G(400)/G(300) >= 1) drops from 95 % to 85 %, since a dressed elastic
conductance can fall. Both priors in amendment 0.5 stay on record for calibration.

**Reading rules (replacing those of amendment 0.5 section 3 for Z1).**
- Z1 FAIL with Z1g PASS: the channel exceeds what the data allow only if the elastic part is flat;
  with the model's own elastic trend it fits. No conclusion against `tbc_inel`.
- Z1 FAIL and Z1g FAIL: the no-fit channel predicts more inelastic transport than the worst-case
  reading allows under the model's own elastic trend. This still holds only for errors shared between
  temperatures (section 2). Stage 3 then flags `tbc_inel` on non-metal pairs as disfavoured, not
  excluded.
- Z1 PASS: as before, no tension.
- The planned question to Gaskins et al. about correlated errors stays useful, but it is no longer
  needed to read Z1. It would go only as a held draft (`--review`).

## 5. Limits

- Readings are from a log-log figure (+-10 on points, +-20 on model curves).
- The profile limit assumes Gaussian, independent errors with 1 sigma = bar / 1.96 and a classical
  channel (I proportional to T). A channel that grows faster is constrained more tightly.
- The closed form is exact for DMM transmission with weights C M_i. Shi et al. write their weights as
  hbar omega D v f0; for a small temperature difference we take this to mean the heat-capacity weight.
  For AGF transmission (t12 + t21 is not 1) the dressed formula is different, and we do not bound it.
- Whether a TDTR fit reports exactly the dressed conductance depends on mean free paths compared with
  the film and the thermal penetration depth. We follow Simons, Zeng and Chen, and Shi et al. in
  assuming it does.
