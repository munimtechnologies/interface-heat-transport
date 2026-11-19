# PLAN amendment 0.5: interfaces without a metal (ZnO/GaN, GaN/AlN), a data-only bound on inelastic transport, and stage 3a pre-registered

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-19. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

Written before any model output of stages 1, 2a, 2b or 3a exists (jobs `tbc-base-v1`, `tbc-anh-v1`
are queued; `tbc-zno-v1` is filed with this amendment). Nothing earlier is changed. New files:
`data/read-nonmetal.json`, `tbc_nonmetal.py` 0.1.0, `test_tbc_nonmetal.py`,
`results/nonmetal/bounds-exploratory.json`, `run_tbc_zno.py` 0.1.0.

## 1. Why interfaces without a metal

At Al/Si and Al/GaN, the measured rise of G with T could come from phonons crossing inelastically or
from the metal's electrons. Amendment 0.4 excluded the bulk-only electron route, but an
interface-specific one is still open. Two non-metals in contact have no electron route at all.
Whatever such an interface does above room temperature is phonon physics alone. So these interfaces
are a clean test of any inelastic phonon channel, including our stage 2b channel (`tbc_inel`).
AlN/GaN is also the one interface the challenge names as an example of its hidden set.

## 2. The readings, and what they allow (exploratory, not pre-registered)

MV2's hands read three papers on 2026-11-19. The reader saw no model result
(`data/read-nonmetal.json`, MW m^-2 K^-1, reading error about +-4 %):

| interface | 100 K | 200 K | 300 K | 350 K | 400 K | 450 K | 500 K | source |
|---|---|---|---|---|---|---|---|---|
| ZnO/GaN (epitaxial) | 201 | 366 | 491 | | 520 | | 502 | Gaskins et al., arXiv:1710.09525 |
| GaN/AlN (bonded wafer) | 117 | 264 | 324 | | | | | UCLA, APL 125, 032104 (2024) |
| GaN/AlN (epitaxial) | | | 483 | | 382 | | | Walwil et al., JAP 137, 095105 |
| AlN/SiC | | | 475 | 551 | | | | Walwil et al. |
| GaN/SiC | | | 307 | 373 | 408 | 447 | | Walwil et al. |

Model curves in Gaskins et al. for ZnO/GaN: AGF 300 (300 K) and 305 (500 K); DMM 240 and 278. The
abstract gives 490 [+150, -110] at room temperature. UCLA's AGF curve gives 346 at 300 K.

**The bound** (`tbc_nonmetal.py`). Write G(T) = E(T) + I(T). E is the elastic part. With fixed
harmonic force constants it cannot fall as T rises, because every mode's heat capacity grows. I is
any channel built from three-phonon rates: a sum over modes of C(nu, T) Gamma(nu, T) with
T-independent weights. For each decay term (1 + n' + n'') and each combination term (n' - n''),
the product with the mode's heat capacity grows by at least the classical factor rho = T2/T1. The
unit tests check this on a grid of 0.05-30 THz, all decay splits and partners, for 300-500, 300-400
and 400-500 K. The minimum is reached in the classical limit. So I(T2) >= rho I(T1), and

    I(T1) <= max(0, G(T2) - G(T1)) / (rho - 1).

The bound is tight when E is flat and I = cT (unit test).

**ZnO/GaN** (hand numbers; `results/nonmetal/bounds-exploratory.json` recomputed by the test):

- Central readings, 300 -> 500 K: G changes by +11, so **I(300 K) <= 16.5 MW m^-2 K^-1, 3.4 % of G**.
  From 400 to 500 K G falls (520 -> 502), so at central values such a channel has no room at all.
- Worst-case reading (G(300) down 4 %, G(500) up 4 %: 471.36 and 522.08): I(300) <= 76.08, 16.1 %
  of G. Using 400 -> 500 K the same way gives I(400) <= 91.52, so I(300) <= 68.64.
- A channel that grows only half as fast as classical three-phonon rates (rho = 4/3) gives 33.0
  (central) and 152.16 (worst case).
- The measurement sits 191 above AGF at 300 K. If that gap were a three-phonon channel, G(500 K)
  would be at least 305 + 191 x 5/3 = 623.3 (central). The reading is 502. The gap to DMM (251)
  would need at least 696.3.

So at central values the ZnO/GaN data leave room for a three-phonon channel of a few percent at
most. Gaskins et al. read the excess over AGF as a sign of "anharmonic phonon interactions
enhancing the TBC". Whatever lifts this interface above AGF, it does not grow with temperature
the way a three-phonon process must. It looks like a missing elastic contribution (for example
the zone-centre modes the authors point to) or something in the measurement.

**This is not decisive.** The published error bars (+31 %/-22 % at room temperature) are far larger
than the reading error. If the errors at 300 and 500 K were independent and each moved to its edge,
the bound would loosen to about the size of G. We use the reading error because TDTR errors at
one interface come mostly from shared inputs (transducer thickness, heat capacities, film
thickness). Those move the whole curve rather than its shape, but we have not checked that for
this data set. Our numbers come from arXiv v1 (2017). The published version (Nano Lett. 18, 7469
(2018)) has a correction that we have not read.

**Other interfaces.** GaN/AlN (Walwil) falls by 21 % from 300 to 400 K (ratio 0.791). No model with
a non-falling E and a growing I can follow that, and Walwil et al. themselves note that it
contradicts predictions. GaN/SiC (ratio 1.456 from 300 to 450 K, against T2/T1 = 1.5) and AlN/SiC
(1.16 from 300 to 350 K, against 1.167) rise almost in proportion to T. Our bound gives no useful
limit there: SiC's high Debye temperature lets E keep growing. UCLA's bonded GaN/AlN sits at 0.936
of their AGF curve at 300 K. So within the GaN family the measured behaviour above room temperature
is flat (ZnO/GaN), falling (GaN/AlN, one group) and close to T-linear (GaN/SiC). A single no-fit
model has to explain all three, or say which measurement it cannot explain and why.

## 3. Stage 3a, pre-registered: job `tbc-zno-v1`

`run_tbc_zno.py` imports stage 2a (`tbc_anh`, `run_tbc_anh`, `tbc_bte`) and stage 2b (`tbc_inel`,
`run_tbc_inel`) and registers wurtzite ZnO in their tables. Nothing is copied. For orb-v3 (primary)
and MACE-MPA-0, it computes three-phonon lifetimes of ZnO, GaN and AlN and the elastic two-stream
conductance G_BTE across [0001] (with isotopes). It adds the stage 2b channel G_inel at 100-700 K
with no fitted factor, so G_total = G_BTE + G_inel, and grades on the primary model:

- **Z0** (ZnO sanity): relaxed, fewer than 1 % of dense-mesh modes below -0.05 THz, highest
  frequency within 15 % of 13.1 THz (E2(high), 438 cm^-1; no LO-TO correction), kappa(300 K)
  in [20, 150] W m^-1 K^-1. Expected PASS (75 %).
- **Z1** (main; the data bound applied to the no-fit channel): G_inel(ZnO/GaN, 300 K) <= 76.08.
  Our pre-data estimate is about 20. ZnO's heat capacity is about 2.8 MJ m^-3 K^-1, the C-weighted
  rate about 5e11 s^-1, t = c/16 about 0.33 A per side, and Q about 0.4. That gives about 9 per side.
  Expected PASS (70 %).
- **Z2** (shape): G_total(500)/G_total(300) for ZnO/GaN <= 1.1076 (the worst-case reading ratio).
  Expected PASS (50 %).
- **Z3** (level): G_total(ZnO/GaN, 300 K) in [380, 640] (the published error bar). Expected FAIL
  low (60 %): AGF gives 300.
- **Z4** (level): G_total(AlN/GaN, 300 K) within 25 % of UCLA's 324, [243, 405]. Expected PASS (55 %).
- **Z5** (trend): G_total(AlN/GaN, 400)/G_total(300) >= 1. Expected PASS (95 %). This records
  that Walwil's fall lies outside this model class.

Also recorded, not graded: the MACE-MPA-0 grades, G at 100 and 200 K against 201 and 366, and the
3ph-only and classical-heat-capacity variants.

**Reading rules (fixed now).**
- Z1 FAIL: the no-fit channel predicts more inelastic transport than ZnO/GaN allows at 300 K, even
  at the worst-case reading. Stage 3 may not then use `tbc_inel` at full strength on non-metal
  pairs. A scale factor fitted to the Al pairs remains forbidden (PLAN.md section 5).
- Z1 PASS with Z3 FAIL low: the excess of the ZnO/GaN measurement over elastic models is neither
  captured by our elastic model nor allowed to be a three-phonon channel. It must then be missing
  elastic physics or come from the measurement. We would then ask Gaskins et al. (with a held
  draft, `--review`) whether their 300-500 K errors are correlated.
- Z2 FAIL: the model rises faster than ZnO/GaN does. Report which part (G_BTE or G_inel) carries
  the rise.
- Stage 3's frozen model takes its non-metal settings from these readings and nothing else.

## 4. Limits

- Readings come from figures (+-4 %), not data files, and from Gaskins et al.'s arXiv v1. The
  published 2018 version and its correction may differ.
- The bound assumes E does not fall with T. Anharmonic softening of real crystals (which our
  model does not include) could lower E by a few percent. A 3 % drop at 500 K (15 MW) would loosen
  the central bound from 16.5 to about 39.
- Three interfaces, three groups, and disagreeing trends within the GaN family. This amendment
  states a constraint and a test; it claims no mechanism.
