# PLAN amendment 0.5.6: stage 3a graded on orb-v3 (ZnO/GaN, AlN/GaN), job tbc-zno-v1

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-19. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

Grades the expectations Z0-Z5 that amendment 0.5 wrote before job `tbc-zno-v1` ran. Nothing in
them is changed. Outputs as returned: `results/nonmetal/zno-v1/` (stdout, tests, versions,
sha256 and the four .npz files).

## 1. What ran

The job ran 72 minutes on 2 CPUs. All unit tests passed (`tests.txt`). The primary model, orb-v3,
finished all three crystals (ZnO, GaN, AlN) and both interfaces. The secondary model, MACE-MPA-0,
finished ZnO and was then killed by the 3 GB memory cap while relaxing GaN (the reported peak,
1.74 GB, is a sampled value; the kernel's kill is what counts). So **every graded expectation can
be graded on orb-v3**, and the MACE-MPA-0 interface numbers, which were recorded rather than
graded, do not exist.

Bulk conductivities at 300 K (W m^-1 K^-1, three-phonon + isotope, in-plane / along c):
orb-v3 ZnO 31.4 / 37.9, GaN 141.1 / 144.6, AlN 169.3 / 171.8; MACE-MPA-0 ZnO 37.3 / 41.7.

Interface conductance from the job (MW m^-2 K^-1; G_BTE = elastic two-stream, G_inel = the stage
2b channel with no fitted factor, our DMM in brackets):

| T (K) | ZnO/GaN G_BTE | G_inel | **G_total** | (DMM) | measured (Gaskins) | AlN/GaN G_BTE | G_inel | **G_total** | (DMM) | measured (UCLA, bonded) |
|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 86 | 0.7 | **87** | 211 | 201 | 103 | 0.1 | **103** | 211 | 117 |
| 200 | 155 | 5.9 | **161** | 290 | 366 | 222 | 3.2 | **225** | 394 | 264 |
| 300 | 190 | 12.9 | **203** | 314 | 491 | 304 | 8.9 | **312** | 482 | 324 |
| 400 | 205 | 19.9 | **225** | 324 | 520 | 351 | 15.2 | **367** | 527 | - |
| 500 | 213 | 26.7 | **240** | 329 | 502 | 379 | 21.4 | **401** | 553 | - |
| 600 | 217 | 33.4 | 251 | 332 | - | 396 | 27.5 | 424 | 568 | - |
| 700 | 220 | 39.9 | 260 | 334 | - | 407 | 33.4 | 441 | 577 | - |

## 2. Grades (orb-v3, as written in amendment 0.5)

| | expectation | result | grade |
|---|---|---|---|
| Z0 | ZnO sanity: relaxed; < 1 % of dense-mesh modes below -0.05 THz; highest frequency within 15 % of 13.1 THz; kappa(300) in [20, 150] | relaxed (conv=True, a = 3.286, c = 5.304 A); kappa 31.4; the mode checks are in `anh-orb-v3-ZnO.npz` and were not printed | **incomplete**: two of four parts pass, two not yet read |
| Z1 | G_inel(ZnO/GaN, 300 K) <= 76.08 | 12.9 | **PASS** (also below the central-value bound 16.5) |
| Z2 | G_total(500)/G_total(300) for ZnO/GaN <= 1.1076 | 240 / 203 = 1.182 | **FAIL** |
| Z3 | G_total(ZnO/GaN, 300 K) in [380, 640] | 203 | **FAIL** (low, as expected at 60 %) |
| Z4 | G_total(AlN/GaN, 300 K) in [243, 405] | 312 | **PASS** |
| Z5 | G_total(AlN/GaN, 400)/G_total(300) >= 1 | 367 / 312 = 1.176 | **PASS** (records that Walwil's fall lies outside this model class) |

Z0's two unread parts will be read from the published .npz by the next job that loads it; until
then Z0 is reported as incomplete, not as passed.

## 3. Reading, by the rules fixed in amendment 0.5

**Z1 PASS with Z3 FAIL low** (rule 2). The no-fit inelastic channel stays well inside what the
ZnO/GaN data allow (12.9 against a worst-case bound of 76.08 and a central bound of 16.5), so
stage 3 may use `tbc_inel` at full strength on non-metal pairs. But the measured ZnO/GaN
conductance is 2.4 times our total at 300 K and well above our DMM, and a three-phonon channel is
not allowed to fill that gap. By the rule written in advance, the gap must be missing elastic
physics or come from the measurement. The rule also says we ask Gaskins et al. whether their
300-500 K errors are correlated. That email is drafted and held; it will go out only with
Sheehan's approval, and not before 2026-12-02 (MV2 has one held author email outstanding).

**Z2 FAIL** (rule 3: say which part carries the rise). Both do: from 300 to 500 K, G_BTE rises
190 -> 213 (ratio 1.121, already above 1.1076 on its own) and G_inel 12.9 -> 26.7. The elastic
rise comes from heat capacity still growing in GaN above room temperature. The measurement is flat
(491 -> 502, 1.022) and falls from 400 to 500 K. As amendment 0.5 noted, no elastic part with fixed
harmonic force constants can fall with T.

**An observation, not graded.** Our ZnO/GaN total is a nearly constant fraction of the
measurement from 100 to 400 K (0.43, 0.44, 0.41, 0.43; 0.48 at 500 K). A missing channel that
grew with T would make that fraction fall. A constant fraction looks instead like elastic
transmission that is too small by a T-independent factor, or like a calibration offset in the
measurement. Our elastic G_BTE is also below our own DMM (190 against 314 at 300 K), which is
unusual for a well-matched pair, and worth checking in stage 3.

**AlN/GaN, the interface the challenge names.** With no fitted parameter, orb-v3 gives 103, 225
and 312 at 100, 200 and 300 K against UCLA's bonded-wafer 117, 264 and 324 (ratios 0.88, 0.85,
0.96). Our DMM gives 211, 394 and 482 (1.8, 1.5 and 1.5 times the measurement). On this one
interface the no-fit model is far closer than DMM. The 100 and 200 K comparisons were recorded,
not graded, and the epitaxial GaN/AlN of Walwil et al. (483 falling to 382 MW m^-2 K^-1 from 300
to 400 K) is not matched in level or trend.

**Bulk conductivities.** orb-v3 again comes out low against measured bulk values (roughly 55-65 %
of commonly quoted room-temperature values, about 50, 230-250 and 285-320 W m^-1 K^-1 for ZnO, GaN
and AlN; these reference values are from memory of the literature, not read for this amendment, as for Si in amendment 0.5.5, A1). The
interface conductance depends on these lifetimes only through G_inel and the two-stream
weighting, so this does not change the grades above, but it is the same systematic softness and
is carried into stage 3.

## 4. What is not re-run

The MACE-MPA-0 half of the job (recorded, not graded) is not re-filed now: the queued job
`tbc-anh-v2` already runs MACE-MPA-0 one material per process at 5 GB, and the orb-v3 grades
above are complete without it.

## 5. Limits

One model (orb-v3) graded; MACE-MPA-0 missing for GaN and AlN. Measurements are by-eye readings
of figures (+-4 %), ZnO/GaN from Gaskins et al.'s arXiv v1 (the 2018 published version has a
correction we have not read). The two GaN/AlN data sets disagree with each other in trend. Z0 is
incomplete. No thermal-transport specialist has checked any of this.
