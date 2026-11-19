# PLAN amendment 0.4: grade of the blind test D1-D4 (a bulk-only electron route at Al interfaces)

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-19. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

This grades section 5 of `PLAN-amendment-0.3.3.md`, which was public before the readings were made.
Nothing in 0.3.3 is changed. New files: `data/read-d1d3.json` (the readings), `tbc_grade_d.py`
0.1.0 (grader), `test_tbc_grade_d.py`, `results/dtest/grade-d1d4.json` (the grade). The numbers
below were worked out by hand. The test recomputes them in code, and this commit was made only
because they agree (0.3 MW m^-2 K^-1 for G_min, 0.002 for slopes).

## 1. The readings (made blind to the thresholds and to every model result)

MV2's hands read both figures on 2026-11-19 (MW m^-2 K^-1):

| interface | 300 K | 400 K | 480 K | 500 K | 600 K | highest T | reading error |
|---|---|---|---|---|---|---|---|
| Al/Si Sample 2 (rough; Li et al. 2022, Fig. 1) | 308 | 329 | | 338 | 347 | 600 K | +-6 |
| Al/sapphire (epitaxial; Cheng et al. 2020) | 346 | 362 | 357 | | | 480 K | +-5 |

Cheng et al.'s own words: "the TBC reaches a plateau because the phonons in the Al side is fully
excited". They also wrote: "The estimated electron-phonon coupling thermal resistance within Al is
about 20% of the overall thermal resistance". Our bulk law G_ep(300 K) = 1455 puts 346/1455 = 24 %
of the resistance there. That is an independent check of the G_ep law used since amendment 0.3.2.

## 2. Arithmetic (G_ep(T) = 1260 (400/T)^0.5; G_ep(300, 400, 480, 600) = 1454.92, 1260, 1150.22, 1028.79)

G_min(T2) = s T2 + series(G_ep(T2), phonon_side(G(T1) - s T1, G_ep(T1))). Sergeev with DFT g:
s = 177.48/300 = 0.59159. DSF-size channel: s = 107/300 = 0.35667.

- **D1** (Sample 2, 400 -> 600 K, Sergeev). s T1 = 236.64; rest 92.36; phonon side 1/(1/92.36 - 1/1260)
  = 99.67; series with 1028.79 = 90.87; plus s T2 = 354.96 -> **G_min = 445.82**, threshold 378.95.
  Measured 347 (ratio 0.778). **EXCLUDED.** (Expected excluded, 75 %.)
- **D2** (same, DSF). rest 186.33; phonon side 218.67; series 180.34; plus 214.00 -> **G_min = 394.34**,
  threshold 335.19. Measured 347 (ratio 0.880). **NOT EXCLUDED.** (Expected excluded, 55 %: this
  expectation failed.)
- **D3** (Al/sapphire, 300 -> 480 K, span 180 K, so graded). Sergeev: rest 168.52; phonon side
  190.60; series with 1150.22 = 163.51; plus 283.96 -> **G_min = 447.47**, threshold 380.35. Measured 357
  (ratio 0.798). **EXCLUDED.** (Expected excluded, 70 %.) DSF: rest 239; phonon side 285.98;
  series 229.03; plus 171.20 -> G_min = 400.23, threshold 340.20. Measured 357. **NOT EXCLUDED.**
- **D4: TRIGGERED.** D1 and D3 both exclude Sergeev's model.

The other three normalisations of Sergeev's model (Lombard et al.'s Table 1, s = 0.88; from lambda,
0.978; their printed Eq. 4, 1.176) are excluded at both interfaces as well. At Sample 2 the channel
alone would already exceed the measured G(400) (352, 391 and 470 against 329). At Al/sapphire the
Table 1 and lambda versions need 503 and 522 at 480 K against 357. The printed Eq. 4 version would
already exceed the measured G(300) (353 against 346).

**Robustness.** G_min grows strictly with s (proved in one line: G_min - G(T1) = s (T2 - T1) minus
a series difference that shrinks as s grows; checked numerically in the test). We moved each
reading by its full reading error in the direction least favourable to each verdict. None of the
four verdicts changes: D1 with (323, 353) gives G_min 440.01 and threshold 374.0. D3 with (341, 362)
gives 442.76 and 376.3. D2 with (335, 341) gives 399.95 and threshold 339.96, so it is still not
excluded, by about 1 MW. D3-DSF with (351, 352) gives 404.82 and threshold 344.1. A G_ep that does
not fall with T (n = 0) would raise every G_min, because the series term would no longer drop. So
the T^-1/2 law is the conservative choice here.

## 3. What it means

**Pre-registered result.** Sergeev's inelastic electron-boundary channel uses only bulk properties
of the metal and grows in proportion to T. Even in its smallest normalisation (DFT g) it requires
these Al interfaces to gain at least 35 % (Sample 2, 400 -> 600 K) and 29 % (sapphire, 300 -> 480 K)
in conductance, whatever the substrate. Al on rough Si gains 5.5 % (400 -> 600 K) and epitaxial Al on sapphire gains 3 % (300 -> 480 K). The model is
excluded at both, with margins of 20-22 % below the predicted minimum. So it cannot be the
explanation of the rise Li et al. measured at sharp Al/Si and Al/GaN. A T-linear channel the
size of De San Feliciano et al.'s 107 MW at 300 K survives the pre-registered 0.85 margin. Our
expectation that it would fail was wrong.

**Exploratory (not pre-registered).** Take central readings and no margin. The largest
interface-independent T-linear channel the flat interfaces allow is s = 0.151 (Sample 2) and
0.148 (sapphire), that is **sigma_e(300 K) <= 44.5 MW m^-2 K^-1**. That is less than half of 107. With the
0.85 margin the bounds are s <= 0.419 and 0.454 (sigma_e(300 K) <= 126 and 136). Now remove the
largest central-value channel (s = 0.148) from Li et al.'s Sample 1. The phonon side left over is
434/467/555/626/771 MW at 300-700 K, a rise of 65 % from 400 to 700 K. Our pre-registered
expectation for elastic models is a rise under 10 % (stage 2a, S2; that job has not run yet). So the degeneracy amendment 0.3.3 found is broken for any
channel that does not depend on the interface. What makes sharp Al/Si and Al/GaN rise must depend on
the interface itself. That can be inelastic phonon transport (Li et al.'s reading), or an electron
route that is open only at sharp interfaces. Bulk-only theory gives neither.

**Consequence for stage 3 (as fixed in 0.3.3 section 4).** Stage 3 carries no Sergeev channel. The
electron side of metal pairs is G_ep(T) alone, unless a no-fit, interface-dependent estimate of a
direct channel is pre-registered and passes a test like this one.

## 4. Limits

- The values were read off published figures (+-5-6 MW), not taken from data files. The published
  measurement error bars are about +-8-10 %, which is what the 0.85 margin allows for.
- Cheng et al. show two samples (Sub100, Sub200). The reading gives one set of values for both. Cheng et al.
  describe the whole high-T behaviour as a plateau.
- G_ep(T) is a bulk law (constant g, three-phonon-limited lattice conductivity). Every film is far
  thicker than the electron-phonon coupling length of Al (about 6 nm).
- Sergeev's model is used as Lombard, Detcheverry and Merabia (arXiv:1501.03176) wrote it. A version
  that depends on the interface's roughness or bonding is not tested here.
- Two interfaces and one metal. Other metals (Au, Pt, Cu) on flat interfaces would test the
  model's metal dependence.
