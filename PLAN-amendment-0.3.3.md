# PLAN amendment 0.3.3: a direct electron route into the substrate, what Li et al.'s data allow, and a blind test

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-19. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

Written **before any stage 1, 2a or 2b model output exists** (jobs `tbc-base-v1`, `tbc-anh-v1`,
`tbc-gep-v1` are queued) and **before the readings that section 5 grades**. Earlier expectations are
unchanged. New code: `tbc_direct.py` 0.1.0, `test_tbc_direct.py` (every expected number worked by
hand), `run_tbc_direct.py` (writes `results/direct/direct-li2022.json`; last line `TBD:`).

## 1. The circuit

Amendment 0.3.2 ended on the main escape left open: a direct channel from the metal's electrons into
the substrate's phonons, running **in parallel** with the usual path (electrons -> metal phonons ->
interface). We solved the steady-state two-temperature model of a semi-infinite metal with both
interfacial channels (sigma_e direct, sigma_p phonon-phonon). The effective conductance, referred to
the extrapolated far-field metal temperature as a TDTR fit on a thick film sees it, is given in closed
form in `tbc_direct.circuit_exact`. To first order in k_p/k_e it is

    G = sigma_e + [1/sqrt(g k_p) + 1/sigma_p]^-1

This is the mixed series-parallel circuit of Lu, Wang and Ruan (Phys. Rev. B 93, 064302, 2016). We
re-derived it, and the credit is theirs. With Al inputs (g = 2.46e17 W m^-3 K^-1, k_p = 8.95, k_e = 228
W m^-1 K^-1; arXiv:2510.27499), sigma_e = 107 and sigma_p = 210 MW m^-2 K^-1, the exact form gives 293.2
and the limit gives 291.0 MW m^-2 K^-1. The 0.8 % difference is negligible here. The steady-state formula in De San
Feliciano et al.'s appendix (their Eq. 53, G ~ G_pp + G_ep minus a small term) leaves out the
sqrt(g k_p) bottleneck on the phonon path. Their 107 MW value comes from a transient fit, not from
that formula, so it is not affected.

## 2. The only no-fit theory of the direct channel: Sergeev's model

Sergeev's inelastic electron-boundary scattering (Phys. Rev. B 58, R10199, 1998), as written by
Lombard, Detcheverry and Merabia (arXiv:1501.03176, Eqs. 1-3), uses only bulk properties of the metal.
It grows in proportion to T. Their printed Eq. 4 drops zeta(3), and their Table 1 does not follow
from it for every metal. So we evaluate Al at 300 K in four ways (MW m^-2 K^-1):

| version | sigma_e(300 K) | slope s (MW m^-2 K^-2) |
|---|---|---|
| tau_ep from DFT g = 2.46e17 (most conservative) | 177.5 | 0.592 |
| tau_ep from lambda<omega^2> = 185.9 meV^2 (implies g = 4.07e17) | 293.5 | 0.978 |
| their printed Eq. 4 (no zeta(3)) | 352.8 | 1.176 |
| their Table 1 | 264 | 0.880 |

## 3. What the measurements allow (Li et al. 2022, Fig. 1, read blind; no model)

With G_ep(T) = 1.26 GW (400/T)^1/2 (amendment 0.3.2) and a direct channel sigma_e = s T:

**Bound A (no assumption on the phonon side):** sigma_e(T) < G(T) at every T, so
s < 0.757 (Al/Si, binding at 700 K) and s < 0.780 (Al/GaN). That means sigma_e(300 K) < 227 and 234.

**Bound B (the phonon side never drops below the ab initio harmonic value of 210 MW at Al/Si,
arXiv:2510.27499).** Then s T <= G(T) - series(G_ep(T), 210). The series term is
183.5/180.0/177.0/174.4/172.1 at 300-700 K. So **s <= 0.506** (binding at 600 K):
**sigma_e(300 K) <= 151.8 MW m^-2 K^-1.**

**Sergeev's model with that floor (Al/Si), predicted G against measured:**

| T (K) | 300 | 400 | 500 | 600 | 700 |
|---|---|---|---|---|---|
| measured | 379 | 400 | 446 | 478 | 530 |
| Sergeev, DFT g (s = 0.592) | 361.0 | 416.6 | 472.8 | 529.4 | 586.2 |
| Sergeev, Table 1 (s = 0.880) | 447.5 | 532.0 | 617.0 | 702.4 | 788.1 |
| DSF-size T-linear channel (107 at 300 K) | 290.5 | 322.7 | 355.3 | 388.4 | 421.7 |

The Table 1, lambda and Eq. 4 versions overshoot by 18-49 % (excluded well beyond the 8-10 % error
bars). The DFT-g version overshoots by 10-11 % at 600-700 K (about 1.3 sigma), which is tension, not
exclusion. The DSF-size channel fits under the measurements.

**Degeneracy.** Remove a T-linear channel from G and invert the series term. The phonon side that
is left (MW m^-2 K^-1):

| s (MW m^-2 K^-2) | pair | 300 | 400 | 500 | 600 | 700 |
|---|---|---|---|---|---|---|
| 0 (amendment 0.3.2) | Al/Si | 512.5 | 586.0 | 738.1 | 892.8 | 1194.9 |
| 0.25 | Al/Si | 384.3 | 393.8 | 448.8 | 481.5 | 565.9 |
| 0.30 | Al/Si | 360.6 | 360.0 | 401.4 | 419.5 | 481.9 |
| 0.357 (DSF size) | Al/Si | 334.5 | 323.4 | 351.0 | 355.1 | 397.3 |
| 0.357 (DSF size) | Al/GaN | 402.0 | 375.6 | 354.5 | 356.9 | 430.2 |
| 0.50 | Al/Si | 271.8 | 237.7 | 237.3 | 215.2 | 221.9 |

**Finding (data plus bulk inputs, no fitted model).** Take a T-linear direct channel of the size De San
Feliciano et al. inferred at 300 K (107 MW). The phonon side left over is then **flat to within
+-5 % from 300 to 600 K** for Al/Si and within -12 % for Al/GaN. Only the 700 K points rise, by 12 %
and 20 %, which is inside the 2.3x-amplified error at 700 K (about 18 %). So the rise Li et al.
measured does not require inelastic phonon transport. An elastic phonon side plus a direct channel of
already-published size reproduces it. Li et al.'s interpretation (inelastic phonons) and a direct
electron route are **degenerate on these data**. Channels with s >= 0.5 would need a phonon side that
falls with T, which elastic and inelastic phonon models both rule out. That gives a third, weaker bound.

The idea that a cross-interface electron route can mimic a high-T rise is not new (silicide/Si work,
arXiv:1609.01776, says so in words). What is new here is the quantitative test against Li et al.'s
Al/Si and Al/GaN data with bulk G_ep(T) and Sergeev's model. We found no paper that does it.

## 4. Consequence for the challenge model (stage 3 design)

Metal/semiconductor data alone cannot separate the two mechanisms. Stage 3 will therefore be graded in two
steps. (i) The **phonon side** (elastic + stage 2b inelastic) is tested first on non-metal pairs
(GaN/AlN 320 and 490; ZnO/GaN ~490 at 300 K), where no electron route exists. (ii) It is then
frozen, and the metal pairs test the electron side. A direct channel enters stage 3 only through a
no-fit estimate that passes section 5. If Sergeev's model fails there, no direct channel is
included and the electron side is G_ep(T) alone.

## 5. Pre-registered blind test: does the same Al show the rise where the interface is different?

Sergeev's channel depends only on the metal. Where the same Al sits on a rougher Si interface (Li et
al.'s Sample 2), or on sapphire (Cheng et al., Commun. Phys. 2020, doi:10.1038/s42005-020-0383-6, arXiv:1906.05484), it
must make G rise with T as well. Given G(T1), the least G(T2) the circuit allows, with a phonon side
that does not fall between T1 and T2, is

    G_min(T2) = s T2 + series(G_ep(T2), phonon_side(G(T1) - s T1, G_ep(T1)))   (tbc_direct.min_rise_prediction)

G_ep(T) is the same Al law as above. The values are read by MV2's hands from the papers' figures,
after this file is public. Each test is graded with the read G(T1); "excluded" means measured
G(T2) < 0.85 G_min(T2).

- **D1** Li Sample 2, T1 = 400 K, T2 = highest measured (about 600 K). Sergeev with DFT g (s = 0.5916).
  With G(400) = 330 the threshold is G_min(600) = 446.8, so excluded if G(600) < 379.8.
  **Expect: excluded (75 %).**
- **D2** Same, DSF-size T-linear channel (s = 0.3567): G_min(600) = 395.3, excluded if < 336.0.
  **Expect: excluded (55 %).**
- **D3** Cheng et al. Al/sapphire. T1 = 300 K (or the nearest measured point at or above it),
  T2 = highest measured. Graded only if T2 - T1 >= 100 K; otherwise UNDETERMINED. Both slopes.
  **Expect: Sergeev DFT-g excluded (70 %).**
- **D4** If D1 and D3 both exclude Sergeev's model: no bulk-only direct channel of Sergeev's size is
  present at Al interfaces. Any direct channel that explains Sample 1's rise must then depend on the
  interface (sharpness, bonding), just as Li et al.'s inelastic reading does. Stage 3 then carries no
  direct channel.

**Disclosure.** Before writing this we had seen two qualitative statements: Li et al. say Sample 2
"saturates" above the Debye temperature, and a search summary said Cheng et al.'s Al/sapphire
conductance "levels off". We had seen no numbers from either above 400 K. The priors above
include this knowledge. The thresholds are what this file fixes in advance.

## Files

`tbc_direct.py` 0.1.0, `test_tbc_direct.py`, `run_tbc_direct.py`. Job `tbc-direct-v1` reproduces
sections 2-3 and 5 (expectation H0: tests pass, and the TBD line matches the tables above to 0.2
MW m^-2 K^-1).
