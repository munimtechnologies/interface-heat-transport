# PLAN amendment 0.3.1: measured inputs recorded, a solver correction, and what the data alone say

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-18. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

Written after the measured conductances were read off Li et al. (2022) and **before any stage 1, 2a
or 2b model output exists** (jobs `tbc-base-v1` and `tbc-anh-v1` are still queued and `tbc-inel-v1`
has not been filed). The grades I0-I4 of amendment 0.3 are unchanged. This amendment does three
things. It records the measured inputs, fixes a solver error that would have mis-graded I1, and adds
recorded (ungraded) diagnostics.

## 1. Measured inputs (data/meas-li2022.json)

Li et al., Nat. Commun. 13, 4901 (2022), Fig. 1. The paper has no source-data file, so the values were
read off the figure by a reader who had seen no model output. Reading uncertainty is about +-5
MW m^-2 K^-1; the published error bars are about +-8-10 %. All three samples were measured up to about
705 K, so the upper temperature of amendment 0.3 stays at 700 K.

| T (K) | Al/Si sample 1 | Al/GaN | Al/Si sample 2 (less sharp) |
|---|---|---|---|
| 300 | 379 | 422 | 310 |
| 400 | 400 | 432 | 330 |
| 500 | 446 | 448 | |
| 600 | 478 | 479 | (highest: 600 K) |
| 700 | 530 | 546 | not measured |

(MW m^-2 K^-1.) Measured rise from 400 to 700 K: Al/Si 130, Al/GaN 114. The authors' phonon-phonon
conductance after removing their electron-phonon term (SI Fig. 10, read +-3 %): Al/Si 585 (400 K) and
1170 (700 K); Al/GaN 688 (400 K) and 1093 (600 K; 700 K not shown). The numbers for g and Lambda_p
are not stated. Their G_ep = sqrt(g Lambda_p) uses g from Waldecker et al. (PRX 2016) and Lambda_p from
Jain and McGaughey (PRB 2016).

## 2. Correction: the scale-factor solver (tbc_inel 0.1.1)

I1 solves for the scale factor s* that makes the model's rise G(700) - G(400) equal the measured one.
Here G = [1/G_ep + 1/(G_BTE + s G_inel)]^-1. With G_ep in series, the rise is **not monotonic in s**.
It grows, peaks, and falls back towards zero as both temperatures approach the G_ep limit.
`solve_scale` 0.1.0 bisected on [0, 1000] and returned infinity whenever the rise at s = 1000 was
below the measured one. With realistic inputs that is almost always the case, even when a small s
matches. I1 would then have read "too-weak" whatever the model computed. The check that found this
used made-up inputs (G_BTE = 400 flat, G_inel = 15 and 26 MW at 400 and 700 K), never model output.
For those inputs the rise is 36 MW at s = 1000, yet s between 10 and 20 gives the 60 MW target
(new unit test `test_solve_scale_series_peak`).

0.1.1 returns the smallest s that matches, found by a log scan from 1e-4 to 1000 followed by
bisection. It returns infinity only if no s in [0, 1000] matches, and `max_rise` reports the largest
rise the model can give. The measured Al/Si rise (130 MW) is close to that ceiling. With G_ep =
1.21 GW and a flat G_BTE of 400 MW, no inelastic strength reaches it: the model's rise peaks at about
127 MW (unit test `test_solve_scale_unreachable`). **Reading rule added:** if s* is infinite for a
pair *because the model's rise peaks below the measured rise* (peak at s < 1000), the reading is
recorded as "no-fit", not "too-weak". The rise and absolute level cannot then be matched by any
strength of this channel with this G_ep. I1 itself is graded exactly as written: it FAILS, because s*
is not in [0.5, 2].

## 3. Recorded diagnostics (not graded), written down before any model output

For each Al pair `run_tbc_inel` 0.1.1 also records:

- **s_pp**, the phonon-side linear scale. Remove the series term from the measurement,
  G_pp = 1/(1/G_meas - 1/G_ep), then s_pp = [dG_pp - dG_BTE] / dG_inel over 400 -> 700 K. It is linear
  and always finite, and it says how much stronger the channel must be to carry the phonon-side rise.
  `reading_pp` applies amendment 0.3's thresholds to it.
- **s_abs400** = [G_pp(400) - G_BTE(400)] / G_inel(400): the strength that matches the *level* at
  400 K. If s_pp and s_abs400 agree, one constant fixes both the level and the slope; if they differ
  a lot, the channel has the wrong shape, not just the wrong size.
- **s* with a falling G_ep**: G_ep(700) = G_ep(400) sqrt(400/700). This is Lambda_p proportional to 1/T,
  which is what Li et al.'s own G_ep implies (section 4).
- The model's peak rise and the s at which it peaks.

Our expectations, recorded now: s_pp > 8 for both pairs (70 %), because our pre-data estimate of
G_inel (10-15 MW at 300 K, so about 0.03-0.05 MW m^-2 K^-2) is far below the slopes in section 4.
The I1 reading is "no-fit" or "too-weak" (80 % together).

## 4. What the measurements alone imply (no model)

Take a phonon-side conductance G_pp(T) = E + aT: an elastic part E that has stopped changing above
400 K (Al's Debye temperature is 428 K) plus any channel proportional to T, which is what three-phonon
processes give in the classical limit. Two temperatures fix E and a. Then E = (7 G_pp(400) -
4 G_pp(700))/3.

| assumption for G_ep | pair | G_pp 400 -> 700 | a (MW m^-2 K^-2) | E | channel share at 400 K |
|---|---|---|---|---|---|
| none (G_ep infinite) | Al/Si | 400 -> 530 | 0.433 | 227 | 43 % |
| none | Al/GaN | 432 -> 546 | 0.380 | 280 | 35 % |
| 1.21 GW, constant (ours) | Al/Si | 597 -> 942 | 1.149 | 137 | 77 % |
| 1.21 GW, constant | Al/GaN | 671 -> 993 | 1.074 | 242 | 64 % |
| Li et al.'s own G_pp (SI Fig. 10) | Al/Si | 585 -> 1170 | 1.95 | -195 | impossible |
| Li et al.'s own G_pp | Al/GaN (400 -> 600) | 688 -> 1093 | 2.03 | -122 | impossible |

With constant G_ep, the phonon-side values in MW m^-2 K^-1 are as follows. Al/Si: 551, 597, 706,
789 and 942 at 300-700 K. Al/GaN: 647, 671, 711, 792 and 993.

What follows from this, with its limits:

1. **For any constant G_ep, a T-proportional channel carries at least 43 % (Al/Si) and 35 %
   (Al/GaN) of the phonon conductance at 400 K.** That share is 4(r - 1)/3 with r = G_pp(700)/G_pp(400),
   and r only grows as G_ep falls, so the minimum is at no series term. The elastic part is at most
   about 227 MW m^-2 K^-1 for Al/Si (no series term) and 288 for Al/GaN (reached at G_ep near 3 GW;
   E peaks where G_pp(700)/G_pp(400) = sqrt(7/4)). E >= 0 needs
   G_ep >= 0.94 GW (Al/Si) and 0.84 GW (Al/GaN). For comparison, MD and machine-learned-potential
   studies put the inelastic share at high temperature at 30-48 % for Si/Al (as summarized in
   arXiv:2401.09212) and about 40 % for GaN/AlN (arXiv:2503.05084).
2. **Li et al.'s own decomposition grows faster than T.** G_pp(700)/G_pp(400) = 2.00 for Al/Si
   (above 7/4) and G_pp(600)/G_pp(400) = 1.59 for Al/GaN (above 6/4). On top of a non-negative,
   T-independent elastic part, no channel proportional to T can do that. If their G_ep is right, the
   phonon side needs a superlinear term (four-phonon-like, ~T^2), an elastic part that grows with T,
   or a G_ep that falls more slowly than theirs. Their implied G_ep is 1.26 GW at 400 K and 0.97 GW
   at 700 K, consistent with Lambda_p proportional to 1/T. **This is not decisive.** A G(700 K) that
   is 8 % lower, within the published error bars, brings the Al/Si ratio to 1.68. We record it and do
   not claim it. We could not read the paper's full text this run, so we do not know whether the
   authors discuss it.
3. The Al/GaN curve is convex: with constant G_ep, G_pp rises by 39, 81 and 202 MW per 100 K from 400
   to 700 K. A classical three-phonon channel in series with a constant G_ep gives a concave G(T),
   not a convex one. This rests mainly on the single 700 K point.

What this means for the study: whatever explains the Li et al. data, it must already carry close
to half or more of the phonon-side conductance at 400 K. That rules out any "small correction to
DMM" picture before a single lifetime has been computed. It also makes the electron-phonon term as
important as the phonon model. Stage 3 should therefore treat G_ep(T) as part of the prediction
(from bulk g and Lambda_p) rather than as a fixed number.

## 5. Files

`tbc_inel.py` 0.1.1 (solver), `test_tbc_inel.py` (three new tests), `run_tbc_inel.py` 0.1.1
(no-fit reading, diagnostics of section 3), `data/meas-li2022.json`. The job `tbc-inel-v1` is filed
once stage 2a's lifetimes are published under results/anh/.
