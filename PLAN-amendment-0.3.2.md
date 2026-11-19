# PLAN amendment 0.3.2: the electron side of Al interfaces, prior art for stage 2b, and Li et al.'s DMM

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-19. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

Written **before any stage 1, 2a or 2b model output exists**. Jobs `tbc-base-v1` and `tbc-anh-v1`
are still queued, and `tbc-inel-v1` has not been filed. Expectations B1-E4 (stage 1), A0-S6 (stage
2a) and I0-I4 (stage 2b) are unchanged. Sources are in
[readings/dnh2011-li2022-notes.md](readings/dnh2011-li2022-notes.md).

## 1. Prior art for stage 2b, and the "tension" with the scaling argument resolved

Duda, Norris and Hopkins (J. Heat Transfer 133, 074501, 2011; "DNH") explain the linear rise of
G with T in the classical limit by "three-phonon scattering events in which phonons propagating
across the interface interact with phonons on the other side of the interface". Their model,
h = h'_S + h'_T, contains a fixed interfacial length lambda and the bulk rates tau^-1 ~ A_j omega^n T,
so h = b + cT.

We had flagged DNH as possibly in tension with amendment 0.2's scaling argument (an elastic
interface between semi-infinite solids is unchanged when every bulk rate is multiplied by the same
factor). It is not. DNH's channel acts within a fixed length lambda of the interface, so it scales
with the rates and breaks that invariance, exactly as our stage 2b channel does. The scaling
argument covers bulk thermalization at an elastic interface. DNH's mechanism, and ours, is
inelastic scattering at the interface itself.

**Credit:** `tbc_inel.py` is a mode-resolved, first-principles version of the mechanism DNH proposed
in 2011. It uses bulk three-phonon rates, densities of states and the bond geometry of the
terminating layers, with nothing fitted to an interface. The README now says so.

## 2. G_ep(T) from bulk inputs: why the escape route of 0.3.1 is closed

Amendment 0.3.1 (section 4.2) found that Li et al.'s phonon-side conductance grows faster than T
and listed three escapes. One was "a G_ep that falls more slowly than theirs". Bulk inputs do not
support that escape:

- Al's lattice conductivity is three-phonon limited. With phonon-phonon scattering only it is about
  5.8 W m^-1 K^-1 at room temperature, and adding phonon-electron scattering leaves it "almost
  unchanged" (Wang, Lu and Ruan, J. Appl. Phys. 119, 225109, 2016). So Lambda_p ~ 1/T above the
  Debye temperature, and faster than that if four-phonon scattering matters.
- g of Al depends only weakly on temperature in DFT (Waldecker et al., PRX 6, 021003, 2016).
- With g = 2.5e17 W m^-3 K^-1 and Lambda_p = 6 W m^-1 K^-1, sqrt(g Lambda_p) = 1.22 GW m^-2 K^-1.
  This matches both our 1.21 GW and the 1.26 GW at 400 K implied by Li et al.'s G_pp. An
  independent set of Al inputs (De San Feliciano et al., arXiv:2510.27499: g = 246e15 W m^-3 K^-1,
  k_p = 8.95 W m^-1 K^-1) gives sqrt(g k_p) = 1.48 GW at 300 K, which is 1.29 GW at 400 K with
  Lambda_p ~ 1/T, within 2 % of Li et al.'s implied value.

So G_ep ~ T^-1/2 is the bulk-input expectation. Under it, the phonon side implied by the
measurements is (MW m^-2 K^-1; `tbc_gep.infer`, all values worked by hand and checked in
`test_tbc_gep.py`):

| T (K) | 300 | 400 | 500 | 600 | 700 |
|---|---|---|---|---|---|
| G_ep (T^-1/2, 1.26 GW at 400 K) | 1455 | 1260 | 1127 | 1029 | 952 |
| Al/Si G_pp | 512.5 | 586.0 | 738.1 | 892.8 | 1194.9 |
| Al/Si rise per 100 K | | 73.5 | 152.1 | 154.7 | 302.1 |
| Al/GaN G_pp | 594.4 | 657.4 | 743.6 | 896.3 | 1279.4 |
| Al/GaN rise per 100 K | | 63.0 | 86.2 | 152.7 | 383.1 |

**Finding (data plus bulk inputs, no model):** with G_ep ~ T^-1/2, the inferred phonon-side
conductance is **convex at every step from 300 to 700 K for both Al/Si and Al/GaN**. Its rise per
100 K grows, and from 400 to 700 K it is far from linear (fitted E + aT gives E = -226 MW for Al/Si
and -172 MW for Al/GaN, which is impossible for an elastic part). A saturated elastic part plus a
classical three-phonon channel (DNH's b + cT, or our stage 2b) gives a linear or concave G_pp. With
a constant G_ep (0.3.1) the curve is not convex at every step (Al/Si increments 46, 109, 84, 153).

**Why this is not decisive.** The phonon side amplifies measurement error by 1/(1 - G/G_ep): 1.5x
at 400 K and 2.3x at 700 K for both pairs. One 8 % error bar on G(700) becomes about 18 % on
G_pp(700). For Al/Si a G(700) 8 % lower (488 MW) gives G_pp(700) = 999 MW, and the 600-700 K step
becomes smaller than the 500-600 K one. The 400-500-600 K Al/Si steps are almost equal (second
difference 2.7 MW). The convexity rests mostly on the 700 K points and on the 300-400 K step.

**What can produce it, if real:**
1. a phonon channel that grows faster than T (four-phonon-like) at the interface;
2. **a direct electron-to-substrate-phonon channel in parallel with the series path.** De San
   Feliciano et al. (arXiv:2510.27499, accepted in PRB) find from first principles that
   phonon-phonon transport cannot explain the measured Al/Si conductance, and infer a direct
   interfacial electron-phonon channel of about one third of the total. If that channel rises with
   T, removing it before inverting the series term would remove much of the convexity;
3. error in the 700 K measurements.

## 3. Li et al.'s own DMM in series with G_ep

Li et al.'s Fig. 1 DMM curve was read as 460/475/490 (Al/Si) and 440/455/470 (Al/GaN) MW m^-2 K^-1
at 300/400/700 K. We take it to be phonon-only (the reader did not say). In series with G_ep ~ T^-1/2:

| pair | DMM + G_ep at 300 K | measured | DMM + G_ep at 700 K | measured |
|---|---|---|---|---|
| Al/Si | 349.5 | 379 (-8 %) | 323.5 | 530 (-39 %) |
| Al/GaN | 337.8 | 422 (-20 %) | 314.7 | 546 (-42 %) |

The elastic, conventional picture falls with T while the measurements rise. By 700 K it misses
about 40 % of G for both pairs. Our stage 1 expectation E2 (conventional DMM + G_ep < 379 at
Al/Si) is consistent with Li et al.'s own DMM. Li et al.'s DMM puts Al/GaN below Al/Si (ratio
0.96), while the measurement has Al/GaN above Al/Si at every temperature (1.11 at 300 K).

## 4. Consequence for stage 3 (design, not yet pre-registered)

A no-fit model of metal/non-metal interfaces must carry three things: the elastic phonon part,
an inelastic phonon channel (stage 2b), and the electron side. The electron side has two parts:
G_ep(T) = sqrt(g Lambda_p(T)) from bulk inputs (`tbc_gep.py`) and a possible direct interfacial
electron-phonon channel. The benchmark set separates them. Non-metal pairs (AlN/GaN, ZnO/GaN)
have no electron route. Au/Si has a negligible direct channel according to De San Feliciano et al.
Al/Si and Al/GaN have both routes. Stage 3 will be pre-registered with these contrasts before it
is run.

## 5. Recorded check: job `tbc-gep-v1`

`run_tbc_gep.py` 0.1.0 reproduces sections 2 and 3 from `data/meas-li2022.json` for three G_ep
laws (constant 1.21 GW; T^-1/2 and T^-1/4 from 1.26 GW at 400 K). It writes
`results/gep/gep-li2022.json`, and its last line starts `TBG:`. Expectation G0: the unit tests pass
and the T^-1/2 values equal the tables above to 0.1 MW m^-2 K^-1, with convex = True for both pairs.
This is a reproducibility check of hand arithmetic, not a test of physics.

## Files

`tbc_gep.py` 0.1.0, `test_tbc_gep.py`, `run_tbc_gep.py`, `readings/dnh2011-li2022-notes.md`.
