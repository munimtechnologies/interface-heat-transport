# PLAN amendment 0.5.2: interfaces up to 1000 K (TiN on diamond, Si and Ge) and stage 3b

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-19. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

Written before any output of `tbc-tin-v1` exists (and while stages 1, 2a, 2b and 3a are still
queued). Nothing earlier changes. New files: `data/read-hight.json`, `tbc_tin.py` 0.1.0,
`test_tbc_tin.py`, `run_tbc_tin.py` 0.1.0, `results/tin/bounds-0.5.2.json`.

## 1. Why these data

The Grand Challenge asks for a no-fit model of G(T). Most interface data stop at 300-500 K, where
elastic and inelastic pictures are hard to tell apart. Khan, Shi, Feser and Wilson (Appl. Phys.
Lett. 125, 041601, 2024; free copy OSTI 2578597) measured epitaxial TiN on diamond, SiC, Si and Ge
by TDTR up to about 1000 K: "the interface conductance increases monotonically without stagnating at
higher temperatures". They read this, and TiN/diamond lying about 25 % above the elastic limit at
300 K, as evidence that inelastic processes matter in all four systems. These are the hardest
benchmarks we know of, and TiN/diamond reaches the highest temperature.

MV2's hands read the figures on 2026-11-19 without sight of any of our files
(`data/read-hight.json`; +-5-15 on clear markers, +-20-30 on the partly hidden TiN/Si points above
480 K, +-2-3 % on curves):

| T (K) | 300 | 400 | 500 | 600 | 700 | 800 | 900 | last point |
|---|---|---|---|---|---|---|---|---|
| TiN/diamond | 449 | 572 | 662 | 740 | 814 | 870 | 910 | 948 at 996 K |
| TiN/SiC | 748 | 842 | 900 | 920 | 968 | | | 999 at 771 K |
| TiN/Si | 507 | 633 | 736 | 870 | 916 | | | 942 at 771 K |
| TiN/Ge | 190 | 228 | 244 | 263 | 281 | | | 299 at 771 K |

Curves drawn in the paper (300 K / 1000 K): diamond panel, elastic limit 364 / 505, DMM 259 / 352,
TiN maximum conductance 1290 / 1510, diamond maximum 3200 / 8200; Si panel, elastic limit
411 / 438, DMM 291 / 311, Si maximum 1065 / 1230, TiN maximum 1280 / 1500.

The same reading took Koh et al. (PRB 102, 205304, 2020): Al/sapphire 315 / 338 / 325 at 300 / 400
/ 498 K (flat, as Cheng et al. 2020 found: 346 / 362 / 357), Co/sapphire 370 / 405 / 460, Ru/sapphire
351 / 340 / 361. Their modal non-equilibrium Landauer model with DMM transmission (Shi's method) gives
Al 334 / 351, Co 372 / 390, Ru 296 / 307 at 300 / 500 K. Within 6 % for Al and Co, it is a strong
elastic benchmark that stage 2a's dressed `tbc_bte` should approach.

## 2. Exploratory numbers from the data alone (`tbc_tin.numbers()`, hand-checked)

Not pre-registered; they explain why stage 3b asks what it asks.

**(a) The room-temperature excess is not evidence for inelastic transport.** A TDTR fit reports a
"dressed" conductance referred to the temperatures next to the interface (amendment 0.5.1). Its
gray form, 1/G_d = 1/G_0 - 1/2R_1 - 1/2R_2, holds for any transmission that obeys detailed
balance. Put in the paper's own curves:

| | 300 K | 1000 K | growth |
|---|---|---|---|
| TiN/diamond elastic limit, as drawn | 364 | 505 | 1.39 |
| **same, dressed** | **453.8** | **629.7** | **1.387** |
| TiN/diamond DMM, dressed | 301.5 | 408.4 | 1.355 |
| TiN/Si elastic limit, dressed | 635.7 | 648.0 | 1.019 |
| TiN/Si DMM, dressed | 388.2 | 403.9 | 1.041 |
| measured TiN/diamond | 449 | 948 (996 K) | 2.11 |
| measured TiN/Si | 507 | 942 (771 K) | 1.86 |

At 300 K the dressed elastic limit (454) equals the measured TiN/diamond value (449). The 25 % excess
over the drawn elastic limit is what non-equilibrium at the interface gives, with no inelastic
scattering. The gray formula is an estimate, not a bound (with spectral transmission the dressing can
be larger or smaller), so this does not show that the 300 K value is elastic. It shows that the 300 K
excess cannot decide the question.

**(b) The rise can decide it.** Between 300 and 996 K the dressed elastic limit grows 1.39-fold;
the measured TiN/diamond conductance grows 2.11-fold. If the elastic part grows no more than that, the
non-elastic part at 996 K is at least 948 - 1.387 x 449 = **325 MW m^-2 K^-1 (34 %)**; for a channel
proportional to T, at least 168 (37 %) already at 300 K. For TiN/Si (dressed elastic growth 1.02)
the same floor is 425 (45 %) at 771 K.

**(c) A model-free floor for the conventional elastic conductance.** With temperature-independent
transmission and no thermal expansion, an elastic conductance is a sum of mode heat capacities with
fixed weights. It can grow at most as fast as the heat capacity of the highest mode that can cross
(the softer side's top frequency; TiN 18 THz, Si 15.6 THz, Ge 9.1 THz). That caps the growth at
1.797 (TiN/diamond, 300-996 K), 1.517 (TiN/Si, 300-771 K) and 1.158 (TiN/Ge, 300-771 K), leaving
non-elastic floors of **141, 173 and 79 MW m^-2 K^-1 (15 %, 18 %, 26 %)**. Moving every reading by
its error against the claim leaves 99, 98 and 47. This holds for the conventional conductance only.
Dressing can grow faster when the radiation limits grow more slowly than G_0, so it is not a bound on
what TDTR reports. Thermal softening is also left out: it lowers frequencies at high temperature, raises the
heat capacities there, and so loosens these floors.

So Khan et al.'s conclusion survives at high temperature, under any of these elastic pictures. What
nobody has shown is whether a **no-fit** inelastic model gets the size of the rise right. That is
stage 3b.

## 3. Stage 3b: job `tbc-tin-v1` and its pre-registered grades

`run_tbc_tin.py` reuses stages 2a and 2b unchanged (`tbc_tin.register_tin` adds rocksalt TiN,
a = 4.2417 A, and diamond-structure Ge, a = 5.6579 A, to their tables, extends the temperatures to
100-1000 K and adds (001) terminating layers: rocksalt f = 1/6, sigma/n = a/2; diamond f = 2/4,
sigma/n = a/4). For orb-v3 (primary) and then MACE-MPA-0, it computes the three-phonon lifetimes of
TiN, C, Si and Ge, the dressed elastic G_BTE and the no-fit three-phonon channel G_inel for
TiN/diamond, TiN/Si and TiN/Ge (001), and G_total = G_BTE + G_inel. TiN/SiC is left out: the
polytype and orientation are only in the paper's Supplementary Information, which was not read.
The orientation (001) for all interfaces is our assumption.

A series electron-phonon term inside TiN (sqrt(g Lambda_p), g = 1e18 W m^-3 K^-1 as quoted by Khan et
al., Lambda_p = the model's TiN phonon conductivity) is reported, **not graded**: g is known only to
an order of magnitude, and the phonon conductivity of a metal is limited by electrons, which our
model leaves out.

Grades (`tbc_tin.grade_t` on the primary model; T_high = 900 K for TiN/diamond, 700 K for TiN/Si
and TiN/Ge, the highest 100 K steps inside each measured range):

| | test | PASS if | prior |
|---|---|---|---|
| T0 | sanity | top phonon frequency TiN 14-23, C 35-45, Si 13.5-17.5, Ge 7.5-10.5 THz; imaginary fraction <= 0.02; all pairs computed at 300 K and T_high | 80 % |
| T1 | the elastic part cannot do it | G_BTE(T_high)/G_BTE(300) below the measured ratio (2.027, 1.807, 1.479) for all three pairs | 85 % |
| T2 | level | G_total within 25 % of the data at both 300 K and T_high, for at least 2 of 3 pairs | 30 % |
| T3 | **rise (main)** | G_total(T_high)/G_total(300) within 15 % of the measured ratio, for at least 2 of 3 pairs | 30 % |
| T4 | beats DMM | at T_high, G_total closer to the data than the model's own DMM, all pairs | 70 % |
| T5 | ordering | G_total(300): TiN/Si > TiN/diamond > TiN/Ge (507 > 449 > 190) | 50 % |
| T6 | soft side | G_total(TiN/Ge) / G_total(TiN/Si) at 300 K in [0.25, 0.55] (measured 0.375) | 55 % |

**Reading rules.**
- T1 PASS and T3 PASS: a no-fit three-phonon interface channel reproduces a factor-two rise to
  900 K. This would be the first such prediction we know of, and the basis of an entry.
- T1 PASS, T3 FAIL with the model rising too little: the channel is too weak at high temperature.
  The missing part then points to higher-order or electron-mediated routes, and stage 2b's
  T-proportional form is disfavoured for these pairs.
- T1 FAIL: our dressed elastic model alone rises as much as the data. Khan et al.'s inference from
  the rise would then need a fully spectral non-equilibrium treatment before it can stand. That would
  be a notable negative result and would be reported as such.
- T2 without T3 (right level, wrong trend), or the reverse, is reported as it is. The challenge needs
  both.

## 4. Limits

- Figure readings (+-5-30) and curve readings (+-2-3 %); supplementary details not read.
- The interface orientation and the ideal, abrupt, epitaxial (001) termination are assumptions.
- The universal potentials are not fitted to TiN phonons; T0 only checks the top frequency.
- The electron-phonon term is a sensitivity, not a prediction. Any direct electron-phonon coupling
  across the interface is not modelled (Khan et al. estimate about 0.1 GW m^-2 K^-1 for TiN/Si).
- Section 2 is exploratory and not graded.
