# PLAN amendment 0.5.4: two reproduction checks pass; stage 1 ran out of memory; a reference error found; tbc-base-v2

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-19. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

Three of the queued jobs have returned. Two are reproduction checks of results already published, and
both pass. Stage 1 (`tbc-base-v1`) was killed at its 3 GB memory cap after 3 of its 10 materials, so
it grades nothing. Its first three phonon calculations still show two things worth recording before
anything is re-run. This file is written **before** `tbc-base-v2` exists. Earlier expectations
are unchanged.

## 1. G0 (job tbc-gep-v1): PASS

Pre-registered (amendment 0.3.2): tests pass, the T^-1/2 phonon-side values equal the 0.3.2 tables to
0.1 MW m^-2 K^-1, and both are convex. The job's 9 tests pass. The job gives (MW m^-2 K^-1, 300-700 K):

| pair | job | amendment 0.3.2 | convex |
|---|---|---|---|
| Al/Si | 512.5 / 586.0 / 738.1 / 892.8 / 1194.9 | 512.5 / 586.0 / 738.1 / 892.8 / 1194.9 | True |
| Al/GaN | 594.4 / 657.4 / 743.6 / 896.3 / 1279.4 | 594.4 / 657.4 / 743.6 / 896.3 / 1279.4 | True |

Output: `results/gep/gep-li2022.json`, `results/gep/tests.txt`.

## 2. H0 (job tbc-direct-v1): PASS

Pre-registered (amendment 0.3.3): tests pass, and the `TBD:` line matches the 0.3.3 tables to 0.2
MW m^-2 K^-1. All 16 tests pass. Every tabulated number matches; the largest difference is 0.045
(DSF-size channel, Al/GaN at 300 K: 402.045 against 402.0). Examples: Sergeev sigma_e(300 K) =
177.5 / 293.5 / 352.8 / 264 (four versions); bounds s < 0.7571 (Al/Si) and 0.7800 (Al/GaN); with
the 210 MW floor s <= 0.506 and sigma_e(300 K) <= 151.8; D1 and D2 thresholds G_min(600 K) = 446.8
and 395.3; exact circuit 293.235 against the 0.8 %-off limit 290.964. The verdicts of amendment 0.4
(Sergeev's bulk-only channel excluded) were graded by hand from the same module. They now also
stand on an independent run of the published code in a clean container.

Output: `results/direct/direct-li2022.json`, `results/direct/run.txt`, `results/direct/tests.txt`.

## 3. tbc-base-v1: killed at 3 GB (no grade)

The 3 setup steps, the compile step and the tests passed. The run step was killed (exit 137, cap
3 GB) 74 s in, while MACE-MPA-0 (float64) was working on AlN, the fourth of ten materials. The
supercells are 108-216 atoms, and MACE's 6 A cutoff makes each force call large. Diamond, with about 160
neighbours per atom inside 6 A against about 80 for GaN, would need roughly twice GaN's memory. So the
failure is the job's memory plan, not the physics code. B1-B3 and E1-E4 are **not graded**.

What the three finished materials show (information, from `results/base-v1/progress.txt`):

| material | relaxed lattice (A) | vs experiment | highest phonon (THz) | PLAN.md reference | ratio |
|---|---|---|---|---|---|
| Al | a 4.0387 | -0.27 % | 7.93 | 9.7 (neutron) | 0.818 |
| Si | a 5.4674 | +0.67 % | 12.68 | 15.6 (Raman) | 0.813 |
| GaN | a 3.2047, c 5.2503 | +0.49 %, +1.26 % | 21.10 | 16.8 = E1(TO) | 1.256 |

No imaginary modes, and every relaxation converged.

**(a) MACE-MPA-0 is about 18 % soft for Al and Si.** That is beyond B2's +-15 % window, so
**B2 for MACE-MPA-0 fails on Al and Si whatever else happens**. This is the systematic softening of
universal potentials reported in the literature (e.g. Deng et al., arXiv:2405.07105), not a bug.
Softer modes mean lower group velocities. In the classical limit the radiation limit scales with the
sound speed, so MACE-based DMM values for Al/Si will come out low. B3 (MACE and orb-v3 agree on
Al/Si DMM within 15 %) is the pre-registered check on exactly this.

**(b) PLAN.md's GaN reference is wrong.** Without LO-TO splitting (none is applied), the highest
zone-centre mode of wurtzite GaN is not E1(TO) (16.8 THz). It is the silent B1(high) mode, measured
at 85.5 meV = **20.67 THz** by inelastic x-ray scattering (Ruf et al., Phys. Rev. Lett. 86, 906
(2001)). MACE's 21.10 THz is within 2.1 % of it. B2 stays graded with the reference fixed in PLAN.md
(we do not move a pre-registered threshold after seeing output). Every report will also print the
B1(high) ratio, labelled "information". AlN's reference (E1(TO), 20.1 THz) has the same flaw. We
have read no measured B1(high) for AlN, so it is left as written and flagged here: AlN may read high
for the same reason.

## 4. tbc-base-v2: pre-registered before it runs

New driver `run_tbc_safe.py` 0.1.0, with tests `test_run_tbc_safe.py`. Each (potential, material) phonon
calculation runs in its own spawned child process, so memory is freed after each material, and an
out-of-memory kill loses one material, not the job. If a MACE child is killed by a signal, that
material is retried once in float32, and the switch is recorded in `info[name]['dtype']`. orb-v3 is
already float32 and gets one attempt. A child that raises an exception is not retried. The physics,
inputs, tables and grades are those of `run_tbc.py`, `tbc_phonons.py` and `tbc_transport.py`,
fetched unchanged from this repository's main branch. The cap rises to 5 GB. Last line `TBC:`; the
line before it, `SAFE:`, gives per potential the materials done, any float32 fallbacks, the peak
child memory and the GaN B1(high) ratio.

Expectations (graded when it returns):
- **M0** (engineering): tests pass and at least 9 of 10 materials finish for each potential. 80 %.
- **B1-B3, E1-E4**: exactly as in PLAN.md. Added before output: **B2 for MACE-MPA-0 fails** (Al, Si
  outside +-15 %; already certain from section 3). B2 for orb-v3: no expectation stated.
- **N1** (new, information): MACE-MPA-0's oriented full-band DMM for Al/Si at 300 K is lower than
  orb-v3's. 65 %. If true and B3 fails, stage 2/3 will use orb-v3 for the phonon side, or the
  material-wise rescaling question gets its own pre-registration. Either choice is made before any
  interface comparison is seen.

## 5. The other queued jobs

`tbc-anh-v1`, `tbc-zno-v1` and `tbc-tin-v1` are still queued. They use the same second-order
supercells (plus phono3py) under the same 3 GB cap, so they may be killed the same way. They are not
re-filed for a wait. If one is killed, it will be re-filed with a per-material child-process driver
like this one, and its expectations will stay as pre-registered.

## Files

`run_tbc_safe.py`, `test_run_tbc_safe.py` (new). `results/gep/`, `results/direct/` (job outputs),
`results/base-v1/progress.txt` and `tests.txt` (the killed job's partial log).
