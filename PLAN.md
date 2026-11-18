# PLAN: a transferable, no-fit model of heat flow across solid interfaces

MV2 of Munim, Inc. (mv2@munimtech.com). Version 0.1, 2026-11-18. Made with AI agents working in
MV2's name; no thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

## 1. The question

Heat crossing a boundary between two solids meets a thermal boundary resistance R_K = 1/G. In
chips, these interface resistances now take a growing share of the total. The two textbook models,
the acoustic mismatch model (AMM, Little 1959) and the diffuse mismatch model (DMM, Swartz and Pohl
1989), use only bulk phonon properties and are known to miss measurements, often by factors of 2-10.
Atomistic simulations (MD, Green's functions) do better for clean crystalline interfaces but are too
slow for routine use and impractical for polycrystalline, amorphous or polymer interfaces.

**Question**: can a computationally cheap framework, with no conductance fitted to any specimen,
predict G(T) across interface classes (epitaxial metal/semiconductor, metal/interlayer/diamond,
metal/polymer) better than AMM and DMM, and if so, which physics does the work?

This is Grand Challenge 05 of Discovered Materials (Matforge Inc.), "Predicting Heat Transport
across Solid Interfaces" (https://discoveredmaterials.com/grand-challenges/interfacial-heat-transport/),
prize up to $100,000. Its public benchmarks (G in MW m^-2 K^-1):

1. Al(111)/Si(111) and Al(111)/GaN(0001), 80-700 K: 379 and 423 at 300 K, a roughly linear rise
   above about 400 K (Li et al., Nat. Commun. 13, 4901 (2022)). Li et al. attribute the rise to
   inelastic phonon transport; their less sharp Al/Si sample saturates instead.
2. Cu/diamond with a 10 nm crystalline TiC interlayer: G rises by a factor 1.48 (Chang et al., ACS
   Appl. Mater. Interfaces 11, 26507 (2019)).
3. Au/PMMA/Si stacks with 2 nm adhesion layers: boundary terms 59 (Au), 115 (Au/Ti), 139 (Au/Ni)
   and the PMMA-thickness dependence (Sandell et al., Nanomaterials 10, 670 (2020)).
A model frozen on these is then run by the organisers on hidden, industrially relevant interfaces
(AlN/GaN is one example) against AMM/DMM, with atomistic results as ground truth.

**The wow, if it works**: a cheap formula-level model that turns bulk phonon data into interface
conductances for materials it was never tuned on, a long-standing open problem in thermal
transport. An honest negative result (which physics a cheap model cannot capture) also counts.

## 2. Stage 1 (this job, tbc-base-v1): honest elastic baselines

Before building anything new, we measure exactly how far the standard elastic models are from the
three benchmarks when they are given the best bulk phonons an open, cheap pipeline can make:

- Bulk phonons of Al, Si, GaN, AlN, Cu, diamond, TiC, Au, Ti, Ni from two open universal
  machine-learned potentials (MACE-MPA-0 medium, MIT; orb-v3-conservative-inf-omat, Apache-2.0)
  with phonopy finite displacements, after a symmetry-constrained relaxation. No LO-TO correction.
- Spectral phonon flux M(nu) toward the interface, isotropic (polycrystal) and along the interface
  normal (Al and Si [111], GaN and AlN [0001]).
- Full-band DMM, the same with the gray Simons/Zeng-Chen non-equilibrium correction, phonon
  radiation limits, Debye AMM and Debye DMM; electron-phonon resistance in the metal in series,
  G_ep = sqrt(g Lambda_p).
- Definitions and formulas: `tbc_transport.py` docstring. Unit tests: `test_tbc_transport.py`.

## 3. Expectations, written before the job ran

Graded from `tbc-summary.json` / the `TBC:` line of run.txt. Primary model MACE-MPA-0; orb-v3 is
the cross-check (a criterion where the two disagree is reported as SPLIT).

- **B1** all unit tests pass.
- **B2** for each model: every material relaxes and phonons, fewer than 1 % of mesh modes below
  -0.05 THz, and the highest frequency within 15 % of the reference (Al 9.7, Cu 7.3, Au 4.6,
  Ni 8.9, Si 15.6, C 39.9 THz; GaN 16.8 and AlN 20.1 THz, highest TO mode).
- **B3** the two models' isotropic full-band DMM for Al/Si at 300 K agree within 15 %.
- **E1** (elastic models saturate) oriented full-band DMM G(700 K)/G(400 K) < 1.10 for both Al/Si
  and Al/GaN. The measurement rises linearly over this range. Expected: PASS (about 95 %).
- **E2** (the measurement sits between the two elastic readings) for Al/Si at 300 K, with the Al
  electron-phonon term in series (g = 2.45e17 W m^-3 K^-1, Lambda_p = 6 W m^-1 K^-1, G_ep about
  1.2 GW m^-2 K^-1): conventional DMM < 379 < Simons-corrected DMM. Expected: PASS (about 55 %).
- **E3** (elastic DMM gets the material ordering wrong) oriented full-band DMM gives
  G(Al/GaN)/G(Al/Si) < 1.12 at 300 K (measured 423/379 = 1.116). Reason: GaN has a gap between its
  Ga-dominated acoustic modes (up to about 8 THz) and its optical modes, so in an elastic model the
  Al phonons between about 8 and 10 THz find few partners in GaN but many in Si. Expected: PASS
  (about 70 %).
- **E4** (elastic DMM overstates the TiC gain) Cu/TiC/diamond over Cu/diamond, isotropic, 300 K,
  with the Cu electron-phonon term (g = 1.0e17, Lambda_p = 2.6) on the Cu/diamond side, a 4 GW
  electronic Cu/TiC contact, 10 nm TiC at 20 W m^-1 K^-1 and TiC's own electron-phonon term
  neglected: ratio > 1.48 for both conventional and corrected DMM. Expected: PASS (about 65 %).

The electron-phonon inputs are approximate literature values for bulk metals, not fitted to any
interface; the job also reports them at twice the coupling.

## 4. What each outcome means for stage 2

- E1 PASS: confirms that an elastic model cannot give the rise, so stage 2 adds an inelastic
  channel. The first candidate is a two-phonon (three-phonon at the interface) channel whose
  strength comes from bulk anharmonicity (Grueneisen parameters from the same potentials), not
  from a fit.
- E2 PASS: the non-equilibrium correction alone brackets the data; the frozen model must then
  decide between the two from physics (sharp vs diffuse interface), not by picking the closer one.
  E2 FAIL low (both below 379): elastic channels are too small at 300 K already.
- E3 PASS: a clear, testable target for the inelastic channel (it must lift Al/GaN more than Al/Si).
- E4: decides whether the interlayer benchmark needs electron-phonon coupling inside TiC or bonding
  physics at the TiC/diamond side.
- B2 FAIL on a material: that material's phonons go to the other model or to DFT data before any
  conductance from it is used.

Polymer stacks (benchmark 3) need a different treatment (weak van der Waals bonding, amorphous
PMMA) and start in stage 3.

## 5. Rules for ourselves

- No interface conductance is fitted to any benchmark specimen, ever; every parameter is a bulk
  property or comes from a calculation. Benchmarks are used only to grade.
- Every new model gets expectations in this file (or a dated amendment) before its job runs.
- Entering the challenge needs Munim, Inc.'s named representative to accept the official rules,
  which grant the sponsor a broad licence to every submission. Nothing is submitted without that.
