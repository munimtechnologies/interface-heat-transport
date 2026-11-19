# interface-heat-transport

Open work by MV2 of Munim, Inc. toward a cheap, transferable model of thermal boundary conductance
G(T) across solid-solid interfaces, with no conductance fitted to any specimen. It targets the
public benchmarks of Discovered Materials' Grand Challenge 05 ("Predicting Heat Transport across
Solid Interfaces"). Made with AI agents working in MV2's name; not yet checked by a
thermal-transport specialist.

- `PLAN.md`: the question, the benchmarks and the expectations written before each job runs.
  Dated amendments (`PLAN-amendment-0.2.md`, `0.3`, `0.3.1`, `0.3.2`, `0.3.3`, `0.4`, `0.5`, `0.5.1`, ...)
  add stages and grades without changing earlier text.
- `tbc_transport.py`: elastic baselines (full-band and Debye DMM, Debye AMM, gray Simons/Zeng-Chen
  correction, phonon radiation limit, electron-phonon series resistance). numpy only.
- `tbc_phonons.py`: bulk phonons from open universal potentials (MACE-MPA-0, MIT; orb-v3, Apache-2.0)
  through phonopy.
- `run_tbc.py`: stage-1 job; results go to `results/` once it returns.
- `tbc_bte.py`: stage 2a. A multiband two-stream Boltzmann model of the interface region with
  mode-resolved bulk lifetimes on both sides of an elastic (DMM) interface. In the gray limit it
  reduces exactly to the Simons/Zeng-Chen corrected DMM. numpy only. Its closest published relative is
  the McKelvey-Shockley flux model of a finite Si/Ge junction with first-principles three-phonon
  rates (arXiv:1808.10713). The conductance it returns is the "dressed" one, referred to local
  temperatures as in Shi, Yang, Fisher and Ruan's non-equilibrium Landauer approach (arXiv:1812.07910).
- `tbc_anh.py`, `run_tbc_anh.py`: stage-2a job. It computes three-phonon and isotope lifetimes
  (phono3py) of Al, Si, GaN, AlN and diamond from the same open potentials, then G(T) from
  `tbc_bte` for Al/Si, Al/GaN and AlN/GaN.
- `tbc_inel.py`, `run_tbc_inel.py`: stage 2b (PLAN-amendment-0.3.md). An inelastic channel across
  the interface whose strength comes only from bulk three-phonon rates, densities of states and the
  bond geometry of the terminating layers. In the classical limit it grows in proportion to T. It runs
  on stage 2a's published lifetimes and is graded against the measured slopes of Li et al. (2022).
  numpy only. **Prior art:** the mechanism was proposed by Duda, Norris and Hopkins (J. Heat Transfer
  133, 074501, 2011): three-phonon scattering of phonons crossing the interface with phonons on the
  other side, which gives G = b + cT. `tbc_inel` is a mode-resolved, first-principles version of it
  with nothing fitted to an interface.
- `tbc_gep.py`, `run_tbc_gep.py`: the electron side of a metal/non-metal interface
  (G_ep = sqrt(g Lambda_p) with its temperature dependence from bulk inputs), and the phonon-side
  conductance a measurement implies (PLAN-amendment-0.3.2.md). Pure Python.
- `tbc_direct.py`, `run_tbc_direct.py`: a direct electron-to-substrate-phonon channel in parallel
  (the two-temperature circuit of Lu, Wang and Ruan, PRB 93, 064302, 2016, re-derived in closed
  form), Sergeev's no-fit model of that channel, and the bounds Li et al.'s data put on it
  (PLAN-amendment-0.3.3.md). Pure Python.
- `tbc_grade_d.py`: grades the blind test D1-D4 of amendment 0.3.3 (PLAN-amendment-0.4.md); the
  grade is in `results/dtest/grade-d1d4.json`. Pure Python.
- `tbc_nonmetal.py`, `run_tbc_zno.py`: stage 3a (PLAN-amendment-0.5.md), interfaces without a metal.
  `tbc_nonmetal` holds the data-only bound on a three-phonon channel from measured G(T)
  (`results/nonmetal/bounds-exploratory.json`) and the pre-registered grades Z0-Z5.
  `run_tbc_zno` is the job: ZnO, GaN and AlN lifetimes, then elastic G_BTE plus the stage 2b channel
  for ZnO/GaN and AlN/GaN. It imports stages 2a and 2b and copies none of their code.
- `tbc_dressed.py`: amendment 0.5.1. The bound under the published error bars, the exact dressed
  DMM conductance (1/G_d = 1/G_0 - 1/(2R_1) - 1/(2R_2)), the corrected bound for an elastic part that
  falls, and the supplementary grades Z1g, Z6, Z7 of stage 3a (`results/nonmetal/bounds-0.5.1.json`).
  Pure Python.
- `data/meas-li2022.json`: the measured Al/Si and Al/GaN conductances (80-700 K) of Li et al.,
  Nat. Commun. 13, 4901 (2022), read off their Fig. 1 without sight of any model output.
  `data/read-d1d3.json`: Li et al.'s rough Al/Si Sample 2 and Cheng et al.'s epitaxial
  Al/sapphire (Commun. Phys. 2020, arXiv:1906.05484), read blind to the thresholds they were graded
  against. `data/read-nonmetal.json`: ZnO/GaN (Gaskins et al., arXiv:1710.09525), bonded GaN/AlN
  (UCLA, APL 2024) and Walwil et al.'s GaN/AlN, AlN/SiC and GaN/SiC (JAP 2025), read blind to every
  model result. `data/read-gaskins2018.json`: the published ZnO/GaN figure (Nano Lett. 2018) with
  its error bars and how they were obtained.
- `readings/`: notes from the papers this work builds on.
- Unit tests: `python -m unittest -v test_tbc_transport test_tbc_bte test_tbc_inel test_tbc_gep test_tbc_direct test_tbc_grade_d test_tbc_nonmetal test_tbc_dressed`.

Status (2026-11-19): one graded result so far (amendment 0.4, below). Stages 1, 2a, 2b and 3a are
pre-registered and their jobs are queued, so there are no model results from them yet.

**Amendment 0.5.1 (correction).** The ZnO/GaN bound of amendment 0.5 does not hold up. The
published error bars (95 % confidence contours of a separate fit at each temperature) allow a
three-phonon channel of up to 229 MW m^-2 K^-1 (47 % of G) at 300 K if the errors are independent;
16.5-19.5 holds only if they are fully shared. The bound also assumed that the elastic part cannot
fall with temperature. That fails for the conductance a TDTR fit reports: for a DMM interface it is
exactly 1/G_d = 1/G_0 - 1/(2R_1) - 1/(2R_2), which falls when the stiffer side's radiation limit
keeps growing after the softer side has saturated (as Shi, Yang, Fisher and Ruan, arXiv:1812.07910,
report for ZnO/GaN). Their non-equilibrium treatment also anticipated our reading that the gap above
AGF is mostly elastic. Three supplementary grades (Z1g, Z6, Z7) now test this in stage 3a, filed
before any output.

**Amendment 0.5 (exploratory bound + stage 3a pre-registered; bound corrected in 0.5.1).**
Interfaces between two non-metals have no electron route. Measured ZnO/GaN conductance is flat from
300 to 500 K (491, 520, 502 MW m^-2 K^-1, read off Gaskins et al.'s figure). With a non-falling
elastic part and the reading error only, a three-phonon channel could carry at most 16.5 (central)
to 76 MW m^-2 K^-1 at 300 K. Job `tbc-zno-v1` tests our no-fit inelastic channel against this
(Z1) and the elastic model against the measured levels (Z3, Z4).

**Amendment 0.4 (graded, pre-registered blind test).** Sergeev's model of a direct route from a
metal's electrons into the substrate's phonons uses only bulk properties of the metal and grows in
proportion to T. Even in its smallest normalisation it requires Al interfaces to gain at least 35 %
(400 -> 600 K) and 29 % (300 -> 480 K) in conductance on any substrate. Measured, the same Al gains 5.5 % on a
rough Si interface (Li et al.'s Sample 2) and 3 % on epitaxial sapphire (Cheng et al.). The model is
excluded at both (D1, D3), in all four normalisations, and the verdicts survive the full reading
error. A T-linear channel of the size De San Feliciano et al. inferred (107 MW m^-2 K^-1 at 300 K)
survives the pre-registered margin; we had expected it to fail (D2). Exploratory: at central
readings, an interface-independent T-linear channel is at most 44.5 MW m^-2 K^-1 at 300 K. Take that
much out of the sharp Al/Si sample and its phonon side still rises 65 % from 400 to 700 K. So the
rise Li et al. measured needs a mechanism that depends on the interface itself.

Amendment 0.3.1 records the measured inputs and fixes a solver error before any model output
existed. With the electron-phonon resistance in series, the fitted rise is not monotonic in the
channel strength, and 0.1.0 would have read "too weak" whatever the model gave. The amendment also
works out what the measurements alone imply. If the rise above 400 K comes from any process that
grows in proportion to T on top of a saturated elastic part, that process already carries at least
43 % (Al/Si) and 35 % (Al/GaN) of the phonon-side conductance at 400 K.

Amendment 0.3.2 adds the electron side. Bulk inputs (Al's lattice conductivity is three-phonon
limited; its electron-phonon coupling barely depends on T) give a series conductance falling as
T^-1/2. Removing it from the measurements leaves a phonon-side conductance that is convex from 300
to 700 K at both Al/Si and Al/GaN. A saturated elastic part plus a classical three-phonon channel
cannot give that shape. The inference amplifies the error bars about 2.3x at 700 K, so this is not
decisive. Li et al.'s own elastic (DMM) curve, put in series with this G_ep, misses about 40 % of
the measured G at 700 K for both pairs.

Amendment 0.3.3 adds a direct route from the metal's electrons into the substrate's phonons, in
parallel with the usual path. It also found that a T-linear direct channel of published size makes
Li et al.'s data degenerate between that route and inelastic phonons, and it pre-registered the
blind test that amendment 0.4 grades.

Stage 2a tests a scaling argument (PLAN-amendment-0.2.md, section 2). For semi-infinite solids
and an elastic interface, multiplying every bulk scattering rate by the same factor leaves the
interface conductance unchanged. In the classical limit three-phonon rates all grow in proportion
to T, so bulk anharmonic thermalization next to the interface cannot by itself produce the linear
rise of G with T measured at Al/GaN. If the job confirms this for real materials, a transferable
model needs an inelastic channel at the interface itself (as in Duda, Norris and Hopkins 2011, and
stage 2b).

Licence: Apache-2.0 (code), CC BY 4.0 (text). Contact: MV2, Munim, Inc., mv2@munimtech.com.
