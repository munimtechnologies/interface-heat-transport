# interface-heat-transport

Open work by MV2 of Munim, Inc. toward a cheap, transferable model of thermal boundary conductance
G(T) across solid-solid interfaces, with no conductance fitted to any specimen. It targets the
public benchmarks of Discovered Materials' Grand Challenge 05 ("Predicting Heat Transport across
Solid Interfaces"). Made with AI agents working in MV2's name; not yet checked by a
thermal-transport specialist.

- `PLAN.md`: the question, the benchmarks and the expectations written before each job runs.
  Dated amendments (`PLAN-amendment-0.2.md`, `PLAN-amendment-0.3.md`, `PLAN-amendment-0.3.1.md`,
  `PLAN-amendment-0.3.2.md`, ...) add stages without changing earlier text.
- `tbc_transport.py`: elastic baselines (full-band and Debye DMM, Debye AMM, gray Simons/Zeng-Chen
  correction, phonon radiation limit, electron-phonon series resistance). numpy only.
- `tbc_phonons.py`: bulk phonons from open universal potentials (MACE-MPA-0, MIT; orb-v3, Apache-2.0)
  through phonopy.
- `run_tbc.py`: stage-1 job; results go to `results/` once it returns.
- `tbc_bte.py`: stage 2a. A multiband two-stream Boltzmann model of the interface region with
  mode-resolved bulk lifetimes on both sides of an elastic (DMM) interface. In the gray limit it
  reduces exactly to the Simons/Zeng-Chen corrected DMM. numpy only. Its closest published relative is
  the McKelvey-Shockley flux model of a finite Si/Ge junction with first-principles three-phonon
  rates (arXiv:1808.10713).
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
- `data/meas-li2022.json`: the measured Al/Si and Al/GaN conductances (80-700 K) of Li et al.,
  Nat. Commun. 13, 4901 (2022), read off their Fig. 1 without sight of any model output.
- `readings/`: notes from the papers this work builds on.
- Unit tests: `python -m unittest -v test_tbc_transport test_tbc_bte test_tbc_inel test_tbc_gep`.

Status (2026-11-19): stages 1, 2a and 2b are pre-registered. The stage 1 and 2a jobs are queued.
Stage 2b runs once stage 2a's lifetimes are published. No model results yet.

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
decisive. If it is real, it points to a superlinear phonon channel or to a direct electron route
into the substrate's phonons, which recent first-principles work (arXiv:2510.27499) puts at about a
third of the Al/Si conductance. Li et al.'s own elastic (DMM) curve, put in series with this G_ep,
misses about 40 % of the measured G at 700 K for both pairs.

Stage 2a tests a scaling argument (PLAN-amendment-0.2.md, section 2). For semi-infinite solids
and an elastic interface, multiplying every bulk scattering rate by the same factor leaves the
interface conductance unchanged. In the classical limit three-phonon rates all grow in proportion
to T, so bulk anharmonic thermalization next to the interface cannot by itself produce the linear
rise of G with T measured at Al/GaN. If the job confirms this for real materials, a transferable
model needs an inelastic channel at the interface itself (as in Duda, Norris and Hopkins 2011, and
stage 2b). Stage 2b builds that channel from bulk data and asks whether it is strong enough. If it
is far too weak, the rise needs either much stronger anharmonicity at the interface or a direct
electron route in metals.

Licence: Apache-2.0 (code), CC BY 4.0 (text). Contact: MV2, Munim, Inc., mv2@munimtech.com.
