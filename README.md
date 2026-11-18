# interface-heat-transport

Open work by MV2 of Munim, Inc. toward a cheap, transferable model of thermal boundary conductance
G(T) across solid-solid interfaces, with no conductance fitted to any specimen. It targets the
public benchmarks of Discovered Materials' Grand Challenge 05 ("Predicting Heat Transport across
Solid Interfaces"). Made with AI agents working in MV2's name; not yet checked by a
thermal-transport specialist.

- `PLAN.md`: the question, the benchmarks and the expectations written before each job runs.
  Dated amendments (`PLAN-amendment-0.2.md`, `PLAN-amendment-0.3.md`, ...) add stages without changing
  earlier text.
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
  numpy only.
- Unit tests: `python -m unittest -v test_tbc_transport test_tbc_bte test_tbc_inel`.

Status (2026-11-18): stages 1, 2a and 2b are pre-registered. The stage 1 and 2a jobs are queued.
Stage 2b runs once stage 2a's lifetimes are published and the measured slopes have been read off
the figure. No results yet.

Stage 2a tests a scaling argument (PLAN-amendment-0.2.md, section 2). For semi-infinite solids
and an elastic interface, multiplying every bulk scattering rate by the same factor leaves the
interface conductance unchanged. In the classical limit three-phonon rates all grow in proportion
to T, so bulk anharmonic thermalization next to the interface cannot by itself produce the linear
rise of G with T measured at Al/GaN. If the job confirms this for real materials, a transferable
model needs an inelastic channel at the interface itself. Stage 2b builds that channel from bulk
data and asks whether it is strong enough. If it is far too weak, the rise needs either much
stronger anharmonicity at the interface or a direct electron route in metals.

Licence: Apache-2.0 (code), CC BY 4.0 (text). Contact: MV2, Munim, Inc., mv2@munimtech.com.
