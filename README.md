# interface-heat-transport

Open work by MV2 of Munim, Inc. toward a cheap, transferable model of thermal boundary conductance
G(T) across solid-solid interfaces, with no conductance fitted to any specimen. It targets the
public benchmarks of Discovered Materials' Grand Challenge 05 ("Predicting Heat Transport across
Solid Interfaces"). Made with AI agents working in MV2's name; not yet checked by a
thermal-transport specialist.

- `PLAN.md`: the question, the benchmarks and the expectations written before each job runs.
- `tbc_transport.py`: elastic baselines (full-band and Debye DMM, Debye AMM, gray Simons/Zeng-Chen
  correction, phonon radiation limit, electron-phonon series resistance). numpy only.
- `tbc_phonons.py`: bulk phonons from open universal potentials (MACE-MPA-0, MIT; orb-v3, Apache-2.0)
  through phonopy.
- `run_tbc.py`: stage-1 job; results go to `results/` once it returns.
- `test_tbc_transport.py`: unit tests (`python -m unittest -v test_tbc_transport`).

Status (2026-11-18): stage 1 (baselines) pre-registered, job queued. No results yet.

Licence: Apache-2.0 (code), CC BY 4.0 (text). Contact: MV2, Munim, Inc., mv2@munimtech.com.
