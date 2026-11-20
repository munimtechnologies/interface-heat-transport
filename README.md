# interface-heat-transport

An open, pre-registered attempt to **predict the thermal boundary conductance G(T) of solid-solid
interfaces without fitting to any specimen**, aimed at Discovered Materials' Grand Challenge 05
(interfacial heat transport). By MV2 of Munim, Inc. (mv2@munimtech.com). The work is made with AI
agents working in MV2's name and has **not yet been checked by a thermal-transport specialist**.
Code Apache-2.0; text and data CC BY 4.0.

## Newest result (2026-11-20): TiN's large deficit below its own DMM was mostly a binning artefact

[PLAN-amendment-0.5.11.md](PLAN-amendment-0.5.11.md) grades job `tbc-fill-v1` against rules
written before it ran ([amendment 0.5.8](PLAN-amendment-0.5.8.md) sec. 4). Our multiband
Boltzmann model had put TiN on diamond, Si and Ge 5-6 times below their own diffuse-mismatch (DMM)
value. The cause was in how frequency bins that the lifetime mesh never sampled were given a
phonon mean free path: in TiN those invented lengths carried two thirds of the model's
conductivity. Making the low-frequency lengths from real conductivities raises G 3.2-6.1 times
(the main test F2 passed); TiN/Si then sits at its DMM (1.00-1.05) and TiN/C and TiN/Ge at
0.62-0.72 of it. Two tests failed: the unsampled bins reach 3.2 THz, not only the gap below the
mesh's lowest phonon (F1), and where the fix stops early they still inflate TiN's conductivity by
13-22 % (F4). The "5-6x below DMM" reading of amendment 0.5.7 is retired. The next stage, a
conserving band rule with no invented lengths (tests V0-V5), is pre-registered in sec. 4.
Raw output: [results/tin/fill/](results/tin/fill/).

## Results so far (each graded against expectations written before its job ran)

- **0.4**: Sergeev's direct electron-to-substrate channel (bulk-only) is ruled out at rough Al/Si
  and epitaxial Al/sapphire; the measured rise of G with temperature needs a mechanism that
  depends on the interface itself.
- **0.5.1**: the ZnO/GaN bound was corrected: the error bars allow a much larger three-phonon
  channel, and the reported conductance can fall with temperature.
- **0.5.2**: for TiN/diamond, the dressed elastic formula explains the room-temperature value but
  not the factor-two rise to about 1000 K.
- **0.5.9**: no-fit values at 300 K: Al/Si 422 against 379 measured, Al/GaN 365 against 423,
  AlN/GaN 304 against 324 MW m^-2 K^-1.
- **0.5.10** (job `tbc-base-v2`): elastic models cannot give the measured rise of G with
  temperature (+2-3 % from 400 to 700 K; measured Al/Si rises from 379 to 530 MW m^-2 K^-1
  between 300 and 700 K); with Al's
  electron-phonon term the non-equilibrium correction brackets Al/Si at 300 K; the elastic DMM
  overstates a TiC interlayer's gain 4-6 times; orb-v3 is the primary phonon potential for new
  stages (top frequencies within 6.1 %).
- **0.5.8 / 0.5.11**: the binning defect above, found while grading and confirmed by a
  pre-registered job.

## Contents

- [PLAN.md](PLAN.md) and its dated amendments (0.2 to 0.5.11): the pre-registered expectations,
  grades and reading rules.
- Modules: `tbc_transport` (elastic baselines), `tbc_phonons` (phonons from open universal
  potentials), `tbc_bte` (multiband Boltzmann interface model), `tbc_anh` (anharmonic lifetimes),
  `tbc_inel` (inelastic three-phonon channel after Duda, Norris and Hopkins 2011), `tbc_gep`
  (electron-phonon coupling), `tbc_direct` (direct electron-to-substrate channel), `tbc_nonmetal`
  and `tbc_dressed` (non-metal interfaces), `tbc_tin` (TiN at high temperature), `tbc_diag` and
  `tbc_fill` (the diagnostics of amendments 0.5.7-0.5.11).
- `results/`: raw job outputs, unedited, with our grades.
- Measured data were digitised from published figures by AI reading agents, blind to the model.

## Status (2026-11-20)

Graded: amendments 0.4, 0.5.1, 0.5.2, 0.5.6-0.5.11. Queued: `tbc-anh-v2`. Next: `tbc-cover-v1`
(amendment 0.5.11 sec. 4). The model is not yet a candidate Grand Challenge entry: no frozen
version meets every benchmark.
