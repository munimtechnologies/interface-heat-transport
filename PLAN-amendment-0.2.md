# PLAN amendment 0.2: stage 2a, bulk anharmonicity next to the interface

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-18. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.
Written before job `tbc-anh-v1` ran and before stage 1 (`tbc-base-v1`) returned. Nothing in
PLAN.md 0.1 changes; this adds stage 2a.

## 1. Why this stage

The measured Al/GaN conductance rises roughly linearly above about 400 K (Li et al. 2022), and an
elastic interface model cannot do that once the heat capacities saturate (stage 1, E1). Two
families of explanations exist:

1. **Inelastic processes at the interface**: phonons of different frequencies combine or split as
   they cross (Hopkins's higher-harmonic and anharmonic inelastic models, 2009-2011; spectral
   NEMD by Feng and co-workers, e.g. PRB 99, 045301 (2019), where inelastic transport carries more
   than half of the conductance).
2. **Bulk anharmonicity next to the interface**: phonon-phonon scattering in the bulk refills the
   modes that cross well with energy from modes that do not (Le et al., PRB 95, 245417 (2017),
   who found in MD that interfacial anharmonicity contributes "a surprisingly small part" and
   bulk anharmonicity most of the rise; Duda et al.'s near-interface scattering DMM).

A cheap transferable model has to know which of the two carries the rise, because only the second
can be built from bulk data alone.

## 2. A scaling argument (to be checked by a specialist)

Take two semi-infinite solids, an elastic interface (each phonon keeps its frequency), and the
linearized Boltzmann transport equation on each side. Multiply the whole collision operator by
a constant lambda. Then x -> x / lambda maps every solution onto a solution with the same heat
flux and the same jump between the two extrapolated bulk temperature profiles. So the interface
conductance G, defined against those profiles, does not change.

In the classical limit every heat capacity is k_B and every three-phonon rate is proportional to T
(each rate is linear in the occupation numbers, which go as k_B T / h nu). The whole collision
operator is then T times a fixed operator, so **G cannot depend on T through bulk three-phonon
thermalization next to an elastic interface**, however strong that thermalization is.
A T-dependence of this kind needs one of: quantum heat capacities, scattering that does not grow
with T (isotopes, defects, boundaries, a finite film), four-phonon and higher processes,
temperature-dependent harmonic properties, or inelastic interface processes (family 1).

This is consistent with Landry and McGaughey's MD (PRB 80, 165304 (2009)): their Si/Ge boundary
resistance is constant at 3.1e-9 m^2 K/W below about 500 K and then falls, which they attribute
to inelastic interface scattering. It is also consistent with finite-domain MD that shows rises
and size effects (e.g. AIP Advances 5, 047102 (2015)), because a finite length breaks the
scaling. It is in tension with readings that credit bulk anharmonicity for the high-temperature
rise at semi-infinite interfaces. We have not found the argument stated
in this form; if it is already in the literature, we will cite it.

The model of section 3 has the property exactly. A unit test checks it (multiplying every length
by 0.1, 3 or 20 changes G by less than 1 %).

## 3. The model (tbc_bte.py)

A multiband two-stream Boltzmann model. Frequency bins of 0.2 THz are the bands. On each side,
band b has a radiation conductance K_b (from the dense phonopy mesh, quantum heat capacity) and a
length l_b = kappa_b / (2 K_b), where kappa_b is the band's share of the bulk conductivity along
the normal from phono3py lifetimes (relaxation-time approximation; three-phonon plus natural
isotope scattering). A stream relaxes toward the local bath temperature over l_b. The heat flux is
conserved, and the interface is diffuse and elastic (DMM, detailed balance per band). G is
referenced to the extrapolated bulk profiles. In the gray limit it reduces exactly to the
Simons/Zeng-Chen corrected DMM (unit test). Limits: isotope scattering is treated as
thermalizing (an upper bound on its mixing effect); Al phonon lifetimes are phonon-phonon only
(electron scattering of Al phonons is omitted, and the electron-phonon resistance is added in
series as in stage 1); no four-phonon scattering; and the third-order force constants use pair
cutoffs of 5 Angstrom (Al) and 4 Angstrom (GaN, AlN).

Potentials: orb-v3-conservative-inf-omat (Apache-2.0) is primary for lifetimes, because its
published Matbench Discovery thermal-conductivity error (kappa_SRME) is lower than MACE-MPA-0's.
MACE-MPA-0 (MIT) is the cross-check. If orb-v3 does not load, MACE-MPA-0 becomes primary and this
is reported.

## 4. Expectations (graded from tba-summary.json / the `TBA:` line)

- **A0** (pipeline): for every material, the sum of our per-bin conductivities (with isotopes)
  reproduces phono3py's own conductivity at 300 K within 3 %. Expected PASS (90 %).
- **A1** (lifetime sanity): phono3py conductivity at 300 K (with isotopes) within Si 100-200,
  diamond 1000-3500, GaN 120-400, AlN 150-450 W m^-1 K^-1. Expected PASS (60 %; universal
  potentials tend to underestimate conductivity).
- **S2** (main): in the oriented geometry with three-phonon plus isotope lifetimes,
  G_BTE(700 K) / G_BTE(400 K) < 1.10 for both Al/Si and Al/GaN. Expected PASS (85 %).
  PASS means that bulk anharmonicity next to an elastic interface does not give the measured rise.
  The frozen model then needs an inelastic interface channel (family 1), which is stage 2b.
  FAIL means the non-scaling ingredients (quantum heat capacities, isotopes) already give a rise,
  and we measure how much of the data they explain before adding anything.
- **S3**: oriented Al/Si at 300 K: conventional DMM < G_BTE < gray-Simons DMM. Expected PASS (65 %).
- **S5**: removing isotope scattering changes G_BTE at 300 K by less than 5 % (Al/Si, Al/GaN,
  oriented). Expected PASS (70 %).
- **S6**: the two potentials' G_BTE for oriented Al/Si at 300 K agree within 15 %. Expected PASS
  (70 %).
- Recorded, not graded: the classical-heat-capacity variant at 400 and 700 K (G should come out
  nearly flat, differing only because phono3py's rates are not exactly proportional to T at these
  temperatures); G_BTE in series with the Al electron-phonon term at 300 K against 379 (Al/Si) and
  423 (Al/GaN); G_BTE for AlN/GaN against the two experiments (320 and 490); the linear-bath
  closed form; and a grid check at half the nodes.

No measured conductance is used anywhere except to print it next to the predictions.
