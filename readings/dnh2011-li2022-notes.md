# Reading notes: Duda, Norris and Hopkins (2011) and Li et al. (2022), plus related prior art

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-19. Made with AI agents working in MV2's name; the
papers were read by an AI reading agent (MV2's hands) on 2026-11-18, without sight of any model
output. Text CC BY 4.0. Quotes are as reported by the reader; plot readings are approximate.

## Duda, Norris and Hopkins, J. Heat Transfer 133, 074501 (2011)

"On the linear temperature dependence of phonon thermal boundary conductance in the classical
limit." A free copy is on the last author's website (patrickehopkins.com).

- **Cause of the linear rise** (quote): the "linear temperature dependence ... is due to three-phonon
  scattering events in which phonons propagating across the interface interact with phonons on the
  other side of the interface". This is anharmonic scattering of crossing phonons with phonons
  across the interface. It is not bulk thermalization next to an elastic interface.
- **Model**: h = h'_S + h'_T, with h'_S = (1/8 pi^2) sum_j integral over k <= k_trans of
  beta hbar omega k^2 (df/dT) (lambda / tau_2) dk. The scattering rate is
  tau^-1 proportional to A_j omega^4 T and A_j omega^2 T, which gives h = b + cT.
- **Molecular dynamics**: none of their own. They compare with Landry and McGaughey's Si/Ge
  Stillinger-Weber MD (PRB 80, 165304, 2009). Domain lengths, mass ratio and length dependence are
  not stated. Varying A_j by +-25 % is "not particularly sensitive".

**What this means for this project.** DNH's channel carries a fixed interfacial length lambda.
Multiplying every scattering rate by s therefore multiplies h'_S by s. Their channel breaks the
scale invariance that PLAN-amendment-0.2 section 2 derives for an elastic interface between
semi-infinite solids. That is why it can give a linear rise. There is no conflict with our
scaling argument. The argument covers only bulk thermalization at an elastic interface, and DNH's
mechanism is an inelastic process within lambda of the interface. **Prior art, stated plainly:**
stage 2b's channel (`tbc_inel.py`: bulk three-phonon rates, densities of states and the bond
geometry of the terminating layers) is a mode-resolved, first-principles version of the mechanism
DNH proposed in 2011. What stage 2b adds is that nothing in it is fitted to an interface.

## Li et al., Nat. Commun. 13, 4901 (2022)

- Temperature dependence (quotes): "interface thermal conductance linearly increases at high
  temperatures"; "Gpp ... has a stronger temperature dependence than G". No quadratic, T^2 or
  four-phonon wording.
- Inelastic share: no number for Al/Si or Al/GaN, only "substantial portion". They cite Landry and
  McGaughey (inelastic dominant above about 500 K) and Feng and Ruan (> 50 % for Si/Ge at room
  temperature).
- G_ep = sqrt(g Lambda_p), with g from Waldecker et al. (PRX 6, 021003, 2016) and Lambda_p from
  Jain and McGaughey (PRB 93, 081206, 2016). Numbers and T-dependence: not stated.
- DMM curve in Fig. 1 (read, approximate): Al/Si 460/475/490 and Al/GaN 440/455/470 MW m^-2 K^-1
  at 300/400/700 K. DMM inputs: Quantum-Espresso DFT, GGA-PBE, 40 Ry, 20x20x20 k-points. The
  reader did not say whether the curve includes the G_ep term; we treat it as phonon-only.
- AGF transmission (quote): "the transmissivity of low-energy phonons (here we define phonons
  lower than 4 THz as low-energy phonons) is higher than the high-energy phonons".

## Related prior art found on 2026-11-19 (abstract level)

- Wang, Lu and Ruan, J. Appl. Phys. 119, 225109 (2016): the first-principles lattice thermal
  conductivity of Al is about 5.8 W m^-1 K^-1 at room temperature with phonon-phonon scattering
  only, and "almost unchanged" when phonon-electron scattering is added. So Lambda_p of Al is
  three-phonon limited and falls roughly as 1/T.
- Waldecker, Bertoni, Ernstorfer and Vorberger, PRX 6, 021003 (2016), arXiv:1507.03743: the
  DFT electron-phonon coupling of Al depends only weakly on electron temperature. Coupling values
  from two-temperature fits are about a factor of two below ab initio values.
- De San Feliciano et al., arXiv:2510.27499 (accepted in PRB): first-principles NEGF for Au/Si,
  Pt/Si and Al/Si. For Al/Si, phonon-phonon transport cannot explain the measured values, and they
  infer a direct interfacial electron-phonon channel of about one third of the total conductance.
  This competes directly with the inelastic phonon reading of the Al/Si rise. From the HTML full
  text (read 2026-11-19): NEGF G_pp for Al/Si about 210 MW m^-2 K^-1 near room temperature, against
  an experimental 350; inferred direct channel G_ep = 107 MW m^-2 K^-1 at 300 K only (no
  temperature dependence reported), placed **in parallel** with the phonon path in a transient
  two-temperature model. Their MD with anharmonicity adds only about 9-12 % over NEGF. Al inputs:
  g = 246e15 W m^-3 K^-1, k_p = 8.95 W m^-1 K^-1, k = 237 W m^-1 K^-1, 100 nm film. Au/Si:
  G_pp 38, measured 43, direct channel 5. Pt/Si: G_pp 65, measured 80-150, direct channel 15-80.
  Their section III.5.4 also describes the coupling as indirect, which conflicts with the "direct"
  wording elsewhere.
