# PLAN amendment 0.3: stage 2b, an inelastic interface channel built from bulk data

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-18. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

This was written before any stage 1 or stage 2a result returned (jobs `tbc-base-v1` and `tbc-anh-v1`
are still queued), before the bulk lifetimes it uses exist, and before the measured temperature
slopes it is graded against have been read off the published figure. Nothing in PLAN.md 0.1 or
amendment 0.2 changes. Results of stages 1 and 2a will be graded in later dated amendments.

## 1. The question

Li et al. (Nat. Commun. 13, 4901 (2022)) measured epitaxial Al/Si and Al/GaN interfaces from 80 to
700 K. G is 379 and 423 MW m^-2 K^-1 at room temperature and keeps rising linearly with T above
about 400 K (Al's Debye temperature is 428 K), where every elastic model flattens. The authors
attribute the rise to inelastic phonon transport but give no predictive model and no fitted
inelastic term. The only models that produce a linear rise without molecular dynamics are the
Hopkins-Duda-Norris family (higher-harmonic and anharmonic inelastic models, 2009-2011). They set
the inelastic transmission from densities of states and occupation numbers, not from how
anharmonic the materials really are.

Stage 2a tests whether bulk anharmonicity next to an elastic interface can produce the rise. Our
scaling argument says it cannot in the classical limit. Stage 2b asks the complementary question:
**if the bonds that cross the interface are exactly as anharmonic as the bulk, how much conductance
do three-phonon processes across them carry, and is it enough to explain the measured slope?**
Every input is a bulk or crystallographic quantity, so the answer is a prediction, not a fit. Both
possible outcomes are informative:

- the right order (scale factor s* between 0.5 and 2, section 4): a transferable, no-fit inelastic
  term exists, which is what Grand Challenge 05 asks for;
- far too weak (s* > 8, beyond any reasonable change of the bond counting below): three-phonon
  exchange across bulk-strength bonds cannot carry the measured rise. Then something else must,
  for example anharmonicity strongly enhanced at the interface, or a direct coupling of the
  metal's electrons to the semiconductor's phonons. The latter cannot act at AlN/GaN, which has no
  free electrons, so the two readings make different predictions.

## 2. Closest prior art (read for this amendment)

- Maassen and co-workers, arXiv:1808.10713: a McKelvey-Shockley (two-flux) Boltzmann model of a
  finite Si/Ge junction with first-principles three-phonon rates. Inelastic scattering is in the bulk
  only, the interface is elastic, and it is studied at 300 K only. It is the closest relative of our
  stage 2a solver `tbc_bte` and is now cited there. It does not state the scale invariance, and it
  defines the jump locally at the interface rather than against extrapolated bulk profiles.
- Le et al., PRB 95, 245417 (2017) (MD, Ar-like solids): interfacial anharmonicity contributes "a
  surprisingly small part" of G, and bulk anharmonicity most of it. Their domains are finite. By
  our scaling argument, a T-dependence from bulk thermalization needs a finite length, quantum
  heat capacities or T-dependent harmonic properties.
- Guo et al., arXiv:2103.10087 (anharmonic NEGF): the net anharmonic effect at the interface is
  "moderate".
- Feng, Zhong, Shi and Ruan, PRB 99, 045301 (2019) (modal NEMD, Si/Ge): inelastic transport over
  50 % of G at room temperature through a "bridging" effect.
- Murakami et al., arXiv:1506.02377: inelastic transmission is large right at the interface but is
  capped by transition regions where phonons relax back to bulk.
- Hopkins, J. Appl. Phys. 106, 013528 (2009); Hopkins, Duda and Norris, J. Heat Transfer (2011,
  AIM); Duda, Norris and Hopkins, J. Heat Transfer 133, 074501 (2011) (linear rise in the classical
  limit). These are DOS- and occupation-based. We could not read the full texts (paywalled PDFs),
  so we do not claim our formula differs from theirs in every detail. It differs in using computed
  three-phonon rates for the strength. One Hopkins-group classical-limit model, as summarized by
  others, makes the temperature dependence of G come from the materials' bulk scattering rather
  than their heat capacities. That is in direct tension with our stage 2a scaling argument for
  semi-infinite solids next to an elastic interface. We will read it in full and resolve the
  tension when stage 2a is graded.

## 3. The model (tbc_inel.py)

On side i, the atoms of the terminating layer have a fraction f_i of their nearest-neighbour bonds
across the interface: fcc(111) 3 of 12, diamond(111) single-bond termination 1 of 4, wurtzite(0001)
1 of 4. The model assumes ideal, bulk-terminated, unreconstructed epitaxy. A side-i phonon of
frequency nu in that layer scatters across at its bulk three-phonon rate Gamma_i(nu, T) (phono3py,
stage 2a) times Q_ij(nu). Q_ij(nu) is the share of its classically weighted two-phonon phase space
in which its energy ends on side j, when each partner on a cross bond sits on either side with
that side's per-atom density of states. Q = 1/2 for two identical sides (unit test), and Q -> 0
when the other side has no states where the energy must go (unit test). Then

    G_inel(T) = 1/2 * sum_i t_i * sum_nu C_i(nu, T) Gamma_i(nu, T) Q_ij(nu),   t_i = f_i sigma_i / n_i,

where sigma_i is the terminating layer's areal atom density and n_i the bulk atom density
(t = 0.58, 0.39 and 0.32 Angstrom for Al(111), Si(111) and GaN(0001)). The 1/2 removes the double
count of each cross event. In the classical limit G_inel is exactly proportional to T (unit test).
The channel adds in parallel to stage 2a's elastic two-stream conductance G_BTE. For Al interfaces,
the electron-phonon term (same inputs as stages 1 and 2a, 1.21 GW m^-2 K^-1) is in series:

    G_model(T) = [1/G_ep + 1/(G_BTE(T) + s G_inel(T))]^-1,  with s = 1 for the prediction.

Inputs: stage 2a's npz files (orb-v3 primary, MACE-MPA-0 cross-check), published to
results/anh/ before this job runs. Frequency grid 0.1 THz.

Known simplifications: the bulk rate stands for the cross-bond anharmonicity (no interfacial
force constants); per-atom DOS is atom-averaged (in GaN the terminating Ga layer is weighted like
N); RTA rates; three-phonon processes only; only the first layer on each side; and the 1/2
double-count factor is a convention. Each of these changes G_inel by a constant factor of order
one, which is why the grades use a scale factor s* with wide bands.

## 4. Expectations (graded from tbi-summary.json / the `TBI:` line)

Measured inputs: G at 400 and 700 K for Al/Si Sample 1 and Al/GaN, read off Li et al. Fig. 1 by a
reader who has not seen any model output. They are recorded in data/meas-li2022.json before the job
is filed. If the measurements stop below 700 K, a dated amendment moves the upper temperature to the
highest measured one before the job is filed.

- **I0** (pipeline): with classical heat capacities, G_inel(700)/G_inel(400) lies in [1.55, 1.80]
  for Al/Si and Al/GaN. (Exactly 1.75 if the rates were proportional to T; phono3py's Bose factors
  lower it slightly.) Expected PASS (90 %).
- **I1** (main: magnitude): solve for the scale factor s* that makes G_model(700) - G_model(400)
  equal the measured rise. PASS if s* is in [0.5, 2] for both Al/Si and Al/GaN. Expected PASS (20 %).
  Our rough estimate before any data: G_inel is about 10-15 MW m^-2 K^-1 at 300 K, mostly from the Al
  side, which is likely too small for the visible rise. **Reading** recorded with I1: "too-weak" if
  s* > 8 for both pairs (our most likely outcome, 45 %), "right-order", "intermediate" or
  "too-strong" (s* < 0.25 for both).
- **I2** (transferability with one universal constant): s* fixed on Al/Si predicts the Al/GaN rise
  within +-35 %. Expected PASS (40 %). This is graded even if I1 fails. A PASS with s* far from 1 would
  mean the shape is right and only the strength is off, so one constant of the framework (not of a
  specimen) would be needed.
- **I3**: G_model(300 K) with s = 1 within 20 % of 379 (Al/Si) and 423 (Al/GaN). Expected PASS (35 %).
- **I4**: AlN/GaN G_model(300 K) with s = 1 within [256, 588] MW m^-2 K^-1 (0.8 x 320 to 1.2 x 490,
  the two experiments). Expected PASS (70 %).
- Recorded, not graded: the unweighted (JDOS) kernel, the classical-heat-capacity variant, each
  side's share, the heat-capacity-weighted Q and rates, the MACE-MPA-0 cross-check, and s* for each
  potential.

If I1 reads "too-weak", stage 3 does not tune this channel up. It tests the electron route
(Al/GaN against AlN/GaN, and the Au/Ti and Au/Ni benchmarks) before anything else.
