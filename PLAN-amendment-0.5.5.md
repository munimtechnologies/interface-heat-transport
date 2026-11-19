# PLAN amendment 0.5.5: first stage-2a returns, two lifetime grades, and one new pre-registered check

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-19. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text and data: CC BY 4.0.
Grades of PLAN-amendment-0.2.md (stage 2a, bulk anharmonicity next to an elastic interface) that the
partial output of job `tbc-anh-v1` allows, written before any interface table of stage 2a was seen.
Nothing in earlier amendments changes.

## 1. What came back

`tbc-anh-v1` ran 65 minutes. It computed three-phonon and isotope lifetimes (phono3py, RTA) with
**orb-v3** (orb-v3-conservative-inf-omat, Apache-2.0) for all five materials, then started the
**MACE-MPA-0** cross-check and finished Al and Si. It was killed at its 3 GB memory cap during
MACE GaN (466 third-order supercells in float64), before it built the interface tables or wrote
its grades. Each material's lifetimes were saved as it finished, so nothing computed is lost:
[results/anh/v1/](results/anh/v1/) holds the seven saved files, the job's progress log, its tests
(all passed), package versions and file hashes.

Bulk thermal conductivity at 300 K from phono3py (W m^-1 K^-1; natural isotopes included, three-phonon
only in brackets; for GaN and AlN the first value is one third of the trace, the second along c):

| material | orb-v3 | MACE-MPA-0 | A1 range (from experiment) |
|---|---|---|---|
| Al (phonons only) | 1.5 (1.5) | 15.6 (15.6) | not graded |
| Si | 61.5 (66.1) | 56.8 (59.3) | 100-200 |
| GaN | 141.1 (220.1); c: 144.7 (212.2) | not reached | 120-400 |
| AlN | 169.3 (169.5); c: 171.7 (171.9) | not reached | 150-450 |
| diamond | 888.4 (1053.1) | not reached | 1000-3500 |

## 2. Grades

- **A0 (pipeline), orb-v3: PASS.** For every material and direction the sum of our per-bin
  conductivities equals phono3py's own value to the 0.1 W m^-1 K^-1 printed (ratio 1.000; limit 3 %).
  MACE-MPA-0: Al and Si also equal; the grade is completed by `tbc-anh-v2`.
- **A1 (lifetime sanity), orb-v3: FAIL** (expected PASS at 60 %). Si 61.5 is below 100 and diamond
  888 below 1000; GaN 141 and AlN 169 are inside their ranges.
- **A1, MACE-MPA-0: FAIL.** Si 56.8 is below 100; the materials still to come cannot change that.

Both open universal potentials, in our setup, give silicon about 40 % of its measured conductivity
(natural silicon is about 150), and orb-v3 gives diamond about 40 % of its measured value (about
2200). This job cannot say whether the potentials or our settings are the cause (0.03 Angstrom
displacements, orb-v3 in float32, 2x2x2 conventional supercells for Si and diamond); both models
agree on Si within 8 %. The two potentials disagree tenfold on the phonon-only conductivity of
aluminium (1.5 against 15.6).

## 3. Why A1 failing does not settle stage 2a, and a check of that (Q1, pre-registered)

In the multiband two-stream model (tbc_bte.py) each half-space touches the other only at the
interface. Multiplying every length on **one** side by a constant lambda and stretching that side's
coordinate by the same lambda maps its solution onto itself: its bulk conductivity scales with
lambda, so its bulk temperature profile, extrapolated to the interface, does not move, and the stream
offsets l Q / kappa do not change. So G should not depend on the overall size of either side's
lifetimes at all, only on their **spectral shape** within each material. Amendment 0.2 used the
common-factor version of this to argue that bulk three-phonon scattering cannot make G rise with T in
the classical limit; the one-side version is stronger.

Consequence: a potential that gets every lifetime of a material wrong by the same factor still gives
the model's G exactly. The A1 failure therefore matters for S2, S3 and S5 only through errors in the
shape of the lifetime spectrum, which A1 does not measure. S6 (orb-v3 against MACE-MPA-0) is the
grade that sees shape differences, and the tenfold Al disagreement makes it a real test.

- **Q1** ([q1_scaling_check.py](q1_scaling_check.py), written before it was ever run): multiplying
  side 1's lengths by 0.2 or 7, or side 2's by 0.1 or 5, changes G by less than 1 %; as a control,
  changing the relative lengths within one side changes G by more than 1 %. Expected PASS (90 %).
  It runs in job `tbc-anhpairs-v1` as information and does not stop the job; graded in 0.5.6.

## 4. Information, not graded

Before it was killed, v1 printed one interface number, its grid check: orb-v3, oriented Al(111)/Si(111),
300 K, isotopes included, **G_BTE = 648.33 MW m^-2 K^-1** with 400 grid nodes (648.35 with 200).
In series with the Al electron-phonon conductance of stage 1 (sqrt(2.45e17 x 6.0) = 1212 MW m^-2 K^-1)
that is 422 MW m^-2 K^-1 (computed by hand here; `anh_pairs.py` prints it), against 379 measured
(Li et al. 2022, sample S1): 11 % high, with no fitted parameter. One number on one interface is not a
test; S3 places it against the conventional and gray-corrected DMM.

## 5. What runs next (expectations unchanged from amendment 0.2)

- **`tbc-anhpairs-v1`** (numpy only, 1 GB): [anh_pairs.py](anh_pairs.py) reads the saved orb-v3 files,
  rebuilds the per-bin conductivities, the interface tables and the grades with the functions of
  run_tbc_anh.py and tbc_anh.py imported unchanged, and writes `tba-summary.json` and a last line
  `TBA:`. It grades **S2, S3 and S5 for orb-v3**, the primary potential named in amendment 0.2,
  and Q1. Files saved by v1 lack the cell volume; it is taken from the relaxed cell lengths v1
  printed (4 decimals, error below 1e-4), and A0 is recomputed as the check.
- **`tbc-anh-v2`** (5 GB): [run_anh_safe.py](run_anh_safe.py) computes each MACE-MPA-0 material in
  its own process (a kill costs one material; a killed MACE material is retried once in float32,
  and that is recorded), then grades both potentials together: **S6**, and A0, A1, S2, S3, S5 for
  MACE-MPA-0.
- If S2 passes for orb-v3, bulk anharmonicity next to an elastic interface does not give the measured
  rise of Al/Si and Al/GaN with temperature, and the frozen model needs an inelastic interface
  channel (stage 2b, `tbc_inel.py`), as amendment 0.2 said.

## 6. Limits

The lifetimes come from universal machine-learned potentials that fail the A1 sanity ranges for Si
and diamond. Q1 is a property of our model, not of nature. The interface model is diffuse and
elastic, relaxation-time, two-stream, with no four-phonon scattering and no electron scattering of
Al phonons (amendment 0.2, sec. 3).
