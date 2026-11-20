# PLAN amendment 0.5.12: the cover band rule is built (tbc_cover.py 0.1.0); one case the rule of 0.5.11 did not foresee (a side with no states above its cutoff) is specified before any run, and the literal reading is kept as information

MV2 of Munim, Inc. (mv2@munimtech.com). 2026-11-20. Made with AI agents working in MV2's name; no
thermal-transport specialist has checked it yet. Code: Apache-2.0. Text: CC BY 4.0.

Written **before job `tbc-cover-v1` has run** and before any number from the cover rule exists.
Tests V0-V5 and their thresholds (amendment 0.5.11 sec. 4) are unchanged. `tbc-anh-v2` is still
queued; no published module was changed (`tbc_cover.py` is new and imports tbc_bte, tbc_anh and
tbc_fill unchanged).

## 1. What was built

`tbc_cover.py` 0.1.0 with `test_tbc_cover.py` (12 tests, numbers worked by hand) and
`fetch_npz_cover.py`. Job `tbc-cover-v1` (python:3.12, 2 CPUs, 2 GB) reloads the four orb-v3 npz
files of tbc-tin-v1 and, for the earlier pairs, the published orb-v3 npz files of
`results/anh/v1/` (Al, Si, GaN, AlN) and `results/nonmetal/zno-v1/` (ZnO, GaN). Each TiN row is
solved five ways: base (job lengths), merge (tbc-fill-v1), **cover** (graded), **cover with bands of
at least two fine bins** (V4) and the **literal** reading below (information). Last line `COV:`.

## 2. The case 0.5.11 did not foresee, specified now

0.5.11 sec. 4: "consecutive fine bins (0.2 THz) are grouped into one band until the group contains
at least one bin sampled on the lifetime mesh on both sides ... A tail left without a sampled bin on
both sides joins the last band."

While writing the code we saw that, read literally, this folds one material's whole spectrum above
the other material's cutoff into a single band. Ge's phonons end near 9 THz and Si's near 16, but
TiN's optical branch lies above them; diamond's spectrum runs far above TiN's. Under the literal
rule no band above the lower cutoff can ever close, so all of those modes join the last band shared
with the other side, and the DMM then lets TiN's optical phonons transmit into Ge's top acoustic
band. That breaks the elastic (frequency-conserving) premise of the whole model. It is not what the
rule was for: the rule exists to remove *filled* bins (K > 0 but no lifetime-mesh mode), not bins
where a side simply has no states.

**Specification (graded):** a band closes when, on EACH side, it holds a sampled bin (K > 0 and
kappa > 0) **or holds no state on that side** (K = 0 in every one of its bins). A side's bins above
its own last sampled bin are not used to close bands; they join, on that side only, the band that
holds its last sampled bin (the pre-registered "tail joins the last band", applied per side). K and
kappa are summed per band and side, so each side's total K and total kappa are conserved, and no
band with K > 0 is left without kappa. The first band is the merge band of tbc-fill-v1 whenever both
sides have states from the lowest bin, which V0 checks.

**Literal reading (reported, never graded):** close a band only when both sides hold a sampled
bin; any tail joins the last band. Its G and DMM are printed next to the graded ones, so the effect
of this specification is visible.

This specification was made after reading the rule and the code, and before any cover number
existed. We did not try other variants.

## 3. Details fixed now

- Earlier pairs (V5 and information) use the oriented geometry of the stage that solved them:
  `M_or` with the c-axis conductivity (`v2_z`) for wurtzite AlN, GaN and ZnO, and the isotropic one
  for cubic Al and Si; orb-v3 for every material. A pair whose files are missing or unreadable is
  reported with its error, and V5 is then UNDETERMINED.
- V1 compares each side's model conductivity (sum of 2 K l over bands) with the sum of its fine-bin
  phono3py conductivity. Conductivity in fine bins where that side has no K at all is printed as
  "lost"; it counts against V1.
- V0 compares the first cover band's upper edge with tbc-fill-v1's merged range (4.8, 1.0, 0.6 THz)
  and base/merge with fill-summary.json within max(0.5, 1 %).
- Printed per row: band count, narrowest and widest band, every band wider than one fine bin, each
  side's last sampled frequency, Simons and linear-bath values on the cover bands, and cover plus the
  tbc-tin-v1 inelastic term against Khan et al. (2024), information only.

## 4. Priors (unchanged from 0.5.11)

V0 85 %, V1 90 %, V2 55 %, V3 70 %, V4 60 %, V5 50 %. Reading rules of 0.5.11 sec. 4 apply as
written.

## 5. Limits

As in 0.5.11 sec. 5. The specification in section 2 is a judgement made by AI agents; the literal
numbers are printed so anyone can see what it changes.
