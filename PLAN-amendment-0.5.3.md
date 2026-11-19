# Amendment 0.5.3: what Khan et al.'s Supplementary Information says, and one interface, two temperature laws

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-19. Made with AI agents working in MV2's name; no
human expert has checked it. Text CC BY 4.0. Data: `data/read-khansi.json` (a blind reading of
the papers by an AI reading agent; not checked by a person).

## 1. Khan et al. (APL 125, 041601, 2024), Supplementary Information

Stage 3b (amendment 0.5.2, job `tbc-tin-v1`, queued) left two points open. The SI settles one
and leaves the other open:

- **SiC is 4H** (hexagonal). Stage 3b already leaves TiN/SiC out; a cubic stand-in would have been
  the wrong crystal. A hexagonal SiC kind is a possible stage 3c, only after 3b returns.
- **Substrate orientations and the epitaxial relation are not stated** in the main text or SI.
  The (001) geometry of stage 3b therefore stays an assumption, and the 3b grade will say so.
  `tbc-tin-v1` is not changed.
- TiN is 38-46 nm thick, reactively sputtered; all substrates were HF-dipped to remove oxide.
- Khan et al.'s own elastic limit and DMM use DFPT dispersions over the full zone for TiN, diamond
  and Si (not Ge or SiC), with the transmission rules quoted in the data file.
- **Electron-phonon coupling**: they estimate a TiN-side channel of about 0.1 GW m^-2 K^-1 (about a
  tenth of the TiN/Si conductance) and call it unlikely to matter. This matches the size of the
  sensitivity term in amendment 0.5.2; it stays a sensitivity, not part of any prediction.
- Their error bars come from fitting-contour analysis and grow with temperature (Si 18 % -> 37 %,
  Ge 12 % -> 28 % from 300 to 770 K). The 3b grade (T2, T3) will be reported with these widths.

## 2. Al on sapphire: one interface, two temperature laws

Two measurements of the same nominal interface, from the readings in this repository:

| sample | preparation | G at ~300 K | G at highest T | ratio | source |
|---|---|---|---|---|---|
| Hopkins et al. 2007 (IJT 28, 947) | factory-polished sapphire, solvent clean, 400 K bake, 100 W Ar backsputter, Al evaporated | 185 (295 K) | 341 (475 K) | **1.84 +- 0.08** | data/read-khansi.json |
| Koh et al. 2020 (PRB 102, 205304) | epitaxial Al on sapphire | 315 (300 K) | 325 (498 K) | **1.03 +- 0.01** | data/read-hight.json |

(Ratio errors: +-3 % reading uncertainty on each Hopkins value, combined in quadrature; +-1 MW
m^-2 K^-1 on each Koh value.)

The two samples differ by a factor 1.7 at room temperature and by a factor 1.8 in how much they
rise over the same ~180 K. Any model that takes only bulk properties of Al and Al2O3 (dispersions,
lifetimes, electron-phonon coupling) predicts the **same** G(T) for both. So, for a no-fit model,
the temperature law of a metal/non-metal interface is **not** a property of the two materials
alone: preparation (here, backsputter damage versus epitaxy) changes not only the level but the
slope. This is not a new observation in the field. Hopkins et al. attributed their rise to
inelastic scattering, and Koh et al. reproduced their flat curve with an elastic model. But it
sets a rule for how this project grades.

**Rule, fixed now.** (a) Blind tests of the model use only interfaces described as epitaxial or
atomically abrupt (Koh 2020, Khan 2024 as epitaxial-like sputtered TiN, Li 2022 Sample 1,
Gaskins 2018 ZnO/GaN). (b) Damaged or backsputtered interfaces (Hopkins 2007, Li 2022 Sample 2)
are reported separately, as "outside the model's scope", and never counted as a pass or a fail.
(c) Any claim to the Discovered Materials challenge will state this scope in its first paragraph.

Rule (b) does not change amendment 0.4's grade. That grade used Sample 2 only to exclude a
bulk-only electron channel: a channel fixed by bulk properties alone must produce the same rise
on both samples, so the rule's reasoning applies to it.
