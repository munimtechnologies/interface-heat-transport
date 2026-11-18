# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""tbc_anh: bulk three-phonon (and isotope) lifetimes of the benchmark materials of Grand Challenge
05 from open universal machine-learned potentials, via phono3py (relaxation-time approximation),
and the per-frequency-bin quantities the interface model tbc_bte needs:
- kappa_b(T): each frequency bin's share of the bulk thermal conductivity along a direction;
- K_b(T): each bin's phonon radiation conductance toward the interface (from a dense phonopy mesh
  of the same second-order force constants).
No measured interface conductance enters anything here.
"""
import time

import numpy as np

import tbc_phonons as tp
import tbc_transport as tt

VERSION = "mv2 tbc-anh 0.1.0"

# third-order supercell (multiples of the conventional / hexagonal cell built by tbc_phonons), pair
# cutoff for third-order displacements (Angstrom, None = all pairs) and the phono3py mesh
# (of the primitive cell). Second-order supercells are those of tbc_phonons.MATERIALS.
ANH = {
    "Al": dict(sc3=(3, 3, 3), cut=5.0, mesh=(18, 18, 18), lattice="cubic"),
    "Si": dict(sc3=(2, 2, 2), cut=None, mesh=(16, 16, 16), lattice="cubic"),
    "C": dict(sc3=(2, 2, 2), cut=None, mesh=(16, 16, 16), lattice="cubic"),
    "GaN": dict(sc3=(3, 3, 2), cut=4.0, mesh=(12, 12, 8), lattice="hex"),
    "AlN": dict(sc3=(3, 3, 2), cut=4.0, mesh=(12, 12, 8), lattice="hex"),
}
# Bulk thermal conductivity at 300 K (W/m/K) for the A1 sanity ranges (PLAN-amendment-0.2):
# natural-isotope experiments, Si ~150, diamond ~2200, GaN ~230, AlN ~300.
KAPPA_RANGE_300 = {"Si": (100.0, 200.0), "C": (1000.0, 3500.0), "GaN": (120.0, 400.0),
                   "AlN": (150.0, 450.0)}
DISP = 0.03  # Angstrom


def _forces(cells, calc):
    from ase import Atoms
    out = []
    for s in cells:
        if s is None:
            out.append(None)
            continue
        a = Atoms(s.symbols, cell=s.cell, scaled_positions=s.scaled_positions, pbc=True)
        a.calc = calc
        f = np.array(a.get_forces(), dtype=float)
        out.append(f - f.mean(axis=0))
    n = next(len(x) for x in out if x is not None)
    return np.array([x if x is not None else np.zeros((n, 3)) for x in out])


def _get(obj, name):
    v = getattr(obj, name, None)
    if v is None and hasattr(obj, "get_" + name):
        v = getattr(obj, "get_" + name)()
    return None if v is None else np.array(v)


def anharmonic(name, calc, temps, log=print):
    """Relax, compute fc2 and fc3, run RTA conductivity with isotopes. Returns a dict with the
    irreducible-mesh data (frequencies THz, |v|^2/3 and v_z^2 in m^2/s^2, weights, gamma_3ph and
    gamma_iso in THz), phono3py's kappa tensor, the relaxed cell and the full fc2."""
    from phono3py import Phono3py
    from phonopy.structure.atoms import PhonopyAtoms
    t0 = time.time()
    m = tp.MATERIALS[name]
    s = ANH[name]
    atoms, conv = tp.relax(tp.build(name), calc)
    unit = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=np.array(atoms.cell),
                        scaled_positions=atoms.get_scaled_positions())
    ph3 = Phono3py(unit, supercell_matrix=np.diag(s["sc3"]), primitive_matrix="auto",
                   phonon_supercell_matrix=np.diag(m["sc"]), log_level=0)
    ph3.generate_displacements(distance=DISP, cutoff_pair_distance=s["cut"])
    if not ph3.phonon_supercells_with_displacements:
        ph3.generate_fc2_displacements(distance=DISP)
    n3 = sum(1 for c in ph3.supercells_with_displacements if c is not None)
    log(f"  {name}: relaxed (conv={conv}) a,b,c={['%.4f' % x for x in atoms.cell.lengths()]}; "
        f"fc3 supercells {n3}, fc2 supercells {len(ph3.phonon_supercells_with_displacements)}")
    ph3.forces = _forces(ph3.supercells_with_displacements, calc)
    ph3.phonon_forces = _forces(ph3.phonon_supercells_with_displacements, calc)
    t_forces = time.time() - t0
    ph3.produce_fc3()
    ph3.produce_fc2(is_compact_fc=False)
    ph3.mesh_numbers = list(s["mesh"])
    ph3.init_phph_interaction()
    ph3.run_thermal_conductivity(temperatures=list(temps), is_isotope=True, write_kappa=False)
    tc = ph3.thermal_conductivity
    freq = _get(tc, "frequencies")
    gv = _get(tc, "group_velocities") * 100.0          # THz Angstrom -> m/s
    w = _get(tc, "grid_weights").astype(float)
    gam = _get(tc, "gamma")[0]                         # (num_temp, num_gp, num_band)
    giso = _get(tc, "gamma_isotope")
    giso = np.zeros_like(gam[0]) if giso is None else giso[0]
    kappa = _get(tc, "kappa")[0]                       # (num_temp, 6)
    vol = float(ph3.primitive.volume) * 1e-30
    return dict(name=name, relax_converged=bool(conv),
                cell_lengths=[float(x) for x in atoms.cell.lengths()],
                freq=freq, v2_iso=np.sum(gv * gv, axis=-1) / 3.0, v2_z=gv[..., 2] ** 2,
                weights=w, gamma3=gam, gamma_iso=giso, kappa_ph3=kappa, vcell=vol,
                n_fc3=n3, unit=unit, sc2=m["sc"], fc2=np.array(ph3.fc2),
                lattice=s["lattice"], seconds_forces=round(t_forces, 1),
                seconds=round(time.time() - t0, 1))


def kappa_bins(r, T_index, T, edges_hz, direction, with_isotope=True, classical=False):
    """Per-bin conductivity (W/m/K) along 'direction' ('iso' = one third of the trace, which is
    exact for cubic crystals along any axis; 'z' = along the hexagonal c axis)."""
    f = r["freq"] * 1e12
    g = r["gamma3"][T_index] + (r["gamma_iso"] if with_isotope else 0.0)
    v2 = r["v2_iso"] if direction == "iso" else r["v2_z"]
    ok = (f > 1e9) & (g > 0.0)
    tau = np.where(ok, 1.0 / (4.0 * np.pi * np.where(ok, g, 1.0) * 1e12), 0.0)
    c = np.full(f.shape, tt.KB) if classical else tt.mode_heat_capacity(np.clip(f, 0.0, None), T)
    contrib = (r["weights"][:, None] * c * v2 * tau) / (r["weights"].sum() * r["vcell"])
    return np.histogram(f[ok], bins=edges_hz, weights=contrib[ok])[0]


def kappa_ph3_dir(r, T_index, direction):
    k = r["kappa_ph3"][T_index]
    return float(np.mean(k[:3])) if direction == "iso" else float(k[2])


def flux_mesh(r, mesh=None):
    """Dense full-mesh frequencies and velocities from the fc2 of 'anharmonic' (phonopy)."""
    from phonopy import Phonopy
    mesh = mesh or tp.MATERIALS[r["name"]]["mesh"]
    ph = Phonopy(r["unit"], supercell_matrix=np.diag(r["sc2"]), primitive_matrix="auto")
    ph.force_constants = r["fc2"]
    ph.run_mesh([mesh] * 3, is_time_reversal=False, is_mesh_symmetry=False,
                with_group_velocities=True)
    mm = ph.mesh
    return dict(freq=np.array(mm.frequencies, dtype=float),
                vel=np.array(mm.group_velocities, dtype=float) * 100.0,
                weights=np.array(mm.weights, dtype=float),
                vcell=float(ph.primitive.volume) * 1e-30)


def K_bins(M_fine, edges_fine, edges_coarse, T, classical=False):
    """Bin radiation conductances K_b (W/m^2/K) from a fine spectral flux M(nu)."""
    import tbc_bte as tb
    c, w = tt.centers_widths(edges_fine)
    cap = np.full(c.shape, tt.KB) if classical else tt.mode_heat_capacity(c, T)
    return tb.rebin(cap * M_fine * w, edges_fine, edges_coarse)
