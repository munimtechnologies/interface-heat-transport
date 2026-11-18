# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""tbc_phonons: bulk phonons (mesh frequencies and group velocities) of the benchmark materials of
Grand Challenge 05 from open universal machine-learned interatomic potentials, via phonopy.

Potentials (both licensed for commercial use):
- "mace-mpa-0": MACE-MPA-0 medium (MIT), ACEsuit/mace-foundations release mace_mpa_0.
- "orb-v3": orb-v3-conservative-inf-omat (Apache-2.0), orbital-materials/orb-models.

No measured interface conductance enters anything here: only bulk crystal structures.
"""
import time

import numpy as np

VERSION = "mv2 tbc-phonons 0.1.0"

MACE_MPA0_URL = ("https://github.com/ACEsuit/mace-foundations/releases/download/mace_mpa_0/"
                 "mace-mpa-0-medium.model")

# name: (builder args, supercell of the conventional/hexagonal cell, mesh)
# Lattice parameters are room-temperature experimental values (starting points; relaxed below).
MATERIALS = {
    "Al": dict(kind="fcc", el="Al", a=4.0495, sc=(3, 3, 3), mesh=32),
    "Cu": dict(kind="fcc", el="Cu", a=3.6149, sc=(3, 3, 3), mesh=32),
    "Au": dict(kind="fcc", el="Au", a=4.0782, sc=(3, 3, 3), mesh=32),
    "Ni": dict(kind="fcc", el="Ni", a=3.5240, sc=(3, 3, 3), mesh=32),
    "Si": dict(kind="diamond", el="Si", a=5.4310, sc=(3, 3, 3), mesh=28),
    "C": dict(kind="diamond", el="C", a=3.5668, sc=(3, 3, 3), mesh=28),
    "TiC": dict(kind="rocksalt", el="TiC", a=4.3280, sc=(3, 3, 3), mesh=28),
    "Ti": dict(kind="hcp", el="Ti", a=2.9508, c=4.6855, sc=(4, 4, 3), mesh=24),
    "GaN": dict(kind="wurtzite", el="GaN", a=3.1890, c=5.1850, u=0.377, sc=(4, 4, 3), mesh=24),
    "AlN": dict(kind="wurtzite", el="AlN", a=3.1120, c=4.9820, u=0.382, sc=(4, 4, 3), mesh=24),
}

# Highest phonon frequency (THz) for the sanity check. Polar crystals: highest TO mode, because no
# non-analytic (LO-TO) correction is applied. Values from Raman / neutron data (rounded):
# Al 9.7 (neutron), Cu 7.3, Au 4.6, Ni 8.9, Si 15.6 (520 cm^-1 Raman), C 39.9 (1332 cm^-1),
# GaN E1(TO) 559 cm^-1 = 16.8, AlN E1(TO) 671 cm^-1 = 20.1.
REF_MAX_THZ = {"Al": 9.7, "Cu": 7.3, "Au": 4.6, "Ni": 8.9, "Si": 15.6, "C": 39.9, "GaN": 16.8,
               "AlN": 20.1}


def build(name):
    from ase.build import bulk
    m = MATERIALS[name]
    if m["kind"] in ("fcc", "diamond", "rocksalt"):
        return bulk(m["el"], m["kind"], a=m["a"], cubic=True)
    if m["kind"] == "hcp":
        return bulk(m["el"], "hcp", a=m["a"], c=m["c"])
    if m["kind"] == "wurtzite":
        return bulk(m["el"], "wurtzite", a=m["a"], c=m["c"], u=m["u"])
    raise ValueError(name)


def load_calculator(model):
    """Return (ase calculator, displacement in Angstrom)."""
    if model == "mace-mpa-0":
        from mace.calculators import mace_mp
        last = None
        # MPA-0 medium by URL, then by name; never silently fall back to another model.
        for spec in (MACE_MPA0_URL, "medium-mpa-0"):
            try:
                calc = mace_mp(model=spec, default_dtype="float64", device="cpu")
                print(f"MACE loaded from {spec}", flush=True)
                return calc, 0.01
            except Exception as e:  # try the next spelling
                last = e
        raise RuntimeError(f"MACE-MPA-0 could not be loaded: {last}")
    if model == "orb-v3":
        from orb_models.forcefield import pretrained
        try:
            from orb_models.forcefield.inference.calculator import ORBCalculator
        except ImportError:  # older orb-models layout
            from orb_models.forcefield.calculator import ORBCalculator
        out = pretrained.orb_v3_conservative_inf_omat(device="cpu", precision="float32-highest")
        if isinstance(out, tuple):
            orbff, adapter = out
            return ORBCalculator(orbff, atoms_adapter=adapter, device="cpu"), 0.02
        return ORBCalculator(out, device="cpu"), 0.02
    raise ValueError(model)


def relax(atoms, calc, fmax=2e-3, steps=400):
    """Relax cell and positions keeping the space group."""
    from ase.constraints import FixSymmetry
    from ase.optimize import BFGS
    try:
        from ase.filters import FrechetCellFilter as CellFilter
    except ImportError:  # older ase
        from ase.constraints import ExpCellFilter as CellFilter
    atoms = atoms.copy()
    atoms.calc = calc
    atoms.set_constraint(FixSymmetry(atoms))
    opt = BFGS(CellFilter(atoms), logfile=None)
    converged = bool(opt.run(fmax=fmax, steps=steps))
    atoms.set_constraint()
    return atoms, converged


def phonon_mesh(atoms, calc, sc, mesh, disp):
    """Force constants by finite displacements, then frequencies (THz), group velocities (m/s),
    weights and primitive-cell volume (m^3) on a full (unreduced) mesh."""
    from ase import Atoms
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    unit = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=np.array(atoms.cell),
                        scaled_positions=atoms.get_scaled_positions())
    ph = Phonopy(unit, supercell_matrix=np.diag(sc), primitive_matrix="auto")
    ph.generate_displacements(distance=disp)
    forces = []
    for s in ph.supercells_with_displacements:
        a = Atoms(s.symbols, cell=s.cell, scaled_positions=s.scaled_positions, pbc=True)
        a.calc = calc
        f = np.array(a.get_forces(), dtype=float)
        forces.append(f - f.mean(axis=0))
    ph.forces = np.array(forces)
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    ph.run_mesh([mesh] * 3, is_time_reversal=False, is_mesh_symmetry=False,
                with_group_velocities=True)
    try:
        m = ph.mesh
        freq, gv, w = m.frequencies, m.group_velocities, m.weights
    except AttributeError:  # older phonopy
        d = ph.get_mesh_dict()
        freq, gv, w = d["frequencies"], d["group_velocities"], d["weights"]
    freq = np.array(freq, dtype=float)
    gv = np.array(gv, dtype=float) * 100.0  # THz Angstrom -> m/s
    w = np.array(w, dtype=float)
    vol = float(ph.primitive.volume) * 1e-30
    return dict(freq=freq, vel=gv, weights=w, vcell=vol, n_disp=len(forces),
                n_prim=len(ph.primitive.symbols))


def material_phonons(name, calc, disp):
    t0 = time.time()
    m = MATERIALS[name]
    atoms, conv = relax(build(name), calc)
    res = phonon_mesh(atoms, calc, m["sc"], m["mesh"], disp)
    f = res["freq"]
    res.update(
        name=name, relax_converged=conv,
        cell_lengths=[float(x) for x in atoms.cell.lengths()],
        n_atoms_per_m3=float(res["n_prim"] / res["vcell"]),
        max_thz=float(f.max()), min_thz=float(f.min()),
        imag_fraction=float(np.mean(f < -0.05)),
        seconds=round(time.time() - t0, 1),
    )
    return res
