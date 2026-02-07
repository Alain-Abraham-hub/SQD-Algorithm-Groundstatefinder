"""Initial parameters for quantum ansatze from classical CCSD calculations."""

from __future__ import annotations

from typing import Tuple

import numpy as np

try:
    from pyscf import cc, gto, scf
except ImportError as exc:
    raise ImportError(
        "PySCF is required for CCSD calculations. "
        "Install it with `pip install pyscf`."
    ) from exc


def get_ccsd_amplitudes(
    *,
    bond_length: float = 1.0977,
    basis: str = "6-31g",
    charge: int = 0,
    spin: int = 0,
    unit: str = "Angstrom",
    symmetry: bool | str = "Dooh",
) -> Tuple[np.ndarray, np.ndarray]:
    """Calculate CCSD t1 and t2 amplitudes for N2 molecule.

    The t1 amplitudes correspond to single excitation operators and t2 amplitudes
    correspond to double excitation operators. These can be used as initial
    parameters for quantum ansatze like UCCSD.

    Args:
        bond_length: N–N bond length.
        basis: Basis set (default 6-31G).
        charge: Total molecular charge.
        spin: 2 * S (0 for singlet).
        unit: Distance unit for geometry.
        symmetry: Symmetry group (default "Dooh").

    Returns:
        Tuple of (t1, t2) amplitudes where:
            - t1: single excitation amplitudes, shape (n_occ, n_virt)
            - t2: double excitation amplitudes, shape (n_occ, n_occ, n_virt, n_virt)
    """
    # Build molecule geometry
    half = bond_length / 2.0
    geometry = f"N 0.0 0.0 {-half}; N 0.0 0.0 {half}"

    # Create molecule
    mol = gto.Mole()
    mol.atom = geometry
    mol.basis = basis
    mol.charge = charge
    mol.spin = spin
    mol.unit = unit
    mol.symmetry = symmetry
    mol.build()

    # Run Hartree-Fock
    if spin == 0:
        mf = scf.RHF(mol)
    else:
        mf = scf.ROHF(mol)
    mf.kernel()

    # Run CCSD
    mycc = cc.CCSD(mf)
    mycc.kernel()

    # Extract amplitudes
    t1 = mycc.t1
    t2 = mycc.t2

    return t1, t2


def flatten_ccsd_amplitudes(
    t1: np.ndarray, t2: np.ndarray
) -> np.ndarray:
    """Flatten CCSD amplitudes into a 1D parameter array for ansatz initialization.

    Args:
        t1: Single excitation amplitudes, shape (n_occ, n_virt).
        t2: Double excitation amplitudes, shape (n_occ, n_occ, n_virt, n_virt).

    Returns:
        1D numpy array containing all amplitudes flattened.
    """
    return np.concatenate([t1.flatten(), t2.flatten()])


if __name__ == "__main__":
    # Example: Get CCSD amplitudes for N2
    t1, t2 = get_ccsd_amplitudes()
    
    print(f"t1 shape: {t1.shape}")
    print(f"t2 shape: {t2.shape}")
    print(f"\nt1 amplitudes (single excitations):\n{t1}")
    print(f"\nt2 amplitudes (double excitations, showing first 2x2x2x2 block):\n{t2[:2, :2, :2, :2]}")
    
    # Flatten for use in quantum circuit
    params = flatten_ccsd_amplitudes(t1, t2)
    print(f"\nTotal number of parameters: {len(params)}")
