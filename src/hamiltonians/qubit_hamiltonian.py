"""Qubit mapping transformations for fermionic Hamiltonians."""

from __future__ import annotations

try:
    from openfermion.ops import FermionOperator, QubitOperator
    from openfermion.transforms import jordan_wigner
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "OpenFermion is required for Jordan-Wigner transformation. "
        "Install it with `pip install openfermion`."
    ) from exc

from .fermionic_hamiltonian import build_n2_fermionic_operator


def map_to_qubit_operator(fermion_operator: FermionOperator) -> QubitOperator:
    """Map fermionic operator to qubit operator using Jordan-Wigner.

    Args:
        fermion_operator: FermionOperator from fermionic_hamiltonian.py.

    Returns:
        QubitOperator representing the transformed Hamiltonian.
    """
    return jordan_wigner(fermion_operator)


def build_qubit_hamiltonian(
    *,
    bond_length: float = 1.0977,
    basis: str = "6-31g",
    charge: int = 0,
    spin: int = 0,
    active_space: tuple[int, int] = (10, 8),
    unit: str = "Angstrom",
    symmetry: bool | str = "Dooh",
) -> QubitOperator:
    """Build the complete qubit Hamiltonian for N2.
    
    This is a convenience function combining fermionic operator generation
    and Jordan-Wigner transformation in a single call.

    Args:
        bond_length: N–N bond length.
        basis: Basis set (default 6-31G).
        charge: Total molecular charge.
        spin: 2 * S (0 for singlet).
        active_space: (n_active_electrons, n_active_orbitals).
        unit: Distance unit for geometry.
        symmetry: Symmetry group.

    Returns:
        QubitOperator representing the N2 Hamiltonian in qubit space.
    """
    fermion_op = build_n2_fermionic_operator(
        bond_length=bond_length,
        basis=basis,
        charge=charge,
        spin=spin,
        active_space=active_space,
        unit=unit,
        symmetry=symmetry,
    )
    return map_to_qubit_operator(fermion_op)


if __name__ == "__main__":
    qubit_op = build_qubit_hamiltonian()
    print(f"Qubit Hamiltonian terms: {len(qubit_op.terms)}")
