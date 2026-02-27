"""Qubit mapping transformations for fermionic Hamiltonians."""

from __future__ import annotations

try:
    from openfermion.transforms import jordan_wigner
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "OpenFermion is required for Jordan-Wigner transformation. "
        "Install it with `pip install openfermion`."
    ) from exc


def map_to_qubit_operator(fermion_operator: "FermionOperator") -> "QubitOperator":
    """Map fermionic operator to qubit operator using Jordan-Wigner.

    Args:
        fermion_operator: FermionOperator from fermionic_hamiltonian.py.

    Returns:
        QubitOperator representing the transformed Hamiltonian.
    """
    return jordan_wigner(fermion_operator)
