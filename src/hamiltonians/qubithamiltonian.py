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


if __name__ == "__main__":
    fermion_op = build_n2_fermionic_operator()
    qubit_op = map_to_qubit_operator(fermion_op)
