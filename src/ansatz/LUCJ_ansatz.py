"""UCJ ansatz implementation using ffsim for N2 ground state calculation."""

from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING

try:
    import ffsim
except ImportError as exc:
    raise ImportError(
        "ffsim is required for UCJ ansatz. "
        "Install it with `pip install ffsim`."
    ) from exc

try:
    from qiskit import QuantumCircuit, QuantumRegister
except ImportError as exc:
    raise ImportError(
        "Qiskit is required for quantum circuit construction. "
        "Install it with `pip install qiskit`."
    ) from exc

from ..hamiltonians.fermionic_hamiltonian import build_n2_fermionic_operator
from ..hamiltonians.qubit_hamiltonian import map_to_qubit_operator
from .initial_parameters import get_ccsd_amplitudes

if TYPE_CHECKING:
    from openfermion.ops import QubitOperator


def build_ucj_circuit(
    *,
    bond_length: float = 1.0977,
    basis: str = "6-31g",
    charge: int = 0,
    spin: int = 0,
    n_reps: int = 2,
    active_space: tuple[int, int] = (10, 8),
    unit: str = "Angstrom",
    symmetry: bool | str = "Dooh",
) -> tuple[QuantumCircuit, QubitOperator]:
    """Build a UCJ ansatz quantum circuit for N2 molecule.

    Args:
        bond_length: N–N bond length.
        basis: Basis set (default 6-31G).
        charge: Total molecular charge.
        spin: 2 * S (0 for singlet).
        n_reps: Number of UCJ operator repetitions (layers).
        active_space: (n_active_electrons, n_active_orbitals).
        unit: Distance unit for geometry.
        symmetry: Symmetry group.

    Returns:
        Tuple of (quantum_circuit, qubit_hamiltonian) where:
            - quantum_circuit: Qiskit circuit with UCJ ansatz
            - qubit_hamiltonian: QubitOperator for energy measurements
    """
    n_active_electrons, num_orbitals = active_space

    # Get CCSD amplitudes for initial parameters
    t1, t2 = get_ccsd_amplitudes(
        bond_length=bond_length,
        basis=basis,
        charge=charge,
        spin=spin,
        unit=unit,
        symmetry=symmetry,
    )

    # Define electron numbers (assuming spin-balanced system)
    num_elec_a = n_active_electrons // 2
    num_elec_b = n_active_electrons // 2
    nelec = (num_elec_a, num_elec_b)

    # Define interaction pairs for UCJ operator
    # Alpha-alpha: pairs of adjacent orbitals
    alpha_alpha_indices = [(p, p + 1) for p in range(num_orbitals - 1)]
    # Alpha-beta: pairs at specific positions
    alpha_beta_indices = [(p, p) for p in range(0, num_orbitals, 4)]

    # Create UCJ operator from CCSD amplitudes
    ucj_op = ffsim.UCJOpSpinBalanced.from_t_amplitudes(
        t2=t2,
        t1=t1,
        n_reps=n_reps,
        interaction_pairs=(alpha_alpha_indices, alpha_beta_indices),
    )

    # Create quantum circuit
    qubits = QuantumRegister(2 * num_orbitals, name="q")
    circuit = QuantumCircuit(qubits)

    # Prepare Hartree-Fock state as the reference state
    prep_gate = ffsim.qiskit.PrepareHartreeFockJW(num_orbitals, nelec)
    circuit.append(prep_gate, qubits)

    # Apply the UCJ operator to the reference state
    ucj_gate = ffsim.qiskit.UCJOpSpinBalancedJW(ucj_op)
    circuit.append(ucj_gate, qubits)

    # Optimize circuit by merging adjacent blocks
    circuit = ffsim.qiskit.PRE_INIT.optimize(circuit)

    # Add measurements
    circuit.measure_all()

    # Build the qubit Hamiltonian
    fermion_op = build_n2_fermionic_operator(
        bond_length=bond_length,
        basis=basis,
        charge=charge,
        spin=spin,
        active_space=active_space,
        unit=unit,
        symmetry=symmetry,
    )
    qubit_hamiltonian = map_to_qubit_operator(fermion_op)

    return circuit, qubit_hamiltonian


if __name__ == "__main__":
    # Build the UCJ circuit and Hamiltonian for N2
    print("Building UCJ ansatz circuit for N2...")
    circuit, hamiltonian = build_ucj_circuit()

    print(f"\nCircuit depth: {circuit.depth()}")
    print(f"Number of qubits: {circuit.num_qubits}")
    print(f"Number of operations: {len(circuit)}")
    print(f"\nHamiltonian terms: {len(hamiltonian.terms)}")
    
    # Optionally visualize circuit (uncomment to use)
    # circuit.decompose().draw("mpl", scale=0.5, fold=-1)
