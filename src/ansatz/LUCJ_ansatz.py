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

from ..hamiltonians.qubit_hamiltonian import build_qubit_hamiltonian
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
    build_hamiltonian: bool = True,
) -> tuple[QuantumCircuit, QubitOperator | None]:
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
        build_hamiltonian: If False, skip Hamiltonian building (faster, for visualization only).

    Returns:
        Tuple of (quantum_circuit, qubit_hamiltonian) where:
            - quantum_circuit: Qiskit circuit with UCJ ansatz
            - qubit_hamiltonian: QubitOperator for energy measurements (None if build_hamiltonian=False)
    """
    n_active_electrons, num_orbitals = active_space

    # Get CCSD amplitudes for initial parameters
    t1, t2 = get_ccsd_amplitudes(
        bond_length=bond_length,
        basis=basis,
        charge=charge,
        spin=spin,
        active_space=active_space,
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
    circuit = ffsim.qiskit.PRE_INIT.run(circuit)

    # Add measurements
    circuit.measure_all()

    # Build the qubit Hamiltonian (can be slow for large systems)
    qubit_hamiltonian = None
    if build_hamiltonian:
        qubit_hamiltonian = build_qubit_hamiltonian(
            bond_length=bond_length,
            basis=basis,
            charge=charge,
            spin=spin,
            active_space=active_space,
            unit=unit,
            symmetry=symmetry,
        )

    return circuit, qubit_hamiltonian


def draw_circuit(
    circuit: QuantumCircuit,
    output: str = "mpl",
    filename: str = "circuit.png",
    decompose: bool = True,
    fold: int = 20,
    scale: float = 0.7,
) -> None:
    """Draw the quantum circuit in various formats.

    Args:
        circuit: Qiskit QuantumCircuit to visualize.
        output: Output format:
            - "mpl": matplotlib image (recommended, saves to file)
            - "text": ASCII art in terminal
            - "latex": LaTeX source code
        filename: Output filename (for "mpl" or "latex" modes).
        decompose: Whether to decompose gates to primitive operations.
        fold: Wrap circuit diagram at this many columns (-1 for no folding).
        scale: Scale factor for matplotlib drawings.
    """
    circuit_to_draw = circuit.decompose() if decompose else circuit
    
    if output == "text":
        print("\n" + "="*70)
        print("CIRCUIT (ASCII)")
        print("="*70)
        print(circuit_to_draw.draw(output='text', fold=fold))
        
    elif output == "mpl":
        try:
            import matplotlib.pyplot as plt
            import matplotlib
            matplotlib.use('Agg')  # Use non-interactive backend
            import sys
            import io
            
            # Suppress text output from draw()
            old_stdout = sys.stdout
            sys.stdout = io.StringIO()
            
            try:
                fig = circuit_to_draw.draw(
                    output="mpl",
                    scale=scale,
                    fold=fold,
                    style={'backgroundcolor': '#EEEEEE'}
                )
            finally:
                sys.stdout = old_stdout
            
            if filename:
                fig.savefig(filename, dpi=300, bbox_inches="tight", facecolor='white')
                print(f"✓ Circuit visualization saved to '{filename}'")
                print(f"  Open the file to view the circuit diagram.")
            plt.close(fig)
        except ImportError:
            print("⚠ Matplotlib not installed. Install with: pip install matplotlib")
            print("  Falling back to text mode...")
            draw_circuit(circuit, output="text", fold=fold)
            
    elif output == "latex":
        try:
            latex_source = circuit_to_draw.draw(output='latex_source')
            if filename:
                with open(filename, 'w') as f:
                    f.write(latex_source)
                print(f"✓ LaTeX source saved to '{filename}'")
                print(f"  Compile with pdflatex or use in your LaTeX document")
            else:
                print(latex_source)
        except Exception as e:
            print(f"⚠ Error generating LaTeX: {e}")
    
    else:
        print(f"⚠ Unknown output format '{output}'. Use 'text', 'mpl', or 'latex'.")


if __name__ == "__main__":
    # Build the UCJ circuit (skip Hamiltonian for faster visualization)
    print("Building UCJ ansatz circuit for N2...")
    circuit, hamiltonian = build_ucj_circuit(build_hamiltonian=False)

    print(f"\nCircuit depth: {circuit.depth()}")
    print(f"Number of qubits: {circuit.num_qubits}")
    print(f"Number of operations: {len(circuit)}")
    
    # Draw and save the circuit as an image (recommended)
    print("\nGenerating circuit visualization...")
    draw_circuit(circuit, output="mpl", filename="n2_ucj_circuit.png", fold=20, scale=0.8)
    
    # To also build the Hamiltonian (slow), use:
    # circuit, hamiltonian = build_ucj_circuit(build_hamiltonian=True)
    # print(f"\nHamiltonian terms: {len(hamiltonian.terms)}")
