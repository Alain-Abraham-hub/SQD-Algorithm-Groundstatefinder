#!/usr/bin/env python
"""Visualize the UCJ ansatz circuit for N2."""

from src.ansatz.LUCJ_ansatz import build_ucj_circuit

# Build the circuit
print("Building UCJ ansatz circuit for N2...")
circuit, hamiltonian = build_ucj_circuit()

# Text representation
print("\n" + "="*70)
print("CIRCUIT TEXT REPRESENTATION")
print("="*70)
print(circuit.decompose())

# ASCII art representation (more compact)
print("\n" + "="*70)
print("CIRCUIT ASCII (DECOMPOSED)")
print("="*70)
print(circuit.decompose().draw(fold=-1))

# Try to save as image if matplotlib is available
try:
    circuit_image = circuit.decompose().draw("mpl", scale=0.5, fold=-1)
    circuit_image.savefig("circuit_visualization.png", dpi=150, bbox_inches="tight")
    print("\n✓ Circuit visualization saved to 'circuit_visualization.png'")
except ImportError:
    print("\n⚠ Matplotlib not available; skipping image save")

# Print circuit statistics
print("\n" + "="*70)
print("CIRCUIT STATISTICS")
print("="*70)
print(f"Number of qubits: {circuit.num_qubits}")
print(f"Circuit depth: {circuit.depth()}")
print(f"Number of operations: {len(circuit)}")
print(f"Number of Hamiltonian terms: {len(hamiltonian.terms)}")
