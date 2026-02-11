#!/usr/bin/env python
"""Generate a clean image of the UCJ circuit without text clutter."""

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

from src.ansatz.LUCJ_ansatz import build_ucj_circuit

print("Building UCJ ansatz circuit for N2...")
circuit, _ = build_ucj_circuit(build_hamiltonian=False)

print(f"Circuit depth: {circuit.depth()}")
print(f"Number of qubits: {circuit.num_qubits}")
print(f"Number of operations: {len(circuit)}")

# Generate the figure
print("\nGenerating circuit visualization...")
fig = circuit.decompose().draw(
    output="mpl",
    scale=0.6,
    fold=25,
    style={'backgroundcolor': '#FFFFFF'}
)

# Save the figure
output_file = "n2_ucj_circuit.png"
fig.savefig(output_file, dpi=300, bbox_inches="tight", facecolor='white')
plt.close(fig)

print(f"✓ Circuit visualization saved to '{output_file}'")
print("  Open the file to view a clean diagram of your quantum circuit.")
