# SQD-Algorithm-Groundstatefinder

A quantum framework for finding molecular ground states using the **LUCJ ansatz** and **SQD method**. Combines quantum circuit execution on IBM Quantum hardware with classical post-processing to estimate ground state energies.

## Key Features

- LUCJ ansatz with CCSD parameter initialization
- IBM Quantum backend integration
- Subspace Quantum Deflation (SQD) for configuration recovery
- PySCF for ab-initio electronic structure
- Modular architecture (Hamiltonians, Ansatz, Measurement)

## Installation

```bash
# Clone and install
git clone <repository-url>
cd SQD-Algorithm-Groundstatefinder
pip install -r requirements.txt

# Optional: Configure IBM Quantum
export QISKIT_IBM_TOKEN="your_token"
export QISKIT_IBM_INSTANCE="ibm-q/open/main"
```

## Usage

Build circuit and run on quantum hardware:
```python
from src.ansatz.LUCJ_ansatz import build_ucj_circuit
from src.ansatz.run_on_ibm import run_ansatz_on_ibm

circuit, hamiltonian = build_ucj_circuit(bond_length=1.0977, n_reps=2)
counts = run_ansatz_on_ibm(circuit, shots=1024)
```

Recover ground state using SQD:
```python
from src.measure.subspace_hamiltonian import run_sqd_from_counts

energies, occupancies, _ = run_sqd_from_counts("ibm_counts.json", iterations=5)
print(f"Ground state energy: {energies[-1]:.6f}")
```

## Project Structure

```
SQD-Algorithm-Groundstatefinder/
├── README.md
├── LICENSE
├── requirements.txt
├── notebooks/
│   └── intro.ipynb              # Tutorial and examples
├── src/
│   ├── __init__.py
│   ├── hamiltonians/
│   │   ├── fermionic_hamiltonian.py    # Build H from ab-initio
│   │   └── qubit_hamiltonian.py        # Transform to qubits
│   ├── ansatz/
│   │   ├── LUCJ_ansatz.py              # Circuit design
│   │   ├── initial_parameters.py       # CCSD initialization
│   │   ├── run_on_ibm.py               # Execute on quantum
│   │   └── generate_circuit_image.py   # Visualization
│   └── measure/
│       └── subspace_hamiltonian.py     # SQD algorithm
└── ibm_counts.json              # Sample measurement results
```

## Algorithm Overview

```
1. Build molecular Hamiltonian (PySCF) → Fermionic operator
2. Transform to qubit Hamiltonian (Jordan-Wigner)
3. Prepare LUCJ ansatz circuit (CCSD parameters)
4. Execute on IBM Quantum hardware
5. Process bitstring counts with SQD method
6. Recover ground state energy
```

## Key Concepts

**LUCJ Ansatz**: Parameterized quantum circuit combining unitary coupled-cluster and Jastrow factors, initialized from classical CCSD.

**SQD Method**: Post-processing algorithm that recovers electronic configurations from quantum bitstrings and refines ground state estimates.

**Active Space**: Approximation where only the most important electrons receive quantum treatment (e.g., 10 electrons in 8 orbitals for N₂).

## Tutorial

See [notebooks/intro.ipynb](notebooks/intro.ipynb) for detailed examples and tutorials.

## References

- [Qiskit Documentation](https://docs.quantum.ibm.com/)
- [OpenFermion](https://quantumai.google/reference/python/openfermion)
- [PySCF](https://pyscf.org/)
- [qiskit-addon-sqd](https://github.com/Qiskit/qiskit-addon-sqd)

## License

Licensed under the MIT License - see [LICENSE](LICENSE)

## Contributing

Contributions welcome! Please submit issues and pull requests.

