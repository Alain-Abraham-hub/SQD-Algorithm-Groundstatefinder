# SQD-Algorithm-Groundstatefinder

A quantum framework for finding molecular ground states using the **LUCJ ansatz** and **SQD method**. Combines quantum circuit execution on IBM Quantum hardware with classical post-processing to estimate ground state energies.

## Key Features

- LUCJ ansatz with CCSD parameter initialization
- IBM Quantum backend integration
- Subspace Quantum Deflation (SQD) for configuration recovery
- PySCF for ab-initio electronic structure
- Modular architecture (Hamiltonians, Ansatz, Measurement)

## Quick Start

Docker is the primary way to run this project.

### Option 1: Use the pre-built image

```bash
docker pull alain1435/sqd-groundstatefinder:latest
```

Run a Python command inside the container:

```bash
docker run --rm -it \
	-v "$(pwd)":/app \
	alain1435/sqd-groundstatefinder:latest \
	python
```

### Option 2: Use Docker Compose

```bash
docker compose run --rm sqd-groundstatefinder
```

Start the notebook service:

```bash
docker compose up jupyter
```

### Option 3: Build locally

```bash
git clone <repository-url>
cd SQD-Algorithm-Groundstatefinder
docker build -t sqd-groundstatefinder .
docker run --rm -it -v "$(pwd)":/app sqd-groundstatefinder python
```

## Local Python Installation

If you do not want to use Docker:

```bash
git clone <repository-url>
cd SQD-Algorithm-Groundstatefinder
pip install -r requirements.txt
```

## Configuration

For IBM Quantum execution, provide these environment variables:

```bash
export QISKIT_IBM_TOKEN="your_token"
export QISKIT_IBM_INSTANCE="ibm-q/open/main"
```

With Docker:

```bash
docker run --rm \
	-v "$(pwd)":/app \
	-e QISKIT_IBM_TOKEN="your_token" \
	-e QISKIT_IBM_INSTANCE="ibm-q/open/main" \
	alain1435/sqd-groundstatefinder:latest \
	python your_script.py
```

## Usage

### Run a script with Docker

```bash
docker run --rm \
	-v "$(pwd)":/app \
	alain1435/sqd-groundstatefinder:latest \
	python visualize_circuit.py
```

### Run the notebook server

```bash
docker compose up jupyter
```

Then open `http://localhost:8888` in your browser.

### Python API example

Build a circuit and run on quantum hardware:

```python
from src.ansatz.LUCJ_ansatz import build_ucj_circuit
from src.ansatz.run_on_ibm import run_ansatz_on_ibm

circuit, hamiltonian = build_ucj_circuit(bond_length=1.0977, n_reps=2)
counts = run_ansatz_on_ibm(circuit, shots=1024)
```

Recover the ground state using SQD:

```python
from src.measure.subspace_hamiltonian import run_sqd_from_counts

energies, occupancies, _ = run_sqd_from_counts("ibm_counts.json", iterations=5)
print(f"Ground state energy: {energies[-1]:.6f}")
```

## Project Structure

```
SQD-Algorithm-Groundstatefinder/
├── Dockerfile
├── docker-compose.yml
├── DOCKER.md
├── README.md
├── LICENSE
├── requirements.txt
├── notebooks/
│   └── intro.ipynb              # Tutorial and examples
├── src/
│   ├── __init.py__
│   ├── hamiltonians/
│   │   ├── fermionic_hamiltonian.py    # Build H from ab-initio
│   │   ├── qubit_hamiltonian.py        # Transform to qubits
│   │   └── qubitmapping.py             # Qubit mapping helpers
│   ├── ansatz/
│   │   ├── LUCJ_ansatz.py              # Circuit design
│   │   ├── initial_parameters.py       # CCSD initialization
│   │   ├── run_on_ibm.py               # Execute on quantum
│   │   └── generate_circuit_image.py   # Visualization
│   └── measure/
│       └── subspace_hamiltonian.py     # SQD algorithm
└── ibm_counts.json              # Sample measurement results
```

## Docker Files

- `Dockerfile`: container image definition with scientific Python dependencies
- `docker-compose.yml`: interactive Python and Jupyter services
- `DOCKER.md`: detailed Docker usage and troubleshooting guide

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

For container-specific usage details, see [DOCKER.md](DOCKER.md).

## References

- [Qiskit Documentation](https://docs.quantum.ibm.com/)
- [OpenFermion](https://quantumai.google/reference/python/openfermion)
- [PySCF](https://pyscf.org/)
- [qiskit-addon-sqd](https://github.com/Qiskit/qiskit-addon-sqd)

## License

Licensed under the MIT License - see [LICENSE](LICENSE)

## Contributing

Contributions welcome! Please submit issues and pull requests.

