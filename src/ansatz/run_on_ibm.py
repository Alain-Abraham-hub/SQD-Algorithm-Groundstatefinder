"""Submit the UCJ circuit to the least-busy IBM Quantum backend using SamplerV2.

This script expects credentials in env vars:
- QISKIT_IBM_TOKEN (required)
- QISKIT_IBM_INSTANCE (optional, e.g. "ibm-q/open/main")
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qiskit import transpile

try:
    from qiskit_ibm_runtime import QiskitRuntimeService
    from qiskit_ibm_runtime import SamplerV2 as Sampler
except ImportError as exc:
    raise ImportError(
        "qiskit-ibm-runtime is required to submit jobs. "
        "Install it with `pip install qiskit-ibm-runtime`."
    ) from exc

from .LUCJ_ansatz import build_ucj_circuit


def _select_least_busy_backend(
    backends: Iterable,
    min_qubits: int,
) -> object:
    candidates = []
    for backend in backends:
        try:
            if getattr(backend, "num_qubits", 0) < min_qubits:
                continue
            status = backend.status()
            if not status.operational:
                continue
            candidates.append((status.pending_jobs, backend.name, backend))
        except Exception:
            continue

    if not candidates:
        raise RuntimeError("No operational backend satisfies the qubit requirement.")

    candidates.sort(key=lambda item: (item[0], item[1]))
    return candidates[0][2]


def _format_bitstring(key: object, num_qubits: int) -> str:
    if isinstance(key, str):
        return key
    try:
        value = int(key)
    except (TypeError, ValueError):
        return str(key)
    return format(value, f"0{num_qubits}b")


def _quasi_to_counts(quasi: dict, num_qubits: int, shots: int) -> dict[str, int]:
    items = []
    weights = []
    for key, prob in quasi.items():
        if prob <= 0:
            continue
        items.append(_format_bitstring(key, num_qubits))
        weights.append(prob)

    if not items:
        return {}

    weights_arr = np.array(weights, dtype=float)
    weights_arr = weights_arr / weights_arr.sum()
    samples = np.random.choice(items, size=shots, p=weights_arr)
    counts = {}
    for sample in samples:
        counts[sample] = counts.get(sample, 0) + 1
    return counts


def _get_field(container: object, names: Iterable[str]) -> object | None:
    for name in names:
        if hasattr(container, name):
            return getattr(container, name)
        if isinstance(container, dict) and name in container:
            return container[name]
    return None


def _extract_counts_or_quasi(result: object, num_qubits: int, shots: int) -> tuple[dict[str, int] | None, dict | None]:
    pub = result[0]
    data = getattr(pub, "data", None)

    meas = _get_field(data, ["meas", "measure", "measurements"])
    if meas is not None:
        for method in ("get_counts", "to_counts"):
            if hasattr(meas, method):
                counts = getattr(meas, method)()
                if isinstance(counts, dict):
                    return counts, None
        if isinstance(meas, dict):
            return meas, None

    quasi = _get_field(data, ["quasi_dists", "quasi_dist", "quasi", "quasi_distribution"])
    if quasi is None:
        quasi = _get_field(pub, ["quasi_dists", "quasi_dist", "quasi", "quasi_distribution"])

    if isinstance(quasi, list):
        quasi = quasi[0] if quasi else None

    if isinstance(quasi, dict):
        return _quasi_to_counts(quasi, num_qubits, shots), quasi

    return None, None


def _plot_counts(counts: dict[str, int], max_bars: int, output_file: str) -> None:
    if not counts:
        print("No counts to plot.")
        return

    ordered = sorted(counts.items(), key=lambda item: item[1], reverse=True)
    if max_bars > 0:
        ordered = ordered[:max_bars]

    labels = [key for key, _ in ordered]
    values = [value for _, value in ordered]

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(labels, values, color="#2A6F97")
    ax.set_xlabel("Bitstring")
    ax.set_ylabel("Frequency")
    ax.set_title("Measurement Frequencies (Sampled)")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Frequency plot saved to '{output_file}'")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the UCJ circuit on IBM Quantum hardware with SamplerV2."
    )
    parser.add_argument("--shots", type=int, default=10_000)
    parser.add_argument("--opt-level", type=int, default=1)
    parser.add_argument("--backend", type=str, default="")
    parser.add_argument("--plot-file", type=str, default="bitstring_frequencies.png")
    parser.add_argument("--max-bars", type=int, default=30)
    parser.add_argument("--counts-file", type=str, default="ibm_counts.json")
    args = parser.parse_args()

    token = os.getenv("QISKIT_IBM_TOKEN")
    if not token:
        raise RuntimeError("Set QISKIT_IBM_TOKEN in your environment.")

    instance = os.getenv("QISKIT_IBM_INSTANCE")
    service = QiskitRuntimeService(
        channel="ibm_quantum_platform",
        token=token,
        instance=instance,
    )

    circuit, _ = build_ucj_circuit(build_hamiltonian=False)

    if args.backend:
        backend = service.backend(args.backend)
    else:
        backend = _select_least_busy_backend(
            service.backends(simulator=False),
            min_qubits=circuit.num_qubits,
        )

    isa_circuit = transpile(
        circuit,
        backend=backend,
        optimization_level=args.opt_level,
    )

    sampler = Sampler(mode=backend)
    sampler.options.dynamical_decoupling.enable = True

    job = sampler.run([isa_circuit], shots=args.shots)

    print(f"Backend: {backend.name}")
    print(f"Job ID: {job.job_id()}")

    result = job.result()
    counts, quasi = _extract_counts_or_quasi(result, circuit.num_qubits, args.shots)

    if quasi is not None:
        print("Quasi distribution:")
        print(quasi)

    if counts is None:
        print("No measurement data found in sampler result.")
        return

    with open(args.counts_file, "w") as handle:
        json.dump(
            {
                "num_qubits": circuit.num_qubits,
                "shots": args.shots,
                "counts": counts,
            },
            handle,
            indent=2,
            sort_keys=True,
        )
    print(f"Counts saved to '{args.counts_file}'")

    _plot_counts(counts, args.max_bars, args.plot_file)


if __name__ == "__main__":
    main()
