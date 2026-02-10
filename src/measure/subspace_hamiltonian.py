"""SQD configuration recovery using IBM bitstrings and a fermionic operator."""

from __future__ import annotations

import argparse
import json
from typing import Iterable

import numpy as np

from qiskit_addon_sqd.configuration_recovery import recover_configurations
from qiskit_addon_sqd.fermion import bitstring_matrix_to_ci_strs, solve_fermion
from qiskit_addon_sqd.subsampling import postselect_and_subsample

from openfermion.transforms import get_interaction_operator

from ..hamiltonians.fermionic_hamiltonian import build_n2_fermionic_operator


def _load_counts(path: str) -> tuple[dict[str, int], int, int]:
	with open(path, "r") as handle:
		payload = json.load(handle)

	counts = payload.get("counts", {})
	if not counts:
		raise ValueError("Counts file is empty or missing 'counts'.")

	num_qubits = int(payload.get("num_qubits", 0))
	if num_qubits <= 0:
		first_key = next(iter(counts))
		num_qubits = len(str(first_key).replace(" ", ""))

	shots = int(payload.get("shots", sum(int(v) for v in counts.values())))
	return counts, num_qubits, shots


def _counts_to_bitstring_matrix(
	counts: dict[str, int],
	num_qubits: int,
) -> tuple[np.ndarray, np.ndarray]:
	bitstrings = []
	weights = []
	total = sum(int(v) for v in counts.values())
	if total <= 0:
		raise ValueError("Total shots must be positive.")

	for bitstring, count in counts.items():
		clean = str(bitstring).replace(" ", "")
		if len(clean) != num_qubits:
			raise ValueError("Bitstring length does not match num_qubits.")
		bitstrings.append([1 if ch == "1" else 0 for ch in clean])
		weights.append(int(count) / total)

	return np.array(bitstrings, dtype=int), np.array(weights, dtype=float)


def _interaction_from_fermion_operator(num_spin_orbitals: int):
	fermion_operator = build_n2_fermionic_operator()
	return get_interaction_operator(fermion_operator, n_orbitals=num_spin_orbitals)


def run_sqd_from_counts(
	counts_file: str,
	*,
	iterations: int = 5,
	n_batches: int = 5,
	samples_per_batch: int = 500,
	max_davidson_cycles: int = 300,
	active_space: tuple[int, int] = (10, 8),
	spin: int = 0,
	rand_seed: int = 24,
) -> tuple[np.ndarray, np.ndarray, list[tuple[np.ndarray, np.ndarray]]]:
	counts, num_qubits, _ = _load_counts(counts_file)
	bitstring_matrix_full, probs_arr_full = _counts_to_bitstring_matrix(
		counts,
		num_qubits,
	)

	n_active_electrons, num_orbitals = active_space
	num_elec_a = n_active_electrons // 2
	num_elec_b = n_active_electrons // 2

	interaction = _interaction_from_fermion_operator(2 * num_orbitals)
	hcore = interaction.one_body_tensor
	eri = interaction.two_body_tensor
	nuclear_repulsion_energy = float(interaction.constant)

	spin_value = spin / 2.0
	spin_sq = spin_value * (spin_value + 1.0)
	open_shell = spin != 0

	rng = np.random.default_rng(rand_seed)

	e_hist = np.zeros((iterations, n_batches))
	s_hist = np.zeros((iterations, n_batches))
	occupancy_hist: list[tuple[np.ndarray, np.ndarray]] = []
	avg_occupancy = None

	for i in range(iterations):
		print(f"Starting configuration recovery iteration {i}")
		if avg_occupancy is None:
			bs_mat_tmp = bitstring_matrix_full
			probs_arr_tmp = probs_arr_full
		else:
			bs_mat_tmp, probs_arr_tmp = recover_configurations(
				bitstring_matrix_full,
				probs_arr_full,
				avg_occupancy,
				num_elec_a,
				num_elec_b,
				rand_seed=rng,
			)

		batches = postselect_and_subsample(
			bs_mat_tmp,
			probs_arr_tmp,
			hamming_right=num_elec_a,
			hamming_left=num_elec_b,
			samples_per_batch=samples_per_batch,
			num_batches=n_batches,
			rand_seed=rng,
		)

		e_tmp = np.zeros(n_batches)
		s_tmp = np.zeros(n_batches)
		occs_tmp = []

		for j in range(n_batches):
			strs_a, strs_b = bitstring_matrix_to_ci_strs(batches[j])
			print(f"  Batch {j} subspace dimension: {len(strs_a) * len(strs_b)}")
			energy_sci, coeffs_sci, avg_occs, spin_val = solve_fermion(
				batches[j],
				hcore,
				eri,
				open_shell=open_shell,
				spin_sq=spin_sq,
				max_davidson=max_davidson_cycles,
			)
			_ = coeffs_sci
			energy_sci += nuclear_repulsion_energy
			e_tmp[j] = energy_sci
			s_tmp[j] = spin_val
			occs_tmp.append(avg_occs)

		avg_occupancy = tuple(np.mean(occs_tmp, axis=0))
		e_hist[i, :] = e_tmp
		s_hist[i, :] = s_tmp
		occupancy_hist.append(avg_occupancy)

	return e_hist, s_hist, occupancy_hist


def main() -> None:
	parser = argparse.ArgumentParser(
		description="Run SQD configuration recovery from IBM bitstring counts."
	)
	parser.add_argument("--counts-file", type=str, default="ibm_counts.json")
	parser.add_argument("--iterations", type=int, default=5)
	parser.add_argument("--n-batches", type=int, default=5)
	parser.add_argument("--samples-per-batch", type=int, default=500)
	parser.add_argument("--max-davidson", type=int, default=300)
	args = parser.parse_args()

	e_hist, s_hist, _ = run_sqd_from_counts(
		args.counts_file,
		iterations=args.iterations,
		n_batches=args.n_batches,
		samples_per_batch=args.samples_per_batch,
		max_davidson_cycles=args.max_davidson,
	)

	print("Energy history (per iteration):")
	print(e_hist)
	print("Spin history (per iteration):")
	print(s_hist)


if __name__ == "__main__":
	main()
