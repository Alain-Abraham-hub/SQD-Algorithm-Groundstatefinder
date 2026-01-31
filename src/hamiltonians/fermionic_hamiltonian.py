"""Fermionic Hamiltonians built from ab-initio integrals."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Tuple

import numpy as np

try:
	from pyscf import ao2mo, gto, mcscf, scf
except ImportError as exc:  # pragma: no cover - handled by runtime env
	raise ImportError(
		"PySCF is required to build the N2 fermionic operator. "
		"Install it with `pip install pyscf`."
	) from exc

try:
	from openfermion.ops import InteractionOperator
	from openfermion.transforms import get_fermion_operator
except ImportError as exc:  # pragma: no cover - handled by runtime env
	raise ImportError(
		"OpenFermion is required to build the FermionOperator. "
		"Install it with `pip install openfermion`."
	) from exc


@dataclass(frozen=True)
class MoleculeSpec:
	"""Molecule specification used for integral generation."""

	bond_length: float = 1.0977
	basis: str = "6-31g"
	charge: int = 0
	spin: int = 0
	unit: str = "Angstrom"
	symmetry: bool | str = "Dooh"


def _build_n2_geometry(bond_length: float) -> str:
	half = bond_length / 2.0
	return f"N 0.0 0.0 {-half}; N 0.0 0.0 {half}"


def _expand_spatial_integrals_to_spin_orbital(
	h1: np.ndarray, h2: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
	"""Expand spatial integrals to spin-orbital tensors.

	Assumes `h2` is in chemist notation (pq|rs). The returned two-body
	tensor is in the convention used by OpenFermion's InteractionOperator:
	1/2 * sum_{pqrs} g[p,q,r,s] a_p^† a_q^† a_r a_s.
	"""

	n_orb = h1.shape[0]
	n_so = 2 * n_orb

	one_body = np.zeros((n_so, n_so))
	two_body = np.zeros((n_so, n_so, n_so, n_so))

	if h2.ndim == 2:
		h2 = ao2mo.restore(1, h2, n_orb)

	# Swap r and s to convert chemist notation into OpenFermion ordering.
	h2_phys = np.transpose(h2, (0, 1, 3, 2))

	for p in range(n_orb):
		for q in range(n_orb):
			value = h1[p, q]
			for spin in (0, 1):
				one_body[2 * p + spin, 2 * q + spin] = value

	for p in range(n_orb):
		for q in range(n_orb):
			for r in range(n_orb):
				for s in range(n_orb):
					value = h2_phys[p, q, r, s]
					if value == 0.0:
						continue
					for spin_p in (0, 1):
						for spin_q in (0, 1):
							two_body[
								2 * p + spin_p,
								2 * q + spin_q,
								2 * r + spin_p,
								2 * s + spin_q,
							] = value

	return one_body, two_body


def build_n2_fermionic_operator(
	*,
	bond_length: float = 1.0977,
	basis: str = "6-31g",
	charge: int = 0,
	spin: int = 0,
	active_space: Tuple[int, int] = (10, 8),
	unit: str = "Angstrom",
	symmetry: bool | str = "Dooh",
) -> "FermionOperator":
	"""Build a FermionOperator for N2 using CASCI integrals.

	Args:
		bond_length: N–N bond length.
		basis: Basis set (default 6-31G).
		charge: Total molecular charge.
		spin: 2 * S (0 for singlet).
		active_space: (n_active_electrons, n_active_orbitals).
		unit: Distance unit for geometry.
		symmetry: Enable/disable symmetry in PySCF.

	Returns:
		FermionOperator for the chosen active space.
	"""

	n_active_electrons, n_active_orbitals = active_space
	if n_active_electrons <= 0 or n_active_orbitals <= 0:
		raise ValueError("active_space must contain positive values.")

	mol_spec = MoleculeSpec(
		bond_length=bond_length,
		basis=basis,
		charge=charge,
		spin=spin,
		unit=unit,
		symmetry=symmetry,
	)

	mol = gto.Mole()
	mol.atom = _build_n2_geometry(bond_length)
	mol.basis = basis
	mol.charge = charge
	mol.spin = spin
	mol.unit = unit
	mol.symmetry = symmetry
	mol.build()

	mf = scf.RHF(mol) if spin == 0 else scf.ROHF(mol)
	mf.kernel()

	mc = mcscf.CASCI(mf, n_active_orbitals, n_active_electrons)

	h1eff = mc.get_h1eff()
	if isinstance(h1eff, tuple):
		h1eff = h1eff[0]

	h2eff = mc.get_h2eff()

	h1eff = np.asarray(h1eff)
	h2eff = np.asarray(h2eff)

	one_body_so, two_body_so = _expand_spatial_integrals_to_spin_orbital(
		h1eff, h2eff
	)

	interaction = InteractionOperator(
		constant=0.0,
		one_body_tensor=one_body_so,
		two_body_tensor=two_body_so,
	)
	fermion_operator = get_fermion_operator(interaction)

	return fermion_operator
