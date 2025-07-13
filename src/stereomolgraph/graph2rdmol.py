# pyright: standard
from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import rdkit.Chem as Chem  # type: ignore

from stereomolgraph.algorithms.bond_orders import connectivity2bond_orders
from stereomolgraph.stereodescriptors import (
    AtropBond,
    Octahedral,
    PlanarBond,
    SquarePlanar,
    Tetrahedral,
    TrigonalBipyramidal,
)

if TYPE_CHECKING:
    from stereomolgraph.graphs.crg import CondensedReactionGraph
    from stereomolgraph.graphs.mg import MolGraph
    from stereomolgraph.graphs.smg import StereoMolGraph


bond_type_dict = {
    0.5: Chem.BondType.HYDROGEN,
    # to be drawn as a dotted line
    # (looks better than other options)
    0: Chem.BondType.UNSPECIFIED,
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE,
    4: Chem.BondType.QUADRUPLE,
    5: Chem.BondType.QUINTUPLE,
    6: Chem.BondType.HEXTUPLE,
    1.5: Chem.BondType.ONEANDAHALF,
    2.5: Chem.BondType.TWOANDAHALF,
    3.5: Chem.BondType.THREEANDAHALF,
    4.5: Chem.BondType.FOURANDAHALF,
    5.5: Chem.BondType.FIVEANDAHALF,
    "AROMATIC": Chem.BondType.AROMATIC,
    "IONIC": Chem.BondType.IONIC,
    "HYDROGEN": Chem.BondType.HYDROGEN,
    "THREECENTER": Chem.BondType.THREECENTER,
    "DATIVEONE": Chem.BondType.DATIVEONE,
    "DATIVE": Chem.BondType.DATIVE,
    "DATIVEL": Chem.BondType.DATIVEL,
    "DATIVER": Chem.BondType.DATIVER,
    "OTHER": Chem.BondType.OTHER,
    "ZERO": Chem.BondType.ZERO,
}


def set_bond_orders(
    graph: MolGraph,
    mol: Chem.rdchem.RWMol,
    idx_map_num_dict: dict[int, int],
    allow_charged_fragments=False,
    charge=0,
) -> Chem.rdchem.RWMol:
    bond_order_mat, atomic_charges, unpaired_electrons = (
        connectivity2bond_orders(
            atom_types=graph.atom_types,
            connectivity_matrix=graph.connectivity_matrix(),
            allow_charged_fragments=allow_charged_fragments,
            charge=charge,
        )
    )

    index_map_num_dict = {i: map_num for i, map_num in enumerate(graph.atoms)}

    map_num_idx_dict = {
        map_num: idx for idx, map_num in idx_map_num_dict.items()
    }

    for bond in graph.bonds:
        atom1, atom2 = bond
        bond_order = bond_order_mat[index_map_num_dict[atom1]][
            index_map_num_dict[atom2]
        ]

        mol.GetBondBetweenAtoms(
            map_num_idx_dict[atom1], map_num_idx_dict[atom2]
        ).SetBondType(bond_type_dict[bond_order])

    for i, atomic_charge in enumerate(atomic_charges):
        if atomic_charge:
            mol.GetAtomWithIdx(i).SetFormalCharge(int(atomic_charge))
    for i, unpaired_e in enumerate(unpaired_electrons):
        if unpaired_e:
            mol.GetAtomWithIdx(i).SetNumRadicalElectrons(unpaired_e)

    return mol


def mol_graph_to_rdmol(
    graph: MolGraph,
    generate_bond_orders=False,
    allow_charged_fragments=False,
    charge=0,
) -> tuple[Chem.rdchem.RWMol, dict[int, int]]:
    mol = Chem.RWMol()

    atom_types_strings = []
    idx_map_num_dict = {}

    for atom in graph.atoms:
        atom_type = graph.get_atom_attribute(atom, "atom_type")
        if atom_type is None:
            raise RuntimeError(atom, graph.atoms, graph.atom_types)
        rd_atom = Chem.Atom(atom_type.symbol)
        rd_atom.SetNoImplicit(True)
        atom_index = mol.AddAtom(rd_atom)
        idx_map_num_dict[atom_index] = atom
        atom_types_strings.append(atom_type.atomic_nr)
        try:
            mol.GetAtomWithIdx(atom_index).SetAtomMapNum(atom)
        except OverflowError as e:
            warnings.warn(f"Atom number number not set: {e}")

    for i in range(graph.n_atoms):
        for j in range(i + 1, graph.n_atoms):
            if graph.has_bond(idx_map_num_dict[i], idx_map_num_dict[j]):
                mol.AddBond(i, j)
                # TODO: check if this is still needed
                # mol.GetBondBetweenAtoms(i, j).SetBondType(
                #    Chem.rdchem.BondType.SINGLE)
    if generate_bond_orders:
        mol = set_bond_orders(
            graph=graph,
            mol=mol,
            idx_map_num_dict=idx_map_num_dict,
            allow_charged_fragments=allow_charged_fragments,
            charge=charge,
        )

    return mol, idx_map_num_dict


def stereo_mol_graph_to_rdmol(
    graph: StereoMolGraph,
    generate_bond_orders=False,
    allow_charged_fragments=False,
    charge=0,
) -> tuple[Chem.rdchem.RWMol, dict[int, int]]:
    """
    Creates a RDKit mol object using the connectivity of the mol graph.
    Stereochemistry is added to the mol object.

    :return: RDKit molecule
    :rtype: Chem.rdchem.RWMol
    """
    rd_tetrahedral = {
        1: Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
        -1: Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
        None: Chem.rdchem.ChiralType.CHI_TETRAHEDRAL,
    }
    mol, idx_map_num_dict = mol_graph_to_rdmol(
        graph,
        generate_bond_orders=generate_bond_orders,
        allow_charged_fragments=allow_charged_fragments,
        charge=charge,
    )

    map_num_idx_dict = {v: k for k, v in idx_map_num_dict.items()}

    for atom in graph.atoms:
        a_stereo = graph.get_atom_stereo(atom)
        atom_idx = map_num_idx_dict[atom]
        rd_atom = mol.GetAtomWithIdx(atom_idx)

        if a_stereo is not None and any(
            a not in graph.atoms for a in a_stereo.atoms
        ):
            raise NotImplementedError(
                "Handling of missing atoms not supported yet"
            )
            for mis_a in [
                a for a in a_stereo.atoms[1:] if a not in graph.atoms
            ]:
                if mis_a not in map_num_idx_dict:
                    # add dummy atom
                    idx = mol.AddAtom(Chem.Atom(0))
                    map_num_idx_dict[mis_a] = idx
                    idx_map_num_dict[idx] = mis_a
                    try:
                        mol.GetAtomWithIdx(idx).SetAtomMapNum(mis_a)
                    except OverflowError as e:
                        warnings.warn(f"Atom number number not set: {e}")

                mol.AddBond(
                    map_num_idx_dict[atom],
                    map_num_idx_dict[mis_a],
                    Chem.BondType.SINGLE,
                )

        # The chirality of an Atom in rdkit is determined by two things:
        # 1. its chiralTag
        # 2. the input order of its bonds
        # (see note below for handling of implicit Hs)
        #
        # For tetrahedral coordination, the chiralTag tells you what
        # direction you have to rotate to get from bond 2 to bond 3 while
        # looking down bond 1. This is pretty much identical to the SMILES
        # representation of chirality.

        if a_stereo is not None and isinstance(a_stereo, Tetrahedral):
            mol.GetAtomWithIdx(atom_idx).SetHybridization(
                Chem.HybridizationType.SP3
            )
            assert len(a_stereo.atoms) == 5
            rd_bonds = [
                idx_map_num_dict[a.GetIdx()]
                for a in mol.GetAtomWithIdx(atom_idx).GetNeighbors()
            ]
            perm = [
                p
                for p in a_stereo._perm_atoms()
                if p[1] == rd_bonds[0] and p[2] == rd_bonds[1]
            ][0]

            if (perm[3], perm[4]) == (rd_bonds[2], rd_bonds[3]):
                mol.GetAtomWithIdx(atom_idx).SetChiralTag(
                    rd_tetrahedral[a_stereo.parity]
                )

            elif (perm[3], perm[4]) == (rd_bonds[3], rd_bonds[2]):
                mol.GetAtomWithIdx(atom_idx).SetChiralTag(
                    rd_tetrahedral[a_stereo.parity]
                )
            else:
                raise RuntimeError(
                    "Central atom was not bonded to all Stereo Atoms"
                )

        elif a_stereo is not None and isinstance(a_stereo, SquarePlanar):
            rd_atom.SetChiralTag(Chem.ChiralType.CHI_SQUAREPLANAR)
            neighbors = tuple(
                [
                    idx_map_num_dict[atom.GetIdx()]
                    for atom in rd_atom.GetNeighbors()
                ]
            )

            if neighbors in {p[1:] for p in a_stereo._perm_atoms()}:
                rd_atom.SetUnsignedProp("_chiralPermutation", 1)
            else:
                rd_atom.SetUnsignedProp("_chiralPermutation", 2)

        elif a_stereo is not None and isinstance(
            a_stereo, TrigonalBipyramidal
        ):
            # the order of the neighbors is important for stereochemistry
            for rd_n in rd_atom.GetNeighbors():
                mol.RemoveBond(rd_n.GetIdx(), atom_idx)

            for a in (1, 3, 4, 5, 2):
                a = a_stereo.atoms[a]
                mol.AddBond(
                    atom_idx,
                    map_num_idx_dict[a],
                )
            rd_atom.SetHybridization(Chem.HybridizationType.SP3D)
            rd_atom.SetChiralTag(Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL)
            if a_stereo.parity == 1:
                rd_atom.SetUnsignedProp("_chiralPermutation", 1)
            elif a_stereo.parity == -1:
                rd_atom.SetUnsignedProp("_chiralPermutation", 2)

        elif a_stereo is not None and isinstance(a_stereo, Octahedral):
            for rd_n in rd_atom.GetNeighbors():
                mol.RemoveBond(rd_n.GetIdx(), atom_idx)

            for a in (1, 5, 6, 3, 4, 2):
                a = a_stereo.atoms[a]
                mol.AddBond(
                    atom_idx,
                    map_num_idx_dict[a],
                )
            rd_atom.SetChiralTag(Chem.ChiralType.CHI_OCTAHEDRAL)
            rd_atom.SetHybridization(Chem.HybridizationType.SP3D2)
            if a_stereo.parity == 1:
                rd_atom.SetUnsignedProp("_chiralPermutation", 1)
            elif a_stereo.parity == -1:
                rd_atom.SetUnsignedProp("_chiralPermutation", 2)

    for b_stereo in (bs for bs in graph.bond_stereo.values() if bs):
        a1, a2 = b_stereo.atoms[2], b_stereo.atoms[3]

        if not all(a in graph.atoms for a in b_stereo.atoms):
            raise NotImplementedError(
                "Handling of missing atoms not supported yet"
            )
            for mis_a in (a for a in b_stereo.atoms if a not in graph.atoms):
                if mis_a not in map_num_idx_dict:
                    # add dummy atom
                    idx = mol.AddAtom(Chem.Atom(0))
                    map_num_idx_dict[mis_a] = idx
                    idx_map_num_dict[idx] = mis_a
                    try:
                        mol.GetAtomWithIdx(idx).SetAtomMapNum(mis_a)
                    except OverflowError as e:
                        warnings.warn(f"Atom number number not set: {e}")

                if mis_a in b_stereo.atoms[:2]:
                    mol.AddBond(
                        map_num_idx_dict[b_stereo.atoms[2]],
                        map_num_idx_dict[mis_a],
                        Chem.BondType.SINGLE,
                    )
                elif mis_a in b_stereo.atoms[4:6]:
                    mol.AddBond(
                        map_num_idx_dict[b_stereo.atoms[3]],
                        map_num_idx_dict[mis_a],
                        Chem.BondType.SINGLE,
                    )
                else:
                    raise RuntimeError("This should not happen")

        rd_a1 = map_num_idx_dict[a1]
        rd_a2 = map_num_idx_dict[a2]
        rd_bond = mol.GetBondBetweenAtoms(rd_a1, rd_a2)
        new_a1 = idx_map_num_dict[rd_bond.GetBeginAtomIdx()]
        new_a2 = idx_map_num_dict[rd_bond.GetEndAtomIdx()]

        assert {a1, a2} == {new_a1, new_a2}

        if isinstance(b_stereo, PlanarBond):
            mol.GetAtomWithIdx(rd_a1).SetHybridization(
                Chem.HybridizationType.SP2
            )
            mol.GetAtomWithIdx(rd_a2).SetHybridization(
                Chem.HybridizationType.SP2
            )

            if b_stereo.parity is None:
                rd_bond.SetStereo(Chem.rdchem.BondStereo.STEREONONE)

            elif (a1, a2) == (new_a1, new_a2):
                rd_bond.SetStereoAtoms(
                    map_num_idx_dict[b_stereo.atoms[0]],
                    map_num_idx_dict[b_stereo.atoms[4]],
                )
                rd_bond.SetStereo(Chem.rdchem.BondStereo.STEREOZ)

            elif (a1, a2) == (new_a2, new_a1):
                rd_bond.SetStereoAtoms(
                    map_num_idx_dict[b_stereo.atoms[4]],
                    map_num_idx_dict[b_stereo.atoms[0]],
                )
                rd_bond.SetStereo(Chem.rdchem.BondStereo.STEREOZ)
            else:
                raise Exception(f"something wrong with {b_stereo}")

            # if no planar bond neigboring set the bond to double
            if False:
                ...
                # TODO: check if this is still needed
                # all(
                # graph.get_bond_stereo(
                #    (b_stereo.atoms[i], b_stereo.atoms[j])
                # )
                # is None
                # for i, j in ((0, 2), (1, 2), (3, 4), (3, 5))
                # if tuple(sorted((b_stereo.atoms[i], b_stereo.atoms[j])))
                #         in graph.bonds):
        #
        #     rd_bond.SetBondType(Chem.BondType.DOUBLE)

        elif isinstance(b_stereo, AtropBond):
            if (a1, a2) == (new_a1, new_a2):
                rd_bond.SetStereoAtoms(
                    map_num_idx_dict[b_stereo.atoms[0]],
                    map_num_idx_dict[b_stereo.atoms[4]],
                )
                if b_stereo.parity == 1:
                    rd_bond.SetStereo(Chem.rdchem.BondStereo.STEREOATROPCW)
                elif b_stereo.parity == -1:
                    rd_bond.SetStereo(Chem.rdchem.BondStereo.STEREOATROPCCW)

            elif (a1, a2) == (new_a2, new_a1):
                rd_bond.SetStereoAtoms(
                    map_num_idx_dict[b_stereo.atoms[4]],
                    map_num_idx_dict[b_stereo.atoms[0]],
                )
                if b_stereo.parity == 1:
                    rd_bond.SetStereo(Chem.rdchem.BondStereo.STEREOATROPCCW)
                elif b_stereo.parity == -1:
                    rd_bond.SetStereo(Chem.rdchem.BondStereo.STEREOATROPCW)
            else:
                raise Exception(f"something wrong with {b_stereo}")

    if generate_bond_orders:
        mol = set_bond_orders(
            graph=graph,
            mol=mol,
            idx_map_num_dict=idx_map_num_dict,
        )

    return mol, idx_map_num_dict


def set_crg_bond_orders(
    graph: CondensedReactionGraph,
    mol: Chem.rdchem.RWMol,
    idx_map_num_dict: dict[int, int],
    generate_bond_orders=False,
    allow_charged_fragments=False,
    charge=0,
) -> tuple[Chem.rdchem.RWMol, dict[int, int]]:
    r, r_idx_map_num_dict = mol_graph_to_rdmol(
        graph.reactant(),
        generate_bond_orders=True,
        allow_charged_fragments=allow_charged_fragments,
        charge=charge,
    )
    p, p_idx_map_num_dict = mol_graph_to_rdmol(
        graph.product(),
        generate_bond_orders=True,
        allow_charged_fragments=allow_charged_fragments,
        charge=charge,
    )
    assert r_idx_map_num_dict == p_idx_map_num_dict == idx_map_num_dict

    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (r_bond := r.GetBondBetweenAtoms(a1, a2)) and (
            p_bond := p.GetBondBetweenAtoms(a1, a2)
        ):
            r_bond_order = r_bond.GetBondTypeAsDouble()
            p_bond_order = p_bond.GetBondTypeAsDouble()
            average = (r_bond_order + p_bond_order) / 2
            bond_order = round(average * 2) / 2
            mol.GetBondBetweenAtoms(a1, a2).SetBondType(
                bond_type_dict[bond_order]
            )

    return mol
