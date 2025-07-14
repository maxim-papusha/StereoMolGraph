# pyright: standard
# typing with rdkit is not fully supported
from __future__ import annotations

from typing import TYPE_CHECKING

import rdkit.Chem as Chem # type: ignore

from stereomolgraph.stereodescriptors import (
    Tetrahedral,
    SquarePlanar,
    TrigonalBipyramidal,
    Octahedral,
    PlanarBond,
    AtomStereo,
)

if TYPE_CHECKING:
    from stereomolgraph.graphs.mg import MolGraph
    from stereomolgraph.graphs.smg import StereoMolGraph


def mol_graph_from_rdmol(
    cls: type[MolGraph], rdmol: Chem.Mol, use_atom_map_number:bool=False
) -> MolGraph:
    """
    Creates a StereoMolGraph from an RDKit Mol object.
    Implicit Hydrogens are added to the graph.
    Stereo information is conserved. Double bonds, aromatic bonds and
    conjugated bonds are interpreted as planar. Atoms with 5 bonding
    partners are assumed to be TrigonalBipyramidal and allow interchange
    of the substituents (berry pseudorotation). Atoms with 6 bonding
    partners are assumed to be octahedral and do not allow interchange of
    the substituents.

    :param rdmol: RDKit Mol object
    :param use_atom_map_number: If the atom map number should be used
                                instead of the atom index
    :return: StereoMolGraph
    """
    # rdmol = Chem.AddHs(rdmol, explicitOnly=True, addCoords=True)

    if use_atom_map_number is False:
        rdmol = Chem.rdmolops.AddHs(rdmol, explicitOnly=True)

    graph = cls()

    if use_atom_map_number:
        id_atom_map = {atom.GetIdx(): atom.GetAtomMapNum() for atom in rdmol.GetAtoms()}
    else:
        id_atom_map = {atom.GetIdx(): atom.GetIdx() for atom in rdmol.GetAtoms()}

    for atom in rdmol.GetAtoms():
        graph.add_atom(id_atom_map[atom.GetIdx()], atom.GetSymbol())

    for bond in rdmol.GetBonds():
        graph.add_bond(
            id_atom_map[bond.GetBeginAtomIdx()],
            id_atom_map[bond.GetEndAtomIdx()],
        )
    return graph


def stereo_mol_graph_from_rdmol(
    cls: type[StereoMolGraph], rdmol: Chem.Mol, use_atom_map_number:bool=False
) -> StereoMolGraph:
    """
    Creates a StereoMolGraph from an RDKit Mol object.
    All hydrogens have to be explicit.
    Stereo information is conserved for tetrahedral atoms and
    double bonds.

    :param rdmol: RDKit Mol object
    :param use_atom_map_number: If the atom map number should be used
                                instead of the atom index, Default: False
    :return: StereoMolGraph
    """
    rd_tetrahedral = {
        Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW: -1,
        Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW: 1,
        Chem.rdchem.ChiralType.CHI_TETRAHEDRAL: None,
    }

    rdmol = Chem.AddHs(rdmol, explicitOnly=False)

    if use_atom_map_number is True:
        if any(atom.GetAtomMapNum() == 0 for atom in rdmol.GetAtoms()):
            raise ValueError("AtomMapNumber has to  be set on all atoms")
        id_atom_map: dict[int, int] = {
            atom.GetIdx(): atom.GetAtomMapNum() for atom in rdmol.GetAtoms()
        }
    else:
        id_atom_map: dict[int, int] = {
            atom.GetIdx(): atom.GetIdx() for atom in rdmol.GetAtoms()
        }

    graph = cls()

    for atom in rdmol.GetAtoms():
        graph.add_atom(id_atom_map[atom.GetIdx()], atom.GetSymbol())

    for bond in rdmol.GetBonds():
        graph.add_bond(
            id_atom_map[bond.GetBeginAtomIdx()],
            id_atom_map[bond.GetEndAtomIdx()],
        )

    for atom in rdmol.GetAtoms():
        atom_idx: int = atom.GetIdx()

        neighbors: tuple[int, ...] = tuple(
            [
                (
                    {b.GetBeginAtomIdx(), b.GetEndAtomIdx()}
                    - {
                        atom_idx,
                    }
                ).pop()
                for b in rdmol.GetAtomWithIdx(atom_idx).GetBonds()
            ]
        )
        neighbors: tuple[int, ...] = tuple(id_atom_map[b] for b in neighbors)
        # idx -> atom map num

        chiral_tag = atom.GetChiralTag()
        hybridization = atom.GetHybridization()
        # rad_elec = atom.GetNumRadicalElectrons()

        if len(neighbors) == 4:
            if chiral_tag in rd_tetrahedral:  #
                atom_stereo: AtomStereo = Tetrahedral(
                    (id_atom_map[atom_idx], *neighbors),
                    rd_tetrahedral[chiral_tag],
                )

                graph.set_atom_stereo(atom_stereo)

            elif hybridization == Chem.HybridizationType.SP3:
                atom_stereo = Tetrahedral((id_atom_map[atom_idx], *neighbors), None)
                graph.set_atom_stereo(atom_stereo)

        if atom.GetChiralTag() == Chem.ChiralType.CHI_SQUAREPLANAR:
            
            sp_order: tuple[int, int, int, int]
            if atom.GetUnsignedProp("_chiralPermutation") == 1:
                sp_order = (0, 1, 2, 3)
            elif atom.GetUnsignedProp("_chiralPermutation") == 2:
                sp_order = (0, 2, 1, 3)
            elif atom.GetUnsignedProp("_chiralPermutation") == 3:
                sp_order = (0, 1, 3, 2)
            else:
                raise RuntimeError("Unknown permutation for SquarePlanar")
            ordered_neighbors = tuple([neighbors[i] for i in sp_order])
            sp_atoms = (id_atom_map[atom_idx], *ordered_neighbors)
            assert len(sp_atoms) == 5
            atom_stereo = SquarePlanar(atoms=sp_atoms, parity=0)
            graph.set_atom_stereo(atom_stereo)

        if atom.GetChiralTag() == Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL:
            perm = atom.GetUnsignedProp("_chiralPermutation")

            # adapted from http://opensmiles.org/opensmiles.html
            atom_order_permutation_dict = {
                (0, 1, 2, 3, 4): 1,
                (0, 1, 3, 2, 4): 2,
                (0, 1, 2, 4, 3): 3,
                (0, 1, 4, 2, 3): 4,
                (0, 1, 3, 4, 2): 5,
                (0, 1, 4, 3, 2): 6,
                (0, 2, 3, 4, 1): 7,
                (0, 2, 4, 3, 1): 8,
                (1, 0, 2, 3, 4): 9,
                (1, 0, 3, 2, 4): 11,
                (1, 0, 2, 4, 3): 10,
                (1, 0, 4, 2, 3): 12,
                (1, 0, 3, 4, 2): 13,
                (1, 0, 4, 3, 2): 14,
                (2, 0, 1, 3, 4): 15,
                (2, 0, 1, 4, 3): 16,
                (3, 0, 1, 2, 4): 17,
                (3, 0, 2, 1, 4): 18,
                (2, 0, 4, 1, 3): 19,
                (2, 0, 3, 1, 4): 20,
            }

            permutation_atom_order_dict = {
                v: k for k, v in atom_order_permutation_dict.items()
            }

            tbp_order = permutation_atom_order_dict[perm]
            neigh_atoms = tuple([neighbors[i] for i in tbp_order])
            tbp_atoms = (id_atom_map[atom_idx], *neigh_atoms)
            assert len(tbp_atoms) == 6
            atom_stereo = TrigonalBipyramidal(tbp_atoms, 1)
            graph.set_atom_stereo(atom_stereo)

        if atom.GetChiralTag() == Chem.ChiralType.CHI_OCTAHEDRAL:
            perm = atom.GetUnsignedProp("_chiralPermutation")

            permutation_atom_order_dict = {
                1: (0, 5, 1, 2, 3, 4),
                2: (0, 5, 1, 4, 3, 2),
                3: (0, 4, 1, 2, 3, 5),
                16: (0, 4, 1, 5, 3, 2),
                6: (0, 3, 1, 2, 4, 5),
                18: (0, 3, 1, 5, 4, 2),
                19: (0, 2, 1, 3, 4, 5),
                24: (0, 2, 1, 5, 4, 3),
                25: (0, 1, 2, 3, 4, 5),
                30: (0, 1, 2, 5, 4, 3),
                4: (0, 5, 1, 2, 4, 3),
                14: (0, 5, 1, 3, 4, 2),
                5: (0, 4, 1, 2, 5, 3),
                15: (0, 4, 1, 3, 5, 2),
                7: (0, 3, 1, 2, 5, 4),
                17: (0, 3, 1, 4, 5, 2),
                20: (0, 2, 1, 3, 5, 4),
                23: (0, 2, 1, 4, 5, 3),
                26: (0, 1, 2, 3, 5, 4),
                29: (0, 1, 2, 4, 5, 3),
                10: (0, 5, 1, 4, 2, 3),
                8: (0, 5, 1, 3, 2, 4),
                11: (0, 4, 1, 5, 2, 3),
                9: (0, 4, 1, 3, 2, 5),
                13: (0, 3, 1, 5, 2, 4),
                12: (0, 3, 1, 4, 2, 5),
                22: (0, 2, 1, 5, 3, 4),
                21: (0, 2, 1, 4, 3, 5),
                28: (0, 1, 2, 5, 3, 4),
                27: (0, 1, 2, 4, 3, 5),
            }

            order = permutation_atom_order_dict[perm]
            neigh_atoms = tuple([neighbors[i] for i in order])
            oct_atoms = (id_atom_map[atom_idx], *neigh_atoms)
            assert len(oct_atoms) == 7
            atom_stereo = Octahedral(oct_atoms, 1)
            graph.set_atom_stereo(atom_stereo)

    for bond in (
        b
        for b in rdmol.GetBonds()
        if b.GetIsConjugated()
        or b.GetBondType() == Chem.rdchem.BondType.DOUBLE
        or b.GetStereo()
        in (Chem.BondStereo.STEREOATROPCW, Chem.BondStereo.STEREOATROPCCW)
    ):
        begin_end_idx: tuple[int, int] = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

        neighbors_begin: list[int] = [
            atom.GetIdx()
            for atom in rdmol.GetAtomWithIdx(begin_end_idx[0]).GetNeighbors()
            if atom.GetIdx() != begin_end_idx[1]
        ]

        neighbors_end: list[int] = [
            atom.GetIdx()
            for atom in rdmol.GetAtomWithIdx(begin_end_idx[1]).GetNeighbors()
            if atom.GetIdx() != begin_end_idx[0]
        ]

        if len({*neighbors_begin, *neighbors_end}) != 4:
            continue
        # TODO: how to deal with double bonds in strained structures?
        # cyclopropane ?

        if len(neighbors_begin) != 2 or len(neighbors_end) != 2:
            continue
        # TODO: how to deal with imines?

        elif bond.GetStereo() in (
            Chem.BondStereo.STEREOATROPCW,
            Chem.BondStereo.STEREOATROPCCW,
        ):
            raise NotImplementedError("Atropisomerism is not implemented yet. ")

        elif (
            bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
            and [a for a in bond.GetStereoAtoms()] != []
        ):
            if bond.GetStereo() == Chem.BondStereo.STEREONONE:
                bond_atoms_idx = (
                    (
                        *neighbors_begin,
                        begin_end_idx[0],
                        begin_end_idx[1],
                        *neighbors_end,
                    )
                )
                
                bond_atoms = tuple([id_atom_map[i] for i in bond_atoms_idx])
                assert len(bond_atoms) == 6
                stereo = PlanarBond(bond_atoms, None)
                invert = None
            else:
                if bond.GetStereo() == Chem.BondStereo.STEREOZ:
                    invert = False
                elif bond.GetStereo() == Chem.BondStereo.STEREOE:
                    invert = True
                else:
                    raise RuntimeError("Unknown Stereo")

                stereo_atoms = [a for a in bond.GetStereoAtoms()]

                if (
                    stereo_atoms[0] in neighbors_begin
                    and stereo_atoms[1] in neighbors_end
                ):
                    bond_atoms_idx = (
                        stereo_atoms[0],
                        *[n for n in neighbors_begin if n != stereo_atoms[0]],
                        begin_end_idx[0],
                        begin_end_idx[1],
                        stereo_atoms[1],
                        *[n for n in neighbors_end if n != stereo_atoms[1]],
                    )

                    bond_atoms = tuple([id_atom_map[a] for a in bond_atoms_idx])

                    # raise Exception(bond_atoms_idx)

                elif (
                    stereo_atoms[0] in neighbors_end
                    and stereo_atoms[1] in neighbors_begin
                ):
                    bond_atoms_idx = (
                        stereo_atoms[0],
                        *[n for n in neighbors_end if n != stereo_atoms[0]],
                        begin_end_idx[1],
                        begin_end_idx[0],
                        stereo_atoms[1],
                        *[n for n in neighbors_begin if n != stereo_atoms[1]],
                    )

                    bond_atoms = tuple([id_atom_map[a] for a in bond_atoms_idx])
                else:
                    raise RuntimeError("Stereo Atoms not neighbors")

                if invert is True:
                    inverted_atoms = tuple([bond_atoms[i] for i in (1, 0, 2, 3, 4, 5)])
                    assert len(inverted_atoms) == 6
                    stereo = PlanarBond(inverted_atoms, 0)
                elif invert is False:
                    assert len(bond_atoms) == 6
                    stereo = PlanarBond(bond_atoms, 0)

        elif bond.GetBondType() == Chem.rdchem.BondType.AROMATIC:
            ri = rdmol.GetRingInfo()
            rings: list[set[int]] = [set(ring) for ring in ri.AtomRings()]
            stereo_atoms = [
                neighbors_begin[0],
                neighbors_begin[1],
                begin_end_idx[0],
                begin_end_idx[1],
                neighbors_end[0],
                neighbors_end[1],
            ]

            common_ring_size_db1: None | int = None
            for ring in rings:
                if all(
                    a in ring
                    for a in (
                        neighbors_begin[0],
                        begin_end_idx[0],
                        begin_end_idx[1],
                        neighbors_end[0],
                    )
                ):
                    if common_ring_size_db1 is None or len(ring) < common_ring_size_db1:
                        common_ring_size_db1 = len(ring)
                if all(
                    a in ring
                    for a in (
                        neighbors_begin[1],
                        begin_end_idx[0],
                        begin_end_idx[1],
                        neighbors_end[1],
                    )
                ):
                    if common_ring_size_db1 is None or len(ring) < common_ring_size_db1:
                        common_ring_size_db1 = len(ring)

            common_ring_size_db2: None|int = None

            for ring in rings:
                if all(
                    a in ring
                    for a in (
                        neighbors_begin[1],
                        begin_end_idx[0],
                        begin_end_idx[1],
                        neighbors_end[0],
                    )
                ):
                    if common_ring_size_db2 is None or len(ring) < common_ring_size_db2:
                        common_ring_size_db2 = len(ring)
            for ring in rings:
                if all(
                    a in ring
                    for a in (
                        neighbors_begin[0],
                        begin_end_idx[0],
                        begin_end_idx[1],
                        neighbors_end[1],
                    )
                ):
                    if common_ring_size_db2 is None or len(ring) < common_ring_size_db2:
                        common_ring_size_db2 = len(ring)

            if common_ring_size_db1 is None and common_ring_size_db2 is None:
                pb_atoms = tuple([id_atom_map[a] for a in stereo_atoms])
                assert len(pb_atoms) == 6
                stereo = PlanarBond(pb_atoms, None)
            elif common_ring_size_db1:
                pb_atoms = tuple(tuple([id_atom_map[a] for a in stereo_atoms]))
                assert len(pb_atoms) == 6
                stereo = PlanarBond(pb_atoms, parity=0)
            elif common_ring_size_db2:
                pb_atoms = tuple(tuple([id_atom_map[stereo_atoms[i]] for i in (0, 1, 2, 3, 4, 5)]))
                assert len(pb_atoms) == 6
                stereo = PlanarBond(pb_atoms, parity=0,)
            elif common_ring_size_db1 is None or common_ring_size_db2 is None:
                raise RuntimeError(
                    "Aromatic Atoms not in ring, "
                    "please check the input structure"
                )
            elif common_ring_size_db2 < common_ring_size_db1:
                pb_atoms = tuple(tuple([id_atom_map[a] for a in stereo_atoms]))
                assert len(pb_atoms) == 6
                stereo = PlanarBond(pb_atoms, parity=0)
            elif common_ring_size_db1 < common_ring_size_db2:
                pb_atoms = tuple([id_atom_map[stereo_atoms[i]] for i in (0, 1, 2, 3, 4, 5)])
                assert len(pb_atoms) == 6
                stereo = PlanarBond(pb_atoms, parity=0)
            else:
                raise RuntimeError("Aromatic Atoms not in ring")

        else:
            stereo_atoms = [
                neighbors_begin[0],
                neighbors_begin[1],
                begin_end_idx[0],
                begin_end_idx[1],
                neighbors_end[0],
                neighbors_end[1],
            ]
            pb_atoms = tuple([id_atom_map[a] for a in stereo_atoms])
            assert len(pb_atoms) == 6
            stereo = PlanarBond(pb_atoms, None)

        # raise Exception(begin_end_idx, )
        bond_atoms = [id_atom_map[i] for i in begin_end_idx]

        graph.set_bond_stereo(stereo) # type: ignore

    return graph
