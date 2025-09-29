# pyright: standard
# typing with rdkit is not fully supported
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import TYPE_CHECKING

import rdkit.Chem as Chem  # type: ignore

from stereomolgraph import AtomId, MolGraph, StereoMolGraph
from stereomolgraph.stereodescriptors import (
    AtomStereo,
    AtropBond,
    Octahedral,
    PlanarBond,
    SquarePlanar,
    Tetrahedral,
    TrigonalBipyramidal,
)

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import ClassVar, Literal


def mol_graph_from_rdmol(
    cls: type[MolGraph], rdmol: Chem.Mol, use_atom_map_number: bool = False
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
        id_atom_map = {
            atom.GetIdx(): atom.GetAtomMapNum() for atom in rdmol.GetAtoms()
        }
    else:
        id_atom_map = {
            atom.GetIdx(): atom.GetIdx() for atom in rdmol.GetAtoms()
        }

    for atom in rdmol.GetAtoms():
        graph.add_atom(id_atom_map[atom.GetIdx()], atom.GetSymbol())

    for bond in rdmol.GetBonds():
        graph.add_bond(
            id_atom_map[bond.GetBeginAtomIdx()],
            id_atom_map[bond.GetEndAtomIdx()],
        )
    return graph


@dataclass
class RDMol2StereoMolGraph:
    resonance: bool = True
    stereo_complete: bool = False
    use_atom_map_number: bool = False

    # aromatic_cis: bool = True # TODO:
    # "aromatic bonds are always planar and cis if no stereochemistry is defined"

    def __call__(self, rdmol: Chem.Mol) -> StereoMolGraph:
        return self.smg_from_rdmol(rdmol)
        
    def smg_from_rdmol(self, rdmol: Chem.Mol) -> StereoMolGraph:
        
        rdmol = Chem.AddHs(rdmol, explicitOnly=False)
        graph = StereoMolGraph()

        id_atom_map: dict[int, AtomId]

        if self.use_atom_map_number is True:
            if any(atom.GetAtomMapNum() == 0 for atom in rdmol.GetAtoms()):
                raise ValueError("AtomMapNumber has to  be set on all atoms")
            id_atom_map = {
                atom.GetIdx(): atom.GetAtomMapNum()
                for atom in rdmol.GetAtoms()
            }
        else:
            id_atom_map = {
                atom.GetIdx(): atom.GetIdx() for atom in rdmol.GetAtoms()
            }

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
                id_atom_map[n.GetIdx()]
                for n in rdmol.GetAtomWithIdx(atom_idx).GetNeighbors()
            )

            chiral_tag = atom.GetChiralTag()
            hybridization = atom.GetHybridization()

            atom_stereo: None | AtomStereo = None

            if chiral_tag in self._rd_tetrahedral:
                stereo_atoms = (id_atom_map[atom_idx], *neighbors)
                assert len(stereo_atoms) == 5

                atom_stereo = Tetrahedral(
                    stereo_atoms,
                    self._rd_tetrahedral[chiral_tag],
                )

            elif chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL or (
                hybridization == Chem.HybridizationType.SP3
                and len(neighbors) == 4
            ):
                short_stereo_atoms = (id_atom_map[atom_idx], *neighbors)
                stereo_atoms = tuple(
                    short_stereo_atoms[i]
                    if i < len(short_stereo_atoms)
                    else None
                    for i in range(5)
                ) # extends with "None" if less than 4 neighbors
                assert len(stereo_atoms) == 5

                if not self.stereo_complete:
                    atom_stereo = Tetrahedral(stereo_atoms, None)
                elif self.stereo_complete:
                    atom_stereo = Tetrahedral(stereo_atoms, parity=1)
                else:
                    raise RuntimeError("This should never happen")

            elif chiral_tag == Chem.ChiralType.CHI_SQUAREPLANAR:
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

            elif chiral_tag == Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL:
                perm = atom.GetUnsignedProp("_chiralPermutation")
                tbp_order = self._tbp_atom_order_permutation_dict[perm]
                neigh_atoms = tuple([neighbors[i] for i in tbp_order])
                tbp_atoms = (id_atom_map[atom_idx], *neigh_atoms)
                assert len(tbp_atoms) == 6
                atom_stereo = TrigonalBipyramidal(tbp_atoms, 1)

            elif chiral_tag == Chem.ChiralType.CHI_OCTAHEDRAL:
                perm = atom.GetUnsignedProp("_chiralPermutation")
                order = self._oct_atom_order_permutation_dict[perm]
                neigh_atoms = tuple([neighbors[i] for i in order])
                oct_atoms = (id_atom_map[atom_idx], *neigh_atoms)
                assert len(oct_atoms) == 7
                atom_stereo = Octahedral(oct_atoms, 1)

            else:
                continue
            
            graph.set_atom_stereo(atom_stereo)

        for bond in (
            b
            for b in rdmol.GetBonds()
            if b.GetIsAromatic()  # TODO: check if include conjugated bonds?
            or b.GetBondType() == Chem.rdchem.BondType.DOUBLE
            or b.GetStereo() in self._rd_atrop.keys()
        ):
            begin_idx: int
            end_idx: int
            begin_idx, end_idx = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

            rd_bond_stereo = bond.GetStereo()

            neighbors_begin: list[int] = [
                atom.GetIdx()
                for atom in rdmol.GetAtomWithIdx(begin_idx).GetNeighbors()
                if atom.GetIdx() != end_idx
            ]

            neighbors_end: list[int] = [
                atom.GetIdx()
                for atom in rdmol.GetAtomWithIdx(end_idx).GetNeighbors()
                if atom.GetIdx() != begin_idx
            ]

            if (
                rd_bond_stereo == Chem.BondStereo.STEREOATROPCW
                or rd_bond_stereo == Chem.BondStereo.STEREOATROPCCW
            ):
                atrop_parity = self._rd_atrop[rd_bond_stereo]
                stereo_atoms = [a for a in bond.GetStereoAtoms()]

                if (
                    stereo_atoms[0] in neighbors_begin
                    and stereo_atoms[1] in neighbors_end
                ):
                    bond_atoms_idx = (
                        stereo_atoms[0],
                        *[n for n in neighbors_begin if n != stereo_atoms[0]],
                        begin_idx,
                        end_idx,
                        stereo_atoms[1],
                        *[n for n in neighbors_end if n != stereo_atoms[1]],
                    )

                    bond_atoms = tuple(
                        [id_atom_map[a] for a in bond_atoms_idx]
                    )

                elif (
                    stereo_atoms[0] in neighbors_end
                    and stereo_atoms[1] in neighbors_begin
                ):
                    bond_atoms_idx = (
                        stereo_atoms[0],
                        *[n for n in neighbors_end if n != stereo_atoms[0]],
                        begin_idx,
                        end_idx,
                        stereo_atoms[1],
                        *[n for n in neighbors_begin if n != stereo_atoms[1]],
                    )

                    bond_atoms = tuple(
                        [id_atom_map[a] for a in bond_atoms_idx]
                    )
                else:
                    raise RuntimeError("Stereo Atoms not neighbors")

                assert len(bond_atoms) == 6
                stereo = AtropBond(bond_atoms, atrop_parity)


            elif (bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
                  and rd_bond_stereo == Chem.BondStereo.STEREONONE
                  and len(neighbors_begin) == 2 == len(neighbors_end)
                ):

                rd_bond_stereo = bond.GetStereo()
                invert: None | bool = None

                bond_atoms_idx = (
                        *neighbors_begin,
                        begin_idx,
                        end_idx,
                        *neighbors_end,
                    )
                
                bond_atoms = tuple(
                [id_atom_map.get(a) for a in bond_atoms_idx]
                )
                assert len(bond_atoms) == 6, bond_atoms
                stereo = PlanarBond(bond_atoms, None)
            
            elif rd_bond_stereo in (Chem.BondStereo.STEREOZ,
                                    Chem.BondStereo.STEREOE):
                
                invert = {Chem.BondStereo.STEREOZ: False,
                          Chem.BondStereo.STEREOE: True}[rd_bond_stereo]

                begin_stereo_atom: int
                end_stereo_atom: int
                begin_stereo_atom, end_stereo_atom = [
                    a for a in bond.GetStereoAtoms()
                ]
                begin_non_stereo_nbr = (
                    None
                    if len(neighbors_begin) == 1
                    else [
                        a
                        for a in neighbors_begin
                        if a != begin_stereo_atom
                    ][0]
                )
                end_non_stereo_nbr = (
                    None
                    if len(neighbors_end) == 1
                    else [
                        a for a in neighbors_end if a != end_stereo_atom
                    ][0]
                )

                bond_atoms_idx = (
                    begin_stereo_atom,
                    begin_non_stereo_nbr,
                    begin_idx,
                    end_idx,
                    end_stereo_atom,
                    end_non_stereo_nbr,
                )
                assert len(bond_atoms_idx) == 6, bond_atoms_idx

                assert len(bond_atoms_idx) == 6, bond_atoms_idx

                bond_atoms = tuple(
                [id_atom_map.get(a) for a in bond_atoms_idx]
                )

                if invert:
                    bond_atoms = tuple(
                    [bond_atoms[i] for i in (1, 0, 2, 3, 4, 5)]
                    )

                assert len(bond_atoms) == 6, bond_atoms
                stereo = PlanarBond(bond_atoms, 0)

            elif bond.GetBondType() == Chem.rdchem.BondType.AROMATIC or (
                bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
                and rd_bond_stereo == Chem.BondStereo.STEREONONE
                and (len(neighbors_begin) == 1 or len(neighbors_end) == 1)
            ):
                print("Aromatic Bond Stereo")
                ri = rdmol.GetRingInfo()
                rings: list[set[int]] = [set(ring) for ring in ri.AtomRings()]

                neighbors_begin0 = neighbors_begin[0]
                neighbors_begin1 = (
                    None if len(neighbors_begin) == 1 else neighbors_begin[1]
                )
                neighbors_end0 = neighbors_end[0]
                neighbors_end1 = (
                    None if len(neighbors_end) == 1 else neighbors_end[1]
                )
                stereo_atoms = [
                    neighbors_begin0,
                    neighbors_begin1,
                    begin_idx,
                    end_idx,
                    neighbors_end0,
                    neighbors_end1,
                ]

                common_ring_size_db1: None | int = None
                for ring in rings:
                    if all(
                        a in ring
                        for a in (
                            neighbors_begin0,
                            begin_idx,
                            end_idx,
                            neighbors_end0,
                        )
                    ):
                        if (
                            common_ring_size_db1 is None
                            or len(ring) < common_ring_size_db1
                        ):
                            common_ring_size_db1 = len(ring)
                    if all(
                        a in ring
                        for a in (
                            neighbors_begin1,
                            begin_idx,
                            end_idx,
                            neighbors_end1,
                        )
                    ):
                        if (
                            common_ring_size_db1 is None
                            or len(ring) < common_ring_size_db1
                        ):
                            common_ring_size_db1 = len(ring)

                common_ring_size_db2: None | int = None

                for ring in rings:
                    if all(
                        a in ring
                        for a in (
                            neighbors_begin1,
                            begin_idx,
                            end_idx,
                            neighbors_end0,
                        )
                    ):
                        if (
                            common_ring_size_db2 is None
                            or len(ring) < common_ring_size_db2
                        ):
                            common_ring_size_db2 = len(ring)
                for ring in rings:
                    if all(
                        a in ring
                        for a in (
                            neighbors_begin0,
                            begin_idx,
                            end_idx,
                            neighbors_end1,
                        )
                    ):
                        if (
                            common_ring_size_db2 is None
                            or len(ring) < common_ring_size_db2
                        ):
                            common_ring_size_db2 = len(ring)

                if (
                    common_ring_size_db1 is None
                    and common_ring_size_db2 is None
                ):
                    pb_atoms = tuple(
                        [id_atom_map.get(a) for a in stereo_atoms]
                    )
                    assert len(pb_atoms) == 6
                    stereo = PlanarBond(pb_atoms, None)
                elif common_ring_size_db1:
                    pb_atoms = tuple(
                        tuple([id_atom_map.get(a) for a in stereo_atoms])
                    )
                    assert len(pb_atoms) == 6
                    stereo = PlanarBond(pb_atoms, parity=0)
                elif common_ring_size_db2:
                    pb_atoms = tuple(
                        tuple(
                            [
                                id_atom_map.get(stereo_atoms[i])
                                for i in (0, 1, 2, 3, 4, 5)
                            ]
                        )
                    )
                    assert len(pb_atoms) == 6
                    stereo = PlanarBond(
                        pb_atoms,
                        parity=0,
                    )
                elif (
                    common_ring_size_db1 is None
                    or common_ring_size_db2 is None
                ):
                    raise RuntimeError(
                        "Aromatic Atoms not in ring, "
                        "please check the input structure"
                    )
                elif common_ring_size_db2 < common_ring_size_db1:
                    pb_atoms = tuple(
                        tuple([id_atom_map.get(a) for a in stereo_atoms])
                    )
                    assert len(pb_atoms) == 6
                    stereo = PlanarBond(pb_atoms, parity=0)
                elif common_ring_size_db1 < common_ring_size_db2:
                    pb_atoms = tuple(
                        [
                            id_atom_map.get(stereo_atoms[i])
                            for i in (0, 1, 2, 3, 4, 5)
                        ]
                    )
                    assert len(pb_atoms) == 6
                    stereo = PlanarBond(pb_atoms, parity=0)
                else:
                    raise RuntimeError(
                        "Aromatic Atoms not in ring, "
                        "please check the input structure"
                    )

            else:
                continue
                stereo_atoms = [
                    neighbors_begin[0],
                    neighbors_begin[1],
                    begin_idx,
                    end_idx,
                    neighbors_end[0],
                    neighbors_end[1],
                ]
                pb_atoms = tuple([id_atom_map[a] for a in stereo_atoms])
                assert len(pb_atoms) == 6
                stereo = PlanarBond(pb_atoms, None)

            graph.set_bond_stereo(stereo)

        return graph

    _rd_tetrahedral: ClassVar[
        Mapping[Chem.rdchem.ChiralType, None | Literal[-1, 1]]
    ] = MappingProxyType(
        {
            Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW: -1,
            Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW: 1,
            Chem.rdchem.ChiralType.CHI_TETRAHEDRAL: None,
        }
    )

    _rd_atrop: ClassVar[Mapping[Chem.BondStereo, None | Literal[-1, 1]]] = (
        MappingProxyType(
            {
                Chem.BondStereo.STEREOATROPCW: 1,
                Chem.BondStereo.STEREOATROPCCW: -1,
            }
        )
    )

    _tbp_atom_order_permutation_dict = MappingProxyType(
        {
            1: (0, 1, 2, 3, 4),
            2: (0, 1, 3, 2, 4),
            3: (0, 1, 2, 4, 3),
            4: (0, 1, 4, 2, 3),
            5: (0, 1, 3, 4, 2),
            6: (0, 1, 4, 3, 2),
            7: (0, 2, 3, 4, 1),
            8: (0, 2, 4, 3, 1),
            9: (1, 0, 2, 3, 4),
            11: (1, 0, 3, 2, 4),
            10: (1, 0, 2, 4, 3),
            12: (1, 0, 4, 2, 3),
            13: (1, 0, 3, 4, 2),
            14: (1, 0, 4, 3, 2),
            15: (2, 0, 1, 3, 4),
            16: (2, 0, 1, 4, 3),
            17: (3, 0, 1, 2, 4),
            18: (3, 0, 2, 1, 4),
            19: (2, 0, 4, 1, 3),
            20: (2, 0, 3, 1, 4),
        }
    )
    "adapted from http://opensmiles.org/opensmiles.html"

    _oct_atom_order_permutation_dict = MappingProxyType(
        {
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
    )
