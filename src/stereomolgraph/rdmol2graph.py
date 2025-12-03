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
    BondStereo
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
    """Convert an RDKit :class:`rdkit.Chem.Mol` to a :class:`StereoMolGraph`.

    All aromatic bonds are considered to be cis. Double bonds in rings are
    assumed to be cis for rings of size <= 7.

    :param stereo_complete: If ``True``, attempt to infer complete stereo
        parities when possible. Defaults to ``False``.
    :param use_atom_map_number: Use RDKit atom map numbers as atom identifiers
        instead of RDKit atom indices. Defaults to ``False``.
    :param lone_pair_stereo: Include stereochemistry that depends on lone
        pairs (if present). Defaults to ``True``.
    :param resonance: Enumerate resonance structures and merge bond stereo
        information from them. Defaults to ``True``.
    """
    stereo_complete: bool = False
    use_atom_map_number: bool = False
    lone_pair_stereo: bool = True
    resonance: bool = True
    _max_resonance_structures: int = 100
    _min_trans_ring_size: int = 7

    def __call__(self, rdmol: Chem.Mol) -> StereoMolGraph:
        smg = self.smg_from_rdmol(rdmol)

        if self.resonance is False:
            return smg
        
        elif self.resonance is True:
            enumerator = (res_mol for res_mol
                          in Chem.ResonanceMolSupplier(
                              rdmol, Chem.KEKULE_ALL,
                              self._max_resonance_structures)
                         if res_mol is not None)

            for res_mol in enumerator:
                res_smg = self.smg_from_rdmol(res_mol)
                for bond, bond_stereo in res_smg.bond_stereo.items():
                    if bond not in smg.bond_stereo:
                        smg.set_bond_stereo(bond_stereo)
        return smg

    def smg_from_rdmol(self, rdmol: Chem.Mol) -> StereoMolGraph:
        if rdmol is None:
            raise ValueError("rdmol is None")
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
                if len(stereo_atoms) == 5:
                    atom_stereo = Tetrahedral(
                        stereo_atoms,
                        self._rd_tetrahedral[chiral_tag],
                    )
                elif len(stereo_atoms) == 4:
                    atom_stereo = Tetrahedral(
                        (*stereo_atoms, None),
                        self._rd_tetrahedral[chiral_tag],
                    )
                else:
                    raise RuntimeError(
                        "Tetrahedral stereo must have 3 or 4 neighbors"
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
                )  # extends with "None" if less than 4 neighbors
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
            
            if not self.lone_pair_stereo and None in atom_stereo.atoms:
                continue
            graph.set_atom_stereo(atom_stereo)

        for bond in (
            b
            for b in rdmol.GetBonds()
            if (rd_bond_stereo := b.GetStereo()) in self._rd_atrop.keys()
            or b.GetIsAromatic()
            or b.GetBondType() == Chem.rdchem.BondType.DOUBLE
        ):
            bond_stereo: BondStereo

            begin_idx: int = bond.GetBeginAtomIdx()
            end_idx: int = bond.GetEndAtomIdx()

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
                bond_stereo = AtropBond(bond_atoms, atrop_parity)

            elif (
                bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
                and rd_bond_stereo == Chem.BondStereo.STEREONONE
                and len(neighbors_begin) == 2 == len(neighbors_end)
            ):
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
                assert len(bond_atoms) == 6

                bond_stereo = PlanarBond(bond_atoms, None)


            elif rd_bond_stereo in (
                Chem.BondStereo.STEREOZ,
                Chem.BondStereo.STEREOE,
            ):
                invert = {
                    Chem.BondStereo.STEREOZ: False,
                    Chem.BondStereo.STEREOE: True,
                }[rd_bond_stereo]

                begin_stereo_atom: int
                end_stereo_atom: int
                begin_stereo_atom, end_stereo_atom = [
                    a for a in bond.GetStereoAtoms()
                ]
                begin_non_stereo_nbr = (
                    None
                    if len(neighbors_begin) == 1
                    else [
                        a for a in neighbors_begin if a != begin_stereo_atom
                    ][0]
                )
                end_non_stereo_nbr = (
                    None
                    if len(neighbors_end) == 1
                    else [a for a in neighbors_end if a != end_stereo_atom][0]
                )

                bond_atoms_idx = (
                    begin_stereo_atom,
                    begin_non_stereo_nbr,
                    begin_idx,
                    end_idx,
                    end_stereo_atom,
                    end_non_stereo_nbr,
                )
                assert len(bond_atoms_idx) == 6

                assert len(bond_atoms_idx) == 6

                bond_atoms = tuple(
                    [id_atom_map.get(a) for a in bond_atoms_idx]
                )
 
                if invert:
                    bond_atoms = tuple(
                        [bond_atoms[i] for i in (1, 0, 2, 3, 4, 5)]
                        )
                assert len(bond_atoms) == 6
                bond_stereo = PlanarBond(bond_atoms, 0)
                
            elif bond.GetBondType() == Chem.rdchem.BondType.AROMATIC or (
                bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
                and rd_bond_stereo == Chem.BondStereo.STEREONONE
                #and (len(neighbors_begin) == 1 or len(neighbors_end) == 1)
            ):
                # Find rings with bond begin_idx-end_idx, sort by aromatic first then size
                rings = [
                    (
                        len(r),
                        all(
                            rdmol.GetBondBetweenAtoms(
                                list(r)[i], list(r)[(i + 1) % len(r)]
                            ).GetIsAromatic()
                            for i in range(len(r))
                        ),
                        r,
                    )
                    for r in rdmol.GetRingInfo().AtomRings()
                    if begin_idx in r
                    and end_idx in r
                    and (
                        lambda rl: abs(rl.index(begin_idx) - rl.index(end_idx))
                        in [1, len(rl) - 1]
                    )(list(r))
                ]
                rings.sort(key=lambda x: (not x[1], x[0]))

                # Get ordered atoms from first ring (assumed cis)

                ring = rings[0][2] if rings else []
                n_begin = [
                    n.GetIdx()
                    for n in rdmol.GetAtomWithIdx(begin_idx).GetNeighbors()
                    if n.GetIdx() != end_idx
                ]
                n_end = [
                    n.GetIdx()
                    for n in rdmol.GetAtomWithIdx(end_idx).GetNeighbors()
                    if n.GetIdx() != begin_idx
                ]

                bond_atoms_idx = [
                    next(
                        (n for n in n_begin if n in ring), None
                    ),  # atom1: in-ring neighbor of begin_idx
                    next(
                        (n for n in n_begin if n not in ring), None
                    ),  # atom2: out-of-ring neighbor of begin_idx
                    begin_idx,  # atom3
                    end_idx,  # atom4
                    next(
                        (n for n in n_end if n in ring), None
                    ),  # atom5: in-ring neighbor of end_idx
                    next(
                        (n for n in n_end if n not in ring), None
                    ),  # atom6: out-of-ring neighbor of end_idx
                ]
                bond_atoms = tuple(
                    [id_atom_map.get(a) for a in bond_atoms_idx]
                )

                assert len(bond_atoms) == 6
                if rings and rings[0][1] is False and rings[0][0] > self._min_trans_ring_size:
                    bond_stereo = PlanarBond(bond_atoms, None)
                else:
                    bond_stereo = PlanarBond(bond_atoms, 0)

            else:
                continue
        
            if not self.lone_pair_stereo and None in bond_stereo.atoms:
                continue
            graph.set_bond_stereo(bond_stereo)

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
