# pyright: standard
from __future__ import annotations

from typing import NamedTuple

import rdkit.Chem.rdDepictor as rdDepictor
from rdkit import Chem
from rdkit.Chem import Draw  # type: ignore

from stereomolgraph import (
    CondensedReactionGraph,
    MolGraph,
    StereoCondensedReactionGraph,
    StereoMolGraph,
)
from stereomolgraph.graph2rdmol import (
    RDKitAtomId,
    condensed_reaction_graph_to_rdmol,
    mol_graph_to_rdmol,
    stereo_condensed_reaction_graph_to_rdmol,
    stereo_mol_graph_to_rdmol,
)
from stereomolgraph.graphs.mg import AtomId
from stereomolgraph.graphs.scrg import Change
from stereomolgraph.stereodescriptors import PlanarBond, Tetrahedral


def default_repr_svg(self: MolGraph) -> str:
    return View2D().svg(self)


def default_view_molgraph(self: MolGraph) -> None:
    View2D()(self)


MolGraph._ipython_display_ = default_view_molgraph
MolGraph._repr_svg_ = default_repr_svg


class _HighlightTuple(NamedTuple):
    atoms_to_highlight: list
    highlight_atom_colors: dict
    bonds_to_highlight: list
    highlight_bond_colors: dict


class View2D(NamedTuple):
    """A class to visualize a MolGraph in 2D using RDKit's MolDraw2DSVG.
    This class can be used in IPython environments to display the graph as an
    SVG image.

    :param height: Height of the SVG image in pixels
    :param width: Width of the SVG image in pixels
    :param show_atom_numbers: Whether to show atom numbers in the visualization
    :param show_h: Whether to show hydrogen atoms in the visualization
    :param generate_bond_orders: Whether to generate bond orders for the
        visualization using
        :func:`~stereomolgraph.algorithms.bond_orders.connectivity2bond_orders`
    :param dummy_atoms: Whether to include dummy atoms in the visualization
    """

    height: int = 300
    width: int = 300
    show_atom_numbers: bool = True
    show_h: bool = True
    generate_bond_orders: bool = False
    dummy_atoms: bool = False
    color_planar_bond_changes: bool = False

    def _to_mol(
        self,
        graph: (
            MolGraph
            | CondensedReactionGraph
            | StereoMolGraph
            | StereoCondensedReactionGraph
        ),
    ) -> tuple[Chem.Mol, _HighlightTuple]:
        match graph:
            case StereoCondensedReactionGraph():
                to_rdmol = stereo_condensed_reaction_graph_to_rdmol
            case CondensedReactionGraph():
                to_rdmol = condensed_reaction_graph_to_rdmol
            case StereoMolGraph():
                to_rdmol = stereo_mol_graph_to_rdmol
            case _:
                to_rdmol = mol_graph_to_rdmol
        mol, idx_map_num_dict = to_rdmol(
            graph,
            generate_bond_orders=self.generate_bond_orders,
        )
        map_num_idx_dict: dict[AtomId, RDKitAtomId] = {
            v: k for k, v in idx_map_num_dict.items()
        }

        if not self.generate_bond_orders:
            for bond in mol.GetBonds():
                bond.SetBondType(Chem.BondType.SINGLE)
                bond.SetIsAromatic(False)

        if self.show_atom_numbers:
            for atom in mol.GetAtoms():
                atom.SetProp("atomNote", str(atom.GetAtomMapNum()))
                atom.SetAtomMapNum(0)
        else:
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)

        atoms_to_highlight = []
        highlight_atom_colors = {}

        bonds_to_highlight = []
        highlight_bond_colors = {}
        formed_bonds = set()
        broken_bonds = set()

        if self.dummy_atoms is False:
            dummy_atoms = [
                atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "*"
            ]
            dummy_atoms.sort(reverse=True)
            for atom in dummy_atoms:
                mol.RemoveAtom(atom)

        if not self.show_h:
            mol = Chem.RemoveHs(mol, implicitOnly=False, sanitize=False)

        rdDepictor.Compute2DCoords(
            mol,
            clearConfs=True,
            sampleSeed=42,
            nSample=100,
            permuteDeg4Nodes=True,
            useRingTemplates=True,
        )
        rdDepictor.StraightenDepiction(mol)

        if isinstance(graph, StereoMolGraph) and not self.generate_bond_orders:
            # RDKit only wedges tetrahedral stereo on single bonds; the
            # default no-bond-order path leaves them as UNSPECIFIED.
            for atom_stereo in graph.atom_stereo.values():
                if (
                    not isinstance(atom_stereo, Tetrahedral)
                    or atom_stereo.parity is None
                ):
                    continue

                atom_idx = map_num_idx_dict[atom_stereo.atoms[0]]
                for bond in mol.GetAtomWithIdx(atom_idx).GetBonds():
                    if bond.GetBondType() == Chem.BondType.UNSPECIFIED:
                        bond.SetBondType(Chem.BondType.SINGLE)

        if isinstance(graph, StereoMolGraph) and not self.generate_bond_orders:
            for db in graph.bond_stereo.values():
                if (
                    isinstance(db, PlanarBond)
                    and isinstance(db.atoms[2], int)
                    and isinstance(db.atoms[3], int)
                ):
                    a1 = map_num_idx_dict[db.atoms[2]]
                    a2 = map_num_idx_dict[db.atoms[3]]
                    rd_bond = mol.GetBondBetweenAtoms(a1, a2)
                    rd_bond.SetBondType(Chem.BondType.AROMATIC)
                    rd_bond.SetIsAromatic(False)

        if isinstance(graph, CondensedReactionGraph):
            formed_bonds = graph.get_formed_bonds()
            broken_bonds = graph.get_broken_bonds()

            for bond in formed_bonds:
                atoms_idx = [map_num_idx_dict[a] for a in bond]
                bond_idx = mol.GetBondBetweenAtoms(*atoms_idx).GetIdx()
                mol.GetBondWithIdx(bond_idx).SetBondType(Chem.rdchem.BondType.HYDROGEN)
                bonds_to_highlight.append(bond_idx)
                highlight_bond_colors[bond_idx] = (0, 0, 1)  # blue

            for bond in broken_bonds:
                atoms_idx = [map_num_idx_dict[a] for a in bond]
                bond_idx = mol.GetBondBetweenAtoms(*atoms_idx).GetIdx()
                mol.GetBondWithIdx(bond_idx).SetBondType(Chem.rdchem.BondType.HYDROGEN)
                bonds_to_highlight.append(bond_idx)
                highlight_bond_colors[bond_idx] = (1, 0, 0)  # red

        if isinstance(graph, StereoCondensedReactionGraph):
            for bond, change_dict in graph.bond_stereo_changes.items():
                atoms_idx = [map_num_idx_dict[a] for a in bond]
                rd_bond = mol.GetBondBetweenAtoms(*atoms_idx)
                bond_idx = rd_bond.GetIdx()

                has_planar_change = any(
                    isinstance(stereo, PlanarBond)
                    for stereo in change_dict.values()
                    if stereo is not None
                )
                if (
                    has_planar_change
                    and not self.generate_bond_orders
                    and bond not in formed_bonds
                    and bond not in broken_bonds
                ):
                    rd_bond.SetBondType(Chem.BondType.AROMATIC)
                    rd_bond.SetIsAromatic(False)

                if change_dict[Change.FLEETING] is not None:
                    bonds_to_highlight.append(bond_idx)
                    highlight_bond_colors[bond_idx] = (1, 0, 1)  # magenta
                elif (
                    change_dict[Change.FORMED] is not None
                    and change_dict[Change.BROKEN] is not None
                ):
                    bonds_to_highlight.append(bond_idx)
                    highlight_bond_colors[bond_idx] = (1, 0, 1)  # magenta
                elif change_dict[Change.FORMED] is not None:
                    bonds_to_highlight.append(bond_idx)
                    highlight_bond_colors[bond_idx] = (0, 0, 1)  # blue
                elif change_dict[Change.BROKEN] is not None:
                    bonds_to_highlight.append(bond_idx)
                    highlight_bond_colors[bond_idx] = (1, 0, 0)  # red

            # make dummy atoms and their bonds grey
        if self.dummy_atoms is True:
            atom_colors = {}
            grey = (0.7, 0.7, 0.7)
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == "*":
                    atom.SetProp("atomLabel", "")
                    atom_colors[atom.GetIdx()] = grey
                    for bond in atom.GetBonds():
                        bonds_to_highlight.append(bond.GetIdx())
                        highlight_bond_colors[bond.GetIdx()] = grey
        ht = _HighlightTuple(
            atoms_to_highlight=atoms_to_highlight,
            highlight_atom_colors=highlight_atom_colors,
            bonds_to_highlight=bonds_to_highlight,
            highlight_bond_colors=highlight_bond_colors,
        )
        return mol, ht

    def svg(
        self,
        graph: (
            MolGraph
            | CondensedReactionGraph
            | StereoMolGraph
            | StereoCondensedReactionGraph
        ),
    ) -> str:
        mol, ht = self._to_mol(graph)

        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(self.width, self.height)

        drawer.drawOptions().useBWAtomPalette()
        drawer.drawOptions().continuousHighlight = False
        drawer.drawOptions().highlightBondWidthMultiplier = 12
        drawer.drawOptions().fillHighlights = False
        drawer.drawOptions().includeRadicals = False

        drawer.DrawMolecule(
            mol,
            highlightAtoms=ht.atoms_to_highlight,
            highlightAtomColors=ht.highlight_atom_colors,
            highlightBonds=ht.bonds_to_highlight,
            highlightBondColors=ht.highlight_bond_colors,
        )

        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg

    def __call__(
        self,
        graph: (
            MolGraph
            | CondensedReactionGraph
            | StereoMolGraph
            | StereoCondensedReactionGraph
        ),
    ):
        # imported here, so that this module does not depend on IPython
        from IPython.display import SVG

        svg = self.svg(graph)
        display(SVG(svg.replace("svg:", "")))  # type: ignore # noqa
