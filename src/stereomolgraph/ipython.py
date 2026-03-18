# pyright: standard
from __future__ import annotations

from typing import NamedTuple

from rdkit import Chem  # type: ignore
from rdkit.Chem import Draw  # type: ignore

from stereomolgraph import (
    CondensedReactionGraph,
    MolGraph,
    StereoCondensedReactionGraph,
    StereoMolGraph,
)
from stereomolgraph.graphs.scrg import Change
from stereomolgraph.stereodescriptors import PlanarBond


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
    dummy_atoms: bool = True
    color_planar_bond_changes: bool = True

    def _to_mol(
        self,
        graph: (
            MolGraph
            | CondensedReactionGraph
            | StereoMolGraph
            | StereoCondensedReactionGraph
        ),
    ) -> tuple[Chem.Mol, _HighlightTuple]:
        mol, idx_map_num_dict = graph._to_rdmol(
            generate_bond_orders=self.generate_bond_orders
        )
        map_num_idx_dict = {v: k for k, v in idx_map_num_dict.items()}

        if not self.generate_bond_orders:
            for bond in mol.GetBonds():
                bond.SetBondType(Chem.BondType.SINGLE)

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

        if self.dummy_atoms is False:
            dummy_atoms = [
                atom.GetIdx()
                for atom in mol.GetAtoms()
                if atom.GetSymbol() == "*"
            ]
            dummy_atoms.sort(reverse=True)
            for atom in dummy_atoms:
                mol.RemoveAtom(atom)

        if not self.show_h:
            mol = Chem.RemoveHs(mol, implicitOnly=False, sanitize=False)

        if isinstance(graph, StereoMolGraph) and not self.generate_bond_orders:
            for db in graph.bond_stereo.values():
                if (isinstance(db, PlanarBond)
                    and isinstance(db.atoms[2], int)
                    and isinstance(db.atoms[3], int)):
                    a1 = map_num_idx_dict[db.atoms[2]]
                    a2 = map_num_idx_dict[db.atoms[3]]
                    rd_bond = mol.GetBondBetweenAtoms(a1, a2)
                    rd_bond.SetBondType(Chem.BondType.AROMATIC)

        if isinstance(graph, CondensedReactionGraph):
            for bond in graph.get_formed_bonds():
                atoms_idx = [map_num_idx_dict[a] for a in bond]
                bond_idx = mol.GetBondBetweenAtoms(*atoms_idx).GetIdx()
                bonds_to_highlight.append(bond_idx)
                mol.GetBondWithIdx(bond_idx).SetBondType(
                    Chem.rdchem.BondType.HYDROGEN
                )
                highlight_bond_colors[bond_idx] = (0, 0, 1)  # blue

            for bond in graph.get_broken_bonds():
                atoms_idx = [map_num_idx_dict[a] for a in bond]
                bond_idx = mol.GetBondBetweenAtoms(*atoms_idx).GetIdx()
                bonds_to_highlight.append(bond_idx)
                mol.GetBondWithIdx(bond_idx).SetBondType(
                    Chem.rdchem.BondType.HYDROGEN
                )
                highlight_bond_colors[bond_idx] = (1, 0, 0)  # red

        if self.color_planar_bond_changes and isinstance(
            graph, StereoCondensedReactionGraph
        ):
            for bond, change_dict in graph.bond_stereo_changes.items():
                if (
                    change_dict[Change.FORMED]
                    and not change_dict[Change.BROKEN]
                ):
                    atoms_idx = [map_num_idx_dict[a] for a in bond]
                    bond_idx = mol.GetBondBetweenAtoms(*atoms_idx).GetIdx()
                    bonds_to_highlight.append(bond_idx)
                    mol.GetBondWithIdx(bond_idx).SetBondType(
                        Chem.rdchem.BondType.AROMATIC
                    )
                    highlight_bond_colors[bond_idx] = (0, 0, 1)  # blue

                if (
                    change_dict[Change.BROKEN]
                    and not change_dict[Change.FORMED]
                ):
                    atoms_idx = [map_num_idx_dict[a] for a in bond]
                    bond_idx = mol.GetBondBetweenAtoms(*atoms_idx).GetIdx()
                    bonds_to_highlight.append(bond_idx)
                    mol.GetBondWithIdx(bond_idx).SetBondType(
                        Chem.rdchem.BondType.AROMATIC
                    )
                    highlight_bond_colors[bond_idx] = (1, 0, 0)  # red

                if change_dict[Change.FLEETING]:
                    atoms_idx = [map_num_idx_dict[a] for a in bond]
                    bond_idx = mol.GetBondBetweenAtoms(*atoms_idx).GetIdx()
                    bonds_to_highlight.append(bond_idx)
                    mol.GetBondWithIdx(bond_idx).SetBondType(
                        Chem.rdchem.BondType.AROMATIC
                    )
                    highlight_bond_colors[bond_idx] = (1, 0, 1)  # magenta

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
        ht = _HighlightTuple(atoms_to_highlight=atoms_to_highlight,
                             highlight_atom_colors=highlight_atom_colors,
                             bonds_to_highlight=bonds_to_highlight,
                             highlight_bond_colors=highlight_bond_colors)
        return mol, ht
    
    def svg(
        self,
        graph: (
            MolGraph
            | CondensedReactionGraph
            | StereoMolGraph
            | StereoCondensedReactionGraph
        )
    ) -> str:
        mol, ht = self._to_mol(graph)
        Chem.rdDepictor.Compute2DCoords(mol, useRingTemplates=True)  # type: ignore
        Chem.rdDepictor.StraightenDepiction(mol)  # type: ignore

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
