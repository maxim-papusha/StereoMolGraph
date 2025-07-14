# pyright: standard
from __future__ import annotations

from rdkit import Chem # type: ignore
from rdkit.Chem import Draw # type: ignore
from typing import NamedTuple

from stereomolgraph import (
        MolGraph,
        CondensedReactionGraph,
        StereoMolGraph,
        StereoCondensedReactionGraph)
from stereomolgraph.stereodescriptors import PlanarBond

def default_repr_svg(graph):
    return View2D().svg(graph)

def default_view_molgraph(self: MolGraph) -> None:
    View2D()(self)

MolGraph._ipython_display_ = default_view_molgraph


class View2D(NamedTuple):
    height: int = 300
    width: int = 300
    show_atom_numbers: bool = True
    show_h: bool = True
    generate_bond_orders: bool = False
    dummy_atoms: bool = True

    def svg(
        self,
        graph: (
            MolGraph
            | CondensedReactionGraph
            | StereoMolGraph
            | StereoCondensedReactionGraph
        ),
    ) -> str:
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

        if isinstance(graph, StereoMolGraph):
            for db in graph.bond_stereo.values():
                if not isinstance(db, PlanarBond):
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
                    Chem.rdchem.BondType.HYDROGEN)
                highlight_bond_colors[bond_idx] = (0, 0, 1)  # blue

            for bond in graph.get_broken_bonds():
                atoms_idx = [map_num_idx_dict[a] for a in bond]
                bond_idx = mol.GetBondBetweenAtoms(*atoms_idx).GetIdx()
                bonds_to_highlight.append(bond_idx)
                mol.GetBondWithIdx(bond_idx).SetBondType(
                    Chem.rdchem.BondType.HYDROGEN)
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

        Chem.rdDepictor.Compute2DCoords(mol, useRingTemplates=True) # type: ignore
        Chem.rdDepictor.StraightenDepiction(mol) # type: ignore

        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(self.width, self.height)

        drawer.drawOptions().useBWAtomPalette()
        drawer.drawOptions().continuousHighlight = False
        drawer.drawOptions().highlightBondWidthMultiplier = 10
        drawer.drawOptions().fillHighlights = False
        drawer.drawOptions().includeRadicals = False

        drawer.DrawMolecule(
            mol,
            highlightAtoms=atoms_to_highlight,
            highlightAtomColors=highlight_atom_colors,
            highlightBonds=bonds_to_highlight,
            highlightBondColors=highlight_bond_colors,
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
