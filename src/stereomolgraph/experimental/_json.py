from __future__ import annotations

import json
from typing import Any

from stereomolgraph import (
    CondensedReactionGraph,
    MolGraph,
    StereoCondensedReactionGraph,
    StereoMolGraph,
)
from stereomolgraph.periodic_table import SYMBOLS
from stereomolgraph.stereodescriptors import (
    AtropBond,
    Octahedral,
    PlanarBond,
    SquarePlanar,
    Tetrahedral,
    TrigonalBipyramidal,
)


STEREO_CLASSES: dict[str, type] = {
    "Tetrahedral": Tetrahedral,
    "TrigonalBipyramidal": TrigonalBipyramidal,
    "Octahedral": Octahedral,
    "SquarePlanar": SquarePlanar,
    "PlanarBond": PlanarBond,
    "AtropBond": AtropBond,
}


class JSONHandler:
    """Serialize and deserialize StereoMolGraph-related graphs to JSON."""

    @staticmethod
    def as_dict(graph: MolGraph) -> dict[str, Any]:
        graph_type = type(graph).__name__

        atoms_dict = sorted(
            (
                atom,
                SYMBOLS[a_type],
            )
            for atom, a_type in zip(graph.atoms, graph.atom_types)
        )

        data: dict[str, Any] = {"Atoms": atoms_dict}

        if isinstance(graph, CondensedReactionGraph):
            formed_bonds = graph.get_formed_bonds()
            broken_bonds = graph.get_broken_bonds()

            bonds_dict = sorted(
                sorted(bond)
                for bond in graph.bonds
                if bond not in formed_bonds | broken_bonds
            )
            data["Bonds"] = bonds_dict
            data["Formed Bonds"] = sorted(
                tuple(sorted(bond)) for bond in formed_bonds
            )
            data["Broken Bonds"] = sorted(
                tuple(sorted(bond)) for bond in broken_bonds
            )
        else:
            bonds_dict = sorted(tuple(sorted(bond)) for bond in graph.bonds)
            data["Bonds"] = bonds_dict

        if isinstance(graph, StereoMolGraph):
            atom_stereos: dict[Any, Any] = {}
            bond_stereos: dict[Any, Any] = {}

            for atom, stereo in graph.atom_stereo.items():
                atom_stereos[atom] = {
                    stereo.__class__.__name__: (stereo.atoms, stereo.parity)
                }
            if atom_stereos:
                data["Atom Stereo"] = atom_stereos

            for bond_fset, stereo in graph.bond_stereo.items():
                bond_key = json.dumps(sorted(bond_fset))
                bond_stereos[bond_key] = {
                    stereo.__class__.__name__: (stereo.atoms, stereo.parity)
                }
            if bond_stereos:
                data["Bond Stereo"] = bond_stereos

        if isinstance(graph, StereoCondensedReactionGraph):
            atom_changes: dict[Any, Any] = {}

            for atom, change_dict in graph.atom_stereo_changes.items():
                atom_change_dict: dict[str, Any] = {}
                for change, stereo in change_dict.items():
                    if stereo is not None:
                        atom_change_dict[change.name] = {
                            stereo.__class__.__name__: (
                                stereo.atoms,
                                stereo.parity,
                            )
                        }
                if atom_change_dict:
                    atom_changes[atom] = atom_change_dict
            if atom_changes:
                data["Atom Stereo Changes"] = atom_changes

            bond_changes: dict[Any, Any] = {}
            for bond_fset, change_dict in graph.bond_stereo_changes.items():
                bond_key = json.dumps(sorted(bond_fset))
                bond_change_dict: dict[str, Any] = {}
                for change, stereo in change_dict.items():
                    if stereo is not None:
                        bond_change_dict[change.name] = {
                            stereo.__class__.__name__: (
                                stereo.atoms,
                                stereo.parity,
                            )
                        }
                if bond_change_dict:
                    bond_changes[bond_key] = bond_change_dict
            if bond_changes:
                data["Bond Stereo Changes"] = bond_changes

        return {graph_type: data}

    @classmethod
    def json_serialize(cls, graph: MolGraph) -> str:
        return json.dumps(cls.as_dict(graph))

    @staticmethod
    def _stereo_from_payload(payload: Any):
        if not payload:
            return None
        class_name, (atoms, parity) = next(iter(payload.items()))
        stereo_cls = STEREO_CLASSES[class_name]
        return stereo_cls(tuple(int(a) for a in atoms), parity)

    @classmethod
    def json_deserialize(cls, payload: str) -> MolGraph:
        graph_type, graph_payload = next(iter(json.loads(payload).items()))
        graph_cls = {
            "MolGraph": MolGraph,
            "StereoMolGraph": StereoMolGraph,
            "CondensedReactionGraph": CondensedReactionGraph,
            "StereoCondensedReactionGraph": StereoCondensedReactionGraph,
        }[graph_type]

        graph = graph_cls()

        for atom_id, atom_type in graph_payload.get("Atoms", []):
            graph.add_atom(int(atom_id), atom_type)

        for bond_entry in graph_payload.get("Bonds", []):
            graph.add_bond(*map(int, bond_entry))

        if isinstance(graph, CondensedReactionGraph):
            for bond_entry in graph_payload.get("Formed Bonds", []):
                graph.add_formed_bond(*map(int, bond_entry))
            for bond_entry in graph_payload.get("Broken Bonds", []):
                graph.add_broken_bond(*map(int, bond_entry))

        if isinstance(graph, StereoMolGraph):
            for entry in graph_payload.get("Atom Stereo", {}).values():
                stereo_obj = cls._stereo_from_payload(entry)
                if stereo_obj is not None:
                    graph.set_atom_stereo(stereo_obj)

            for entry in graph_payload.get("Bond Stereo", {}).values():
                stereo_obj = cls._stereo_from_payload(entry)
                if stereo_obj is not None:
                    graph.set_bond_stereo(stereo_obj)

        if isinstance(graph, StereoCondensedReactionGraph):
            for change_dict in graph_payload.get(
                "Atom Stereo Changes", {}
            ).values():
                broken = cls._stereo_from_payload(change_dict.get("BROKEN"))
                formed = cls._stereo_from_payload(change_dict.get("FORMED"))
                fleeting = cls._stereo_from_payload(
                    change_dict.get("FLEETING")
                )
                if any((broken, formed, fleeting)):
                    graph.set_atom_stereo_change(
                        broken=broken,
                        formed=formed,
                        fleeting=fleeting,
                    )

            for change_dict in graph_payload.get(
                "Bond Stereo Changes", {}
            ).values():
                broken = cls._stereo_from_payload(change_dict.get("BROKEN"))
                formed = cls._stereo_from_payload(change_dict.get("FORMED"))
                fleeting = cls._stereo_from_payload(
                    change_dict.get("FLEETING")
                )
                if any((broken, formed, fleeting)):
                    graph.set_bond_stereo_change(
                        broken=broken,
                        formed=formed,
                        fleeting=fleeting,
                    )

        return graph
