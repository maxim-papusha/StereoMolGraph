from __future__ import annotations

import itertools
import json
from collections import deque, Counter
from collections.abc import Iterable, Iterator
from typing import TYPE_CHECKING, Any

from stereomolgraph import (
    AtomId,
    Bond,
    CondensedReactionGraph,
    MolGraph,
    StereoCondensedReactionGraph,
    StereoMolGraph,
)
from stereomolgraph.algorithms.color_refine import color_refine_smg
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms
from stereomolgraph.graphs.crg import Change
from stereomolgraph.periodic_table import SYMBOLS
from stereomolgraph.stereodescriptors import (
    AtropBond,
    Octahedral,
    PlanarBond,
    SquarePlanar,
    Tetrahedral,
    TrigonalBipyramidal,
    NonRotatableBond
)

BondAutoCls = int
MappingId = int

STEREO_CLASSES: dict[str, type] = {
    "Tetrahedral": Tetrahedral,
    "TrigonalBipyramidal": TrigonalBipyramidal,
    "Octahedral": Octahedral,
    "SquarePlanar": SquarePlanar,
    "PlanarBond": PlanarBond,
    "AtropBond": AtropBond,
    "NonRotatableBond": NonRotatableBond,
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


if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

    from stereomolgraph.graphs.smg import AtomId, Bond, StereoMolGraph


def unique_generator(input_generator: Iterator) -> Iterator:
    """
    A generator that yields unique objects from another generator.

    Args:
        input_generator: A generator yielding hashable objects

    Yields:
        Only the first occurrence of each unique object from the input generator
    """
    seen_hash = set()
    for item in input_generator:
        item_hash = hash(item)
        if item_hash not in seen_hash:
            seen_hash.add(item_hash)
            yield item


def generate_stereoisomers(
    graph: StereoMolGraph,
    enantiomers: bool = True,
    atoms: None | Iterable[AtomId] = None,
    bonds: None | Iterable[Bond] = None,
) -> Iterator[StereoMolGraph]:
    """Generates all unique stereoisomers of a StereoMolGraph by generation of
    all combinations of parities. Only includes stereocenters which have a
    parity of None. If a parity is set, it is not changed.

    If include_enantiomers is True, both enantiomers of a stereoisomer are
    included, if it is False, only one enantiomer is included.

    Args:
        enantiomers: If True, both enantiomers are included
        atoms: Optional subset of atoms to consider for stereoisomerism
        bonds: Optional subset of bonds to consider for stereoisomerism

    Yields:
        StereoMolGraph: Each unique stereoisomer (and enantiomer if requested)
    """
    if atoms is None:
        atom_stereos = (
            stereo.get_isomers()
            for a in graph.atoms
            if ((stereo := graph.get_atom_stereo(a)) and stereo.parity is None)
        )
    else:
        atom_stereos = (
            stereo.get_isomers()
            for a in atoms
            if (stereo := graph.get_atom_stereo(a)) is not None
        )

    if bonds is None:
        bond_stereos = (
            stereo.get_isomers()
            for b in graph.bonds
            if ((stereo := graph.get_bond_stereo(b)) and stereo.parity is None)
        )
    else:
        bond_stereos = (
            stereo.get_isomers()
            for b in bonds
            if (stereo := graph.get_bond_stereo(b)) is not None
        )

    seen_hash: set[int] = set()
    enantiomers_seen_hash: set[int] = set()

    for a_stereos, b_stereos in itertools.product(
        itertools.product(*atom_stereos), itertools.product(*bond_stereos)
    ):
        stereoisomer = graph.copy()
        for a_stereo in a_stereos:
            stereoisomer.set_atom_stereo(a_stereo)
        for b_stereo in b_stereos:
            stereoisomer.set_bond_stereo(b_stereo)

        stereoisomer_hash = hash(stereoisomer)
        if stereoisomer_hash not in seen_hash:
            seen_hash.add(stereoisomer_hash)
            yield stereoisomer

            if not enantiomers:
                enantiomer = stereoisomer.enantiomer()
                enantiomer_hash = hash(enantiomer)
                if enantiomer_hash not in enantiomers_seen_hash:
                    enantiomers_seen_hash.add(enantiomer_hash)
                    yield enantiomer


def generate_fleeting_stereoisomers(
    graph: StereoCondensedReactionGraph,
    enantiomers: bool = True,
    atoms: None | Iterable[AtomId] = None,
    bonds: None | Iterable[Bond] = None,
) -> Iterator[StereoCondensedReactionGraph]:
    """Generates all unique fleeting stereoisomers of a StereoCondensedReactionGraph.

    Only includes stereocenters which have a parity of None for the fleeting change.
    If a parity is set, it is not changed.

    Args:
        graph: The reaction graph to generate isomers from
        enantiomers: If True, both enantiomers are included (default: True)
        atoms: Optional subset of atoms to consider for stereoisomerism
        bonds: Optional subset of bonds to consider for stereoisomerism

    Yields:
        StereoCondensedReactionGraph: Each unique fleeting stereoisomer
    """
    # Get atom stereoisomers
    if atoms is None:
        atom_stereos = (
            stereo.get_isomers()
            for a in graph.atoms
            if (
                (stereo_change_dict := graph.get_atom_stereo_change(a))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        )
    else:
        atom_stereos = (
            stereo.get_isomers()
            for a in atoms
            if (
                (stereo_change_dict := graph.get_atom_stereo_change(a))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        )

    # Get bond stereoisomers
    if bonds is None:
        bond_stereos = (
            stereo.get_isomers()
            for b in graph.bonds
            if (
                (stereo_change_dict := graph.get_bond_stereo_change(b))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        )
    else:
        bond_stereos = (
            stereo.get_isomers()
            for b in bonds
            if (
                (stereo_change_dict := graph.get_bond_stereo_change(b))
                and (stereo := stereo_change_dict[Change.FLEETING])
            )
            and stereo.parity is None
        )

    seen_isomers_hash = set()
    seen_enantiomers_hash = set()

    for a_stereos, b_stereos in itertools.product(
        itertools.product(*atom_stereos), itertools.product(*bond_stereos)
    ):
        stereoisomer = graph.copy()
        for a_stereo in a_stereos:
            stereoisomer.set_atom_stereo_change(fleeting=a_stereo)
        for b_stereo in b_stereos:
            stereoisomer.set_bond_stereo_change(fleeting=b_stereo)

        stereoisomer_hash = hash(stereoisomer)
        if stereoisomer_hash not in seen_isomers_hash:
            seen_isomers_hash.add(stereoisomer_hash)
            yield stereoisomer

            if not enantiomers:
                enantiomer = stereoisomer.enantiomer()
                enantiomer_hash = hash(enantiomer)
                if enantiomer_hash not in seen_enantiomers_hash:
                    seen_enantiomers_hash.add(enantiomer_hash)
                    yield enantiomer


def topological_symmetry_number(graph: StereoMolGraph) -> int:
    """
    Calculated from the number of graph isomorphisms which conserve the
    stereo information.
    topological_symmetry_number = external_symmetry_number * Product(internal_symmetry_numbers)
    TODO: add paper reference
    """

    if any(stereo.parity is None for stereo in graph.stereo.values()):
        raise NotImplementedError(
            "all stereocenters have to be defined"
            " to calculate the symmetry number"
        )
    #colorings = color_refine_smg(graph)
    mappings = vf2pp_all_isomorphisms(
        graph, graph,
        #atom_labels=(colorings, colorings),
        stereo=True
    )
    return deque(enumerate(mappings, 1), maxlen=1)[0][0]


def external_symmetry_number(graph: StereoMolGraph) -> int:
    """
    TODO: add paper reference.
    """
    if any(stereo.parity is None for stereo in graph.stereo.values()):
        raise NotImplementedError(
            "all stereocenters have to be defined"
            " to calculate the symmetry number"
        )
    # colorings = color_refine_smg(graph)
    mappings = tuple(
        vf2pp_all_isomorphisms(
            graph, graph, stereo=True,
            #atom_labels=(colorings, colorings),
        )
    )

    atom_auto_set: dict[int, set[int]] = {
        a: {a} for a in graph.atoms
    }

    bond_auto_set: dict[Bond, set[Bond]] = {
        b: {b} for b in graph.bonds
    }

    for mapping in mappings:
        for a1, a2 in mapping.items():
            if a2 in atom_auto_set[a1]:
                continue
            else:
                atom_auto_set[a1].update(
                    atom_auto_set[a2]
                )
                atom_auto_set[a2] = atom_auto_set[a1]

        for b1 in graph.bonds:
            b2 = Bond({mapping[a] for a in b1})
            if b2 in bond_auto_set[b1]:
                continue
            else:
                bond_auto_set[b1].update(
                    bond_auto_set[b2]
                )
                bond_auto_set[b2] = bond_auto_set[b1]

    atom_auto_int: dict[int, int] = {
        a: min(atom_set) for a, atom_set in atom_auto_set.items()
    }
    bond_auto_int: dict[Bond, int] = {
        bond: hash(frozenset(bond_set)) for bond, bond_set in bond_auto_set.items()
    }

    #int_atom_auto_set: dict[int, set[int]] = {auto_cls_int: atom_auto_set[a]
    #                                          for a, auto_cls_int
    #                                          in atom_auto_int.items()}
    int_bond_auto_set: dict[int, set[Bond]] = {auto_cls_int: bond_auto_set[b]
                                             for b, auto_cls_int
                                             in bond_auto_int.items()}
    
    counter: dict[BondAutoCls,
                  dict[frozenset[NonRotatableBond], set[MappingId]]] = {}

    for bond_auto_cls, bonds in int_bond_auto_set.items():

        counter[bond_auto_cls] = {}

        bond_to_state: dict[Bond, set[NonRotatableBond]] = {}

        for bond in bonds:
            if graph.get_bond_stereo(bond):
                continue
            a1, a2 = tuple(bond)
            nbrs1 = set(graph.bonded_to(a1))
            nbrs1.remove(a2)
            nbrs2 = set(graph.bonded_to(a2))
            nbrs2.remove(a1)
            if len(nbrs1) != 3 or len(nbrs2) != 3:
                continue

            nrb_set: set[NonRotatableBond] = set()
            if (len({atom_auto_int[nbr] for nbr in nbrs1}) == 2 and
                len({atom_auto_int[nbr] for nbr in nbrs2}) == 2):

                unique_left_int = next(v for v, c in Counter(atom_auto_int[x] for x in nbrs1).items() if c == 1)
                unique_right_int = next(v for v, c in Counter(atom_auto_int[x] for x in nbrs2).items() if c == 1)
                unique_left = next(nbr for nbr in nbrs1 if atom_auto_int[nbr] == unique_left_int)
                unique_right = next(nbr for nbr in nbrs2 if atom_auto_int[nbr] == unique_right_int)

                s1 = graph.get_atom_stereo(a1)
                s2 = graph.get_atom_stereo(a2)
                assert isinstance(s1, Tetrahedral)
                assert isinstance(s2, Tetrahedral)

                left = next((stereo_atoms[2:5] if s1.parity == 1
                          else stereo_atoms[2:5][::-1])
                          for stereo_atoms in s1._perm_atoms()
                          if stereo_atoms[1] == a2 and stereo_atoms[2] == unique_left)
            
                right = next((stereo_atoms[2:5] if s2.parity == -1
                          else stereo_atoms[2:5][::-1])
                          for stereo_atoms in s2._perm_atoms()
                          if stereo_atoms[1] == a1 and stereo_atoms[2] == unique_right)

                nrb_set.add(NonRotatableBond(atoms=(*left, a1, a2, *right), parity=0))

            elif (len({atom_auto_int[nbr] for nbr in nbrs1}) == 1 or
                len({atom_auto_int[nbr] for nbr in nbrs2}) == 1):
                s1 = graph.get_atom_stereo(a1)
                s2 = graph.get_atom_stereo(a2)
                assert isinstance(s1, Tetrahedral)
                assert isinstance(s2, Tetrahedral)

                left = next((stereo_atoms[2:5] if s1.parity == 1
                          else stereo_atoms[2:5][::-1])
                          for stereo_atoms in s1._perm_atoms()
                          if stereo_atoms[1] == a2 )
            
                right = next((stereo_atoms[2:5] if s2.parity == -1
                          else stereo_atoms[2:5][::-1])
                          for stereo_atoms in s2._perm_atoms()
                          if stereo_atoms[1] == a1 )

                nrb_set.add(NonRotatableBond(atoms=(*left, a1, a2, *right), parity=0))

                #for perm_nbrs1 in itertools.permutations(nbrs1):
                #    for perm_nbrs2 in itertools.permutations(nbrs2):
                #        nrb = NonRotatableBond.from_bond_and_nrbs(perm_nbrs1, a1, a2, perm_nbrs2)
                #        nrb_set.add(nrb)
            else:
                continue
                raise Exception(len({atom_auto_int[nbr] for nbr in nbrs1}), len({atom_auto_int[nbr] for nbr in nbrs2}))
            
            bond_to_state[bond] = nrb_set

        for _, comb in enumerate(itertools.product(*bond_to_state.values())):

            counter[bond_auto_cls][frozenset(comb)] = set()

    for mapping_id, mapping in enumerate(mappings):

        for bond_auto_cls, bond_state_dict in counter.items():
            for bond_state, mapping_ids in bond_state_dict.items():
                if all(NonRotatableBond(tuple(mapping[a] for a in non_rot_bond.atoms), parity=0) in bond_state
                       for non_rot_bond in bond_state):
                    mapping_ids.add(mapping_id)

        
    # pick one bond state per auto class, maximize the intersection size of mapping ids
    if not counter:
        return 1
    mapping_sets_by_cls = [list(bond_state_dict.values()) for bond_state_dict in counter.values() if bond_state_dict]

    sym_number = 1
    for combo in itertools.product(*mapping_sets_by_cls):
        intersect: set[MappingId] = set.intersection(*combo) if combo else set()
        sym_number = max(sym_number, len(intersect))
        if sym_number == len(mappings):
            break

    return sym_number


            










