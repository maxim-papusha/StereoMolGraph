# pyright: standard
# typing with rdkit is not fully supported
from __future__ import annotations

from dataclasses import dataclass

import rdkit.Chem as Chem
from rdkit.Chem import rdChemReactions

from stereomolgraph import (
    CondensedReactionGraph,
    MolGraph,
    StereoCondensedReactionGraph,
    StereoMolGraph,
)
from stereomolgraph.rdmol2graph import RDMol2StereoMolGraph, mol_graph_from_rdmol

# ---------------------------------------------------------------------------
# Reaction SMILES → StereoCondensedReactionGraph
# ---------------------------------------------------------------------------


@dataclass
class RSMILES2StereoCondensedReactionGraph(RDMol2StereoMolGraph):
    """Convert a reaction SMILES string to a
    :class:`StereoCondensedReactionGraph`.

    **Requirements for the input reaction SMILES:**

    * All hydrogens must be **explicit** (no implicit hydrogens).
    * Every atom, **including every hydrogen**, must carry a non-zero atom
      map number.

    Inherits all stereochemistry-related configuration from
    :class:`RDMol2StereoMolGraph`.  By default ``use_atom_map_number`` is
    ``True`` because atom-mapped SMILES are the expected input.

    :param sanitize: If ``True``, sanitize the parsed molecules.
    """

    use_atom_map_number: bool = True
    sanitize: bool = True

    def __call__(self, rsmiles: str) -> StereoCondensedReactionGraph:
        """Convert *rsmiles* to a :class:`StereoCondensedReactionGraph`.

        :param rsmiles: Reaction SMILES string (``reactants>>products``)
            with explicit hydrogens and atom mapping on every atom.
        :return: A stereo condensed reaction graph.
        :raises ValueError: If any hydrogen is implicit, any atom lacks
            an atom map number, or the reactant and product atom sets do
            not match.
        """
        rdrxn = rdChemReactions.ReactionFromSmarts(rsmiles)
        if self.sanitize:
            rdChemReactions.SanitizeRxn(rdrxn)

        reactant_mols = list(rdrxn.GetReactants())
        product_mols = list(rdrxn.GetProducts())

        if not reactant_mols or not product_mols:
            raise ValueError(
                "Reaction SMILES must have at least one reactant and one product"
            )

        _validate_explicit_hydrogens(reactant_mols, product_mols)
        _enforce_atom_maps_for_reaction(reactant_mols, product_mols)

        # Convert each molecule using the parent's stereo-enabled __call__.
        # We call the unbound method explicitly so that self.__call__
        # (which now expects str) is not dispatched.
        _parent_call = RDMol2StereoMolGraph.__call__
        reactant_graphs = tuple(_parent_call(self, mol) for mol in reactant_mols)
        product_graphs = tuple(_parent_call(self, mol) for mol in product_mols)

        reactant_graph = (
            reactant_graphs[0]
            if len(reactant_graphs) == 1
            else StereoMolGraph.compose(reactant_graphs)
        )
        product_graph = (
            product_graphs[0]
            if len(product_graphs) == 1
            else StereoMolGraph.compose(product_graphs)
        )

        _validate_balanced_atom_sets(reactant_graph, product_graph, rsmiles)

        return StereoCondensedReactionGraph.from_graphs(
            reactant_graph=reactant_graph,
            product_graph=product_graph,
        )


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------


def _validate_explicit_hydrogens(
    reactants: list[Chem.Mol],
    products: list[Chem.Mol],
) -> None:
    """Raise :exc:`ValueError` if any molecule contains implicit hydrogens."""
    for mol in reactants + products:
        for atom in mol.GetAtoms():
            if atom.GetNumImplicitHs() > 0:
                raise ValueError(
                    "All hydrogens must be explicit in the reaction SMILES. "
                    "Found implicit hydrogens on atom "
                    f"{atom.GetSymbol()}{atom.GetIdx()}"
                    + (f" (map {atom.GetAtomMapNum()})" if atom.GetAtomMapNum() else "")
                    + "."
                )


def _enforce_atom_maps_for_reaction(
    reactants: list[Chem.Mol],
    products: list[Chem.Mol],
) -> None:
    """Raise :exc:`ValueError` if any atom (including hydrogen) in the
    reactant or product molecules lacks a non-zero atom map number."""
    for side, mols in [("reactant", reactants), ("product", products)]:
        for mol in mols:
            for atom in mol.GetAtoms():
                if atom.GetAtomMapNum() == 0:
                    raise ValueError(
                        f"Every atom must have a non-zero atom map number. "
                        f"Missing on {side} atom "
                        f"{atom.GetSymbol()}{atom.GetIdx()}."
                    )


def _validate_balanced_atom_sets(
    reactant_graph: MolGraph,
    product_graph: MolGraph,
    rsmiles: str,
) -> None:
    """Check that reactant and product graphs have identical
    ``(atom_id, element)`` sets.  Raise a helpful :exc:`ValueError`
    if they differ (e.g. unbalanced explicit hydrogens)."""
    r_set = set(zip(reactant_graph.atoms, reactant_graph.atom_types))
    p_set = set(zip(product_graph.atoms, product_graph.atom_types))
    if r_set == p_set:
        return
    only_r = r_set - p_set
    only_p = p_set - r_set
    msg = (
        f"Reactant and product atom sets do not match.\n"
        f"  Only in reactants: {sorted(only_r)}\n"
        f"  Only in products:  {sorted(only_p)}\n"
        "Ensure every atom on the reactant side has a "
        "corresponding atom (same map number and element) "
        "on the product side."
    )
    raise ValueError(msg)


def crg_from_rsmiles(
    rsmiles: str,
    *,
    sanitize: bool = True,
) -> CondensedReactionGraph:
    """Create a :class:`CondensedReactionGraph` from a reaction SMILES string.

    **Requirements for the input reaction SMILES:**

    * All hydrogens must be **explicit** (no implicit hydrogens).
    * Every atom, **including every hydrogen**, must carry a non-zero atom
      map number.

    :param rsmiles: Reaction SMILES string (``reactants>>products``).
    :param sanitize: If ``True``, sanitize the parsed molecules.
    :return: A :class:`CondensedReactionGraph` representing the reaction.
    """
    rdrxn = rdChemReactions.ReactionFromSmarts(rsmiles)
    if sanitize:
        rdChemReactions.SanitizeRxn(rdrxn)

    reactant_mols = list(rdrxn.GetReactants())
    product_mols = list(rdrxn.GetProducts())

    if not reactant_mols or not product_mols:
        raise ValueError(
            "Reaction SMILES must have at least one reactant and one product"
        )

    _validate_explicit_hydrogens(reactant_mols, product_mols)
    _enforce_atom_maps_for_reaction(reactant_mols, product_mols)

    reactant_graphs = tuple(
        mol_graph_from_rdmol(MolGraph, mol, use_atom_map_number=True)
        for mol in reactant_mols
    )
    product_graphs = tuple(
        mol_graph_from_rdmol(MolGraph, mol, use_atom_map_number=True)
        for mol in product_mols
    )

    reactant_graph = (
        reactant_graphs[0]
        if len(reactant_graphs) == 1
        else MolGraph.compose(reactant_graphs)
    )
    product_graph = (
        product_graphs[0]
        if len(product_graphs) == 1
        else MolGraph.compose(product_graphs)
    )

    _validate_balanced_atom_sets(reactant_graph, product_graph, rsmiles)

    return CondensedReactionGraph.from_graphs(
        reactant_graph=reactant_graph,
        product_graph=product_graph,
    )


def scrg_from_rsmiles(
    rsmiles: str,
    *,
    sanitize: bool = True,
    stereo_complete: bool = False,
    lone_pair_stereo: bool = True,
    resonance: bool = True,
) -> StereoCondensedReactionGraph:
    """Create a :class:`StereoCondensedReactionGraph` from a reaction SMILES
    string.

    **Requirements for the input reaction SMILES:**

    * All hydrogens must be **explicit** (no implicit hydrogens).
    * Every atom, **including every hydrogen**, must carry a non-zero atom
      map number.

    All stereochemistry-related parameters are forwarded to
    :class:`RDMol2StereoMolGraph`.

    :param rsmiles: Reaction SMILES string (``reactants>>products``).
    :param sanitize: If ``True``, sanitize the parsed molecules.
    :param stereo_complete: If ``True``, attempt to infer complete stereo
        parities when possible.
    :param lone_pair_stereo: Include stereochemistry that depends on lone
        pairs (if present).
    :param resonance: Enumerate resonance structures and merge bond stereo
        information from them.
    :return: A :class:`StereoCondensedReactionGraph` representing the reaction.
    """
    return RSMILES2StereoCondensedReactionGraph(
        sanitize=sanitize,
        stereo_complete=stereo_complete,
        lone_pair_stereo=lone_pair_stereo,
        resonance=resonance,
    )(rsmiles)
