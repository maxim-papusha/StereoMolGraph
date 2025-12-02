from stereomolgraph import (MolGraph,
                            CondensedReactionGraph,
                            StereoMolGraph,
                            StereoCondensedReactionGraph)
from stereomolgraph.algorithms.color_refine import (color_refine_mg,
                                                    color_refine_crg,
                                                    color_refine_smg,
                                                    color_refine_scrg,
                                                    numpy_int_multiset_hash)


def consistent_hash_v1(g: MolGraph
                       | CondensedReactionGraph
                       | StereoMolGraph
                       | StereoCondensedReactionGraph) -> int:
    if g.n_atoms == 0:
        raise NotImplementedError("Hashing for empty graphs is not implemented"
                                  "consistently.")
    atom_labels = g.atom_types

    if isinstance(g, StereoCondensedReactionGraph):
        color_array = color_refine_scrg(g, atom_labels=atom_labels)
    elif isinstance(g, StereoMolGraph):
        color_array = color_refine_smg(g, atom_labels=atom_labels)
    elif isinstance(g, CondensedReactionGraph):
        color_array = color_refine_crg(g, atom_labels=atom_labels)
    elif isinstance(g, MolGraph):
        color_array = color_refine_mg(g, atom_labels=atom_labels)
    else:
        raise TypeError(f"Unsupported graph type: {type(g)}")
    return int(numpy_int_multiset_hash(color_array))