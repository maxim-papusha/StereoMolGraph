"""Isotope-exchange stereoisomer enumeration regression tests."""

import pytest

from stereomolgraph import StereoMolGraph
from stereomolgraph.experimental import generate_isotope_stereoisomers
from stereomolgraph.stereodescriptors import PlanarBond, Tetrahedral


def _molecule(elements, bonds, isotopes=None):
    graph = StereoMolGraph()
    isotopes = isotopes or {}
    for atom, element in enumerate(elements):
        graph.add_atom(atom, element, isotope=isotopes.get(atom, 0))
    for a, b in bonds:
        graph.add_bond(a, b)
    return graph


def _methane():
    return _molecule(["C", "H", "H", "H", "H"], [(0, i) for i in range(1, 5)])


def _chlorofluoromethane():
    return _molecule(
        ["C", "H", "H", "F", "Cl"],
        [(0, i) for i in range(1, 5)],
    )


def _ethanol():
    # CH3-CH2-OH, all hydrogens explicit
    return _molecule(
        ["C", "C", "O", "H", "H", "H", "H", "H", "H"],
        [(0, 1), (1, 2), (0, 3), (0, 4), (0, 5),
         (1, 6), (1, 7), (2, 8)],
    )


def test_methane_mono_through_tetrasubstitution():
    methane = _methane()
    for n in range(1, 5):
        products = list(generate_isotope_stereoisomers(
            methane, "H", 2, min_substitutions=n, max_substitutions=n
        ))
        assert len(products) == 1
        assert products[0].isotopes.count(2) == n
        assert products[0].get_atom_stereo(0) is None
    assert all(isotope == 0 for isotope in methane.isotopes)


def test_enantiomers_from_monodeuterated_chlorofluoromethane():
    graph = _chlorofluoromethane()
    products = list(generate_isotope_stereoisomers(
        graph, "H", 2, max_substitutions=1
    ))
    assert len(products) == 2
    assert all(isinstance(p.get_atom_stereo(0), Tetrahedral) for p in products)
    assert products[0] != products[1]
    assert products[0] == products[1].enantiomer()


def test_dideuterated_chlorofluoromethane_is_not_chiral():
    products = list(generate_isotope_stereoisomers(
        _chlorofluoromethane(), "H", 2,
        min_substitutions=2, max_substitutions=2
    ))
    assert len(products) == 1
    assert products[0].get_atom_stereo(0) is None


def test_ethanol_has_four_monodeuterated_stereoisomers():
    products = list(generate_isotope_stereoisomers(
        _ethanol(), "H", 2, max_substitutions=1
    ))
    assert len(products) == 4
    assert all(p.isotopes.count(2) == 1 for p in products)
    assert len([p for p in products if p.get_atom_stereo(1)]) == 2


def test_fixed_stereo_is_not_inverted_by_exchange():
    # The central carbon has four distinct ligands, and fixed configuration.
    graph = _molecule(
        ["C", "F", "Cl", "Br", "C", "H", "H", "H"],
        [(0, 1), (0, 2), (0, 3), (0, 4),
         (4, 5), (4, 6), (4, 7)],
    )
    fixed = Tetrahedral((0, 1, 2, 3, 4), 1)
    graph.set_atom_stereo(fixed)
    products = list(generate_isotope_stereoisomers(graph, "H", 2,
                                                    max_substitutions=1))
    assert len(products) == 1
    assert products[0].get_atom_stereo(0) == fixed
    assert products[0] != products[0].enantiomer()


def test_preexisting_explicit_isotope_can_be_exchanged():
    graph = _molecule(
        ["C", "H", "H", "H", "H"],
        [(0, i) for i in range(1, 5)],
        isotopes={1: 1, 2: 1, 3: 1, 4: 1},
    )
    products = list(generate_isotope_stereoisomers(
        graph, "H", 2, source_isotope=1, max_substitutions=1
    ))
    assert len(products) == 1
    assert sorted(products[0].isotopes) == [0, 1, 1, 1, 2]


def test_allowed_atoms_limits_exchange_sites():
    products = list(generate_isotope_stereoisomers(
        _ethanol(), "H", 2, max_substitutions=1, atoms=(6, 7)
    ))
    assert len(products) == 2


def test_isotope_induced_ez_stereoisomerism():
    # CH2=CHCl (explicit H); only replace a hydrogen on the CH2 end.
    graph = _molecule(
        ["C", "C", "H", "H", "H", "Cl"],
        [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)],
    )
    graph.set_bond_stereo(PlanarBond((2, 3, 0, 1, 4, 5), None))
    products = list(generate_isotope_stereoisomers(
        graph, "H", 2, max_substitutions=1, atoms=(2, 3)
    ))
    assert len(products) == 2
    assert all(p.get_bond_stereo({0, 1}) is not None for p in products)


def test_requires_explicit_hydrogen_atoms():
    graph = _molecule(["C"], [])
    assert list(generate_isotope_stereoisomers(graph, "H", 2,
                                               min_substitutions=0)) == [graph]


@pytest.mark.parametrize("target", [0, -1, 1.5])
def test_invalid_target_isotope(target):
    with pytest.raises(ValueError):
        list(generate_isotope_stereoisomers(_methane(), "H", target))


def test_input_graph_can_be_frozen():
    graph = _chlorofluoromethane().freeze()
    assert len(list(generate_isotope_stereoisomers(
        graph, "H", 2, max_substitutions=1
    ))) == 2
