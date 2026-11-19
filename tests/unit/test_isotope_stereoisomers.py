"""Isotope-exchange stereoisomer enumeration regression tests."""

import pytest

from stereomolgraph import (
    IsotopeExchange,
    StereoMolGraph,
    generate_isotope_stereoisomers,
)
from stereomolgraph.stereodescriptors import PlanarBond, Tetrahedral


def _single_exchange(graph, atom_type, isotope, *,
                     source_isotope=0, atoms=None, **kwargs):
    """Adapter for the existing single-exchange regression tests."""
    return generate_isotope_stereoisomers(
        graph,
        IsotopeExchange(
            atom_type=atom_type,
            isotopes=(isotope,),
            source_isotope=source_isotope,
            positions=atoms,
        ),
        **kwargs,
    )


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
        products = list(_single_exchange(
            methane, "H", 2, substitutions=n
        ))
        assert len(products) == 1
        assert products[0].isotopes.count(2) == n
        assert products[0].get_atom_stereo(0) is None
    assert all(isotope == 0 for isotope in methane.isotopes)


def test_enantiomers_from_monodeuterated_chlorofluoromethane():
    graph = _chlorofluoromethane()
    products = list(_single_exchange(
        graph, "H", 2, substitutions=1
    ))
    assert len(products) == 2
    assert all(isinstance(p.get_atom_stereo(0), Tetrahedral) for p in products)
    assert products[0] != products[1]
    assert products[0] == products[1].enantiomer()


def test_dideuterated_chlorofluoromethane_is_not_chiral():
    products = list(_single_exchange(
        _chlorofluoromethane(), "H", 2,
        substitutions=2
    ))
    assert len(products) == 1
    assert products[0].get_atom_stereo(0) is None


def test_ethanol_has_four_monodeuterated_stereoisomers():
    products = list(_single_exchange(
        _ethanol(), "H", 2, substitutions=1
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
    products = list(_single_exchange(graph, "H", 2,
                                                    substitutions=1))
    assert len(products) == 1
    assert products[0].get_atom_stereo(0) == fixed
    assert products[0] != products[0].enantiomer()


def test_preexisting_explicit_isotope_can_be_exchanged():
    graph = _molecule(
        ["C", "H", "H", "H", "H"],
        [(0, i) for i in range(1, 5)],
        isotopes={1: 1, 2: 1, 3: 1, 4: 1},
    )
    products = list(_single_exchange(
        graph, "H", 2, source_isotope=1, substitutions=1
    ))
    assert len(products) == 1
    assert sorted(products[0].isotopes) == [0, 1, 1, 1, 2]


def test_allowed_atoms_limits_exchange_sites():
    products = list(_single_exchange(
        _ethanol(), "H", 2, substitutions=1, atoms=(6, 7)
    ))
    assert len(products) == 2


def test_isotope_induced_ez_stereoisomerism():
    # CH2=CHCl (explicit H); only replace a hydrogen on the CH2 end.
    graph = _molecule(
        ["C", "C", "H", "H", "H", "Cl"],
        [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)],
    )
    graph.set_bond_stereo(PlanarBond((2, 3, 0, 1, 4, 5), None))
    products = list(_single_exchange(
        graph, "H", 2, substitutions=1, atoms=(2, 3)
    ))
    assert len(products) == 2
    assert all(p.get_bond_stereo({0, 1}) is not None for p in products)


def test_requires_explicit_hydrogen_atoms():
    graph = _molecule(["C"], [])
    assert list(_single_exchange(graph, "H", 2,
                                               substitutions=0)) == [graph]


@pytest.mark.parametrize("target", [0, -1, 1.5])
def test_invalid_target_isotope(target):
    with pytest.raises(ValueError):
        list(_single_exchange(_methane(), "H", target))


def test_input_graph_can_be_frozen():
    graph = _chlorofluoromethane().freeze()
    assert len(list(_single_exchange(
        graph, "H", 2, substitutions=1
    ))) == 2


def test_multiple_target_isotopes_are_distinct_and_unique():
    products = list(generate_isotope_stereoisomers(
        _methane(),
        (IsotopeExchange("H", (2, 3)),),
        substitutions=2,
    ))
    # CH2D2, CH2DT, CH2T2, all achiral.
    assert len(products) == 3
    assert len({p.freeze() for p in products}) == 3
    assert sorted(
        tuple(sorted(m for m in p.isotopes if m)) for p in products
    ) == [(2, 2), (2, 3), (3, 3)]


def test_two_elements_can_be_exchanged_together():
    products = list(generate_isotope_stereoisomers(
        _methane(),
        (
            IsotopeExchange("H", (2,)),
            IsotopeExchange("C", (13,), positions={0}),
        ),
        substitutions=2,
    ))
    assert len(products) == 2
    assert sorted((p.get_atom_isotope(0), p.isotopes.count(2)) for p in products) == [
        (0, 2), (13, 1)
    ]


def test_unequal_allowed_positions_do_not_prune_two_step_path():
    # T is allowed only at site 1, D at sites 1 or 2.
    # D at 1 blocks T; D at 2 leaves T at 1 possible.
    products = list(generate_isotope_stereoisomers(
        _methane(),
        (
            IsotopeExchange("H", (2,), positions={1, 2}),
            IsotopeExchange("H", (3,), positions={1}),
        ),
        substitutions=2,
    ))
    assert len(products) == 2
    assert sorted(
        tuple(sorted(m for m in p.isotopes if m)) for p in products
    ) == [(2, 2), (2, 3)]


def test_same_product_from_overlapping_rules_is_emitted_once():
    products = list(generate_isotope_stereoisomers(
        _methane(),
        (
            IsotopeExchange("H", (2,), positions={1, 2}),
            IsotopeExchange(1, (2,), positions={2, 3}),
        ),
        substitutions=2,
    ))
    assert len(products) == 1
    assert len({p.freeze() for p in products}) == 1
    assert products[0].isotopes.count(2) == 2


def test_different_allowed_targets_on_equivalent_atoms():
    products = list(generate_isotope_stereoisomers(
        _methane(),
        (
            IsotopeExchange("H", (2,), positions={1}),
            IsotopeExchange("H", (3,), positions={2}),
        ),
        substitutions=2,
    ))
    assert len(products) == 1
    assert sorted(m for m in products[0].isotopes if m) == [2, 3]


def test_allowed_positions_refer_to_atom_ids():
    graph = _methane().relabel_atoms({1: 100, 2: 200, 3: 300, 4: 400})
    products = list(generate_isotope_stereoisomers(
        graph,
        IsotopeExchange(atom_type="H", isotopes=(2,), positions={100, 200}),
        substitutions=2,
    ))
    assert len(products) == 1
    assert products[0].get_atom_isotope(100) == 2
    assert products[0].get_atom_isotope(200) == 2


def test_validate_exchange_rule_and_positions():
    for targets in ((), (0,), (-1,), (2, 0)):
        with pytest.raises(ValueError):
            IsotopeExchange("H", targets)
    with pytest.raises(ValueError):
        IsotopeExchange("H", (2,), source_isotope=2)
    with pytest.raises(ValueError):
        IsotopeExchange("Xx", (2,))
    with pytest.raises(ValueError):
        list(generate_isotope_stereoisomers(
            _methane(), IsotopeExchange("H", (2,), positions={999})
        ))


def test_partial_products_are_not_emitted():
    results = list(generate_isotope_stereoisomers(
        _methane(),
        IsotopeExchange("H", (2,)),
        substitutions=1,
    ))
    assert len(results) == 1
    assert results[0].isotopes.count(2) == 1


def test_graph_method_and_frozen_input():
    graph = _methane().freeze()
    products = list(graph.generate_isotope_stereoisomers(
        IsotopeExchange("H", (2, 3)),
        substitutions=2,
    ))
    assert len(products) == 3  # DD, DT, TT
    assert all(p.isotopes.count(2) + p.isotopes.count(3) == 2 for p in products)
    assert graph.isotopes == (0, 0, 0, 0, 0)


def test_default_returns_only_fully_exchanged_methane():
    parent = _methane()
    products = list(parent.generate_isotope_stereoisomers(
        IsotopeExchange("H", (2,))
    ))
    assert len(products) == 1
    assert products[0].isotopes.count(2) == 4
    assert parent.isotopes == (0, 0, 0, 0, 0)


def test_default_full_exchange_for_multiple_elements():
    products = list(generate_isotope_stereoisomers(
        _methane(), (IsotopeExchange("H", (2,)), IsotopeExchange("C", (13,)))
    ))
    assert len(products) == 1
    assert products[0].get_atom_isotope(0) == 13
    assert products[0].isotopes.count(2) == 4


def test_default_full_exchange_for_two_target_isotopes():
    products = list(generate_isotope_stereoisomers(
        _methane(), IsotopeExchange("H", (2, 3))
    ))
    assert len(products) == 5  # DDDD, DDDT, DDTT, DTTT, TTTT
    assert all(p.isotopes.count(2) + p.isotopes.count(3) == 4 for p in products)


def test_default_full_exchange_only_at_permitted_atom_ids():
    products = list(generate_isotope_stereoisomers(
        _ethanol(), IsotopeExchange("H", (2,), positions={6, 7})
    ))
    assert len(products) == 1
    assert products[0].get_atom_isotope(6) == 2
    assert products[0].get_atom_isotope(7) == 2


def test_zero_substitutions_returns_original_when_explicit():
    parent = _methane()
    products = list(generate_isotope_stereoisomers(
        parent, IsotopeExchange("H", (2,)), substitutions=0
    ))
    assert len(products) == 1
    assert products[0] == parent
    assert products[0].isotopes.count(2) == 0


@pytest.mark.parametrize("substitutions", (-1, 5, 1.5, True))
def test_invalid_exact_substitution_count(substitutions):
    with pytest.raises(ValueError):
        list(generate_isotope_stereoisomers(
            _methane(), IsotopeExchange("H", (2,)),
            substitutions=substitutions,
        ))


def test_no_matching_atoms_yields_no_products_by_default():
    graph = _molecule(["C"], [])
    assert list(generate_isotope_stereoisomers(
        graph, IsotopeExchange("H", (2,))
    )) == []


def test_easy_top_level_imports_have_matching_types():
    from stereomolgraph.experimental import IsotopeExchange as ExperimentalExchange
    from stereomolgraph.experimental import (
        generate_isotope_stereoisomers as experimental_generate,
    )

    assert IsotopeExchange is ExperimentalExchange
    assert generate_isotope_stereoisomers is experimental_generate


def test_renumbered_input_preserves_counts():
    original = _ethanol()
    relabeled = original.relabel_atoms({atom: atom + 100 for atom in original.atoms})
    exchange = IsotopeExchange("H", (2,))
    for molecule in (original, relabeled):
        products = molecule.generate_isotope_stereoisomers(exchange, substitutions=1)
        assert len(list(products)) == 4
