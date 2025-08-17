import pytest

from stereomolgraph.coords import Geometry
from stereomolgraph.graphs import (
    CondensedReactionGraph,
    MolGraph,
    StereoCondensedReactionGraph,
    StereoMolGraph,
)
from stereomolgraph.algorithms.isomorphism import vf2pp_all_isomorphisms


class TestIsomorphismMG:
    _TestClass: type[MolGraph] = MolGraph

    @pytest.fixture
    def water_graph(self, water_geo):
        return self._TestClass.from_geometry(water_geo)

    @pytest.fixture
    def mol_graph(self):
        mol_graph = self._TestClass()
        mol_graph.add_atom(0, atom_type="C")
        mol_graph.add_atom(1, atom_type="H")
        mol_graph.add_atom(2, atom_type="O")
        mol_graph.add_bond(0, 1, bond_order=1)
        return mol_graph
    
    @pytest.fixture
    def enantiomer_graph1(self, enantiomer_geos):
        return self._TestClass.from_geometry(enantiomer_geos[0])

    @pytest.fixture
    def enantiomer_graph2(self, enantiomer_geos):
        return self._TestClass.from_geometry(enantiomer_geos[1])

    def test_get_isomorphic_mappings(self, water_graph, mol_graph):
        assert [] == [
            i for i in vf2pp_all_isomorphisms(mol_graph, water_graph)
        ]

        assert all(
            mapping in ({0: 0, 1: 1, 2: 2}, {0: 0, 2: 1, 1: 2})
            for mapping in (
                i for i in vf2pp_all_isomorphisms(water_graph, water_graph)
            )
        )

    def test_get_isomorphic_mappings_of_enantiomers(
        self, enantiomer_graph1, enantiomer_graph2
    ):
        assert 2 == len(
            list(vf2pp_all_isomorphisms(enantiomer_graph1, enantiomer_graph2))
        )

    def test_get_automorphic_mappings(self, water_graph):
        assert all(
            mapping in ({0: 0, 1: 1, 2: 2}, {0: 0, 2: 1, 1: 2})
            for mapping in (
                i for i in vf2pp_all_isomorphisms(water_graph, water_graph)
            )
        )


class TestIsomorphismCRG(TestIsomorphismMG):
    _TestClass: type[CondensedReactionGraph] = CondensedReactionGraph


class TestIsomorphismSMG(TestIsomorphismMG):
    _TestClass: type[StereoMolGraph] = StereoMolGraph

    @pytest.fixture
    def chiral_product_graph1(self, data_path):
        filepath = (
            data_path
            / "disrot_reaction"
            / "(Z)-(4S)-3,4-Dichlor-2-pentene.xyz"
        )
        chiral_product_geo1 = Geometry.from_xyz_file(filepath)

        return self._TestClass.from_geometry(chiral_product_geo1)

    @pytest.fixture
    def chiral_product_graph2(self, data_path):
        filepath = (
            data_path / "conrot_reaction/(Z)-(4S)-3,4-Dichlor-2-pentene.xyz"
        )
        chiral_product_geo2 = Geometry.from_xyz_file(filepath)
        return self._TestClass.from_geometry(chiral_product_geo2)

    def test_get_atom_stereo_isomorphic_mappings(self, chiral_product_graph1):
        isomorphic_graph = chiral_product_graph1.copy()
        isomorphic_graph.relabel_atoms(
            {0: 1, 1: 0, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7}
        )
        m1 = {i: i for i in range(15)}
        m2 = {**m1, 11: 10, 12: 11, 10: 12}
        m3 = {**m1, 12: 10, 10: 11, 11: 12}
        m4 = {**m1, 7: 6, 8: 7, 6: 8}
        m5 = {**m1, 7: 6, 8: 7, 6: 8, 11: 10, 12: 11, 10: 12}
        m6 = {**m1, 7: 6, 8: 7, 6: 8, 12: 10, 10: 11, 11: 12}
        m7 = {**m1, 8: 6, 6: 7, 7: 8}
        m8 = {**m1, 8: 6, 6: 7, 7: 8, 11: 10, 12: 11, 10: 12}
        m9 = {**m1, 8: 6, 6: 7, 7: 8, 12: 10, 10: 11, 11: 12}
        assert {
            frozenset(mapping.items())
            for mapping in vf2pp_all_isomorphisms(
                chiral_product_graph1, isomorphic_graph, stereo=True
            )
        } == {
            frozenset(mapping.items())
            for mapping in (m1, m2, m3, m4, m5, m6, m7, m8, m9)
        }

    def test_get_isomorphic_mappings2(
        self, chiral_product_graph1, chiral_product_graph2
    ):
        assert chiral_product_graph1.is_isomorphic(chiral_product_graph2)
        assert chiral_product_graph1.is_isomorphic(chiral_product_graph2)
        gen_iso = vf2pp_all_isomorphisms(
            chiral_product_graph1,
            chiral_product_graph2,
            stereo=True,
        )
        assert len(list(gen_iso)) == 9


class TestIsomorphismSCRG(TestIsomorphismCRG, TestIsomorphismSMG):
    _TestClass: type[StereoCondensedReactionGraph] = (
        StereoCondensedReactionGraph
    )
