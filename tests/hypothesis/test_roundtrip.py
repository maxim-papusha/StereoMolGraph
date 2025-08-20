import random
from pathlib import Path
import shelve

from hypothesis import assume, given
import hypothesis.strategies as st
from hypothesis.stateful import Bundle, RuleBasedStateMachine, initialize, rule, invariant, precondition

from stereomolgraph.periodic_table import PERIODIC_TABLE
from stereomolgraph import MolGraph


class MolGraphModifier(RuleBasedStateMachine):
    _TestClass = MolGraph

    def __init__(self):
        super().__init__()
        self.g = self._TestClass()

    @rule(atom_id=st.integers(),
          atom_type=st.sampled_from(list(PERIODIC_TABLE.values())))
    def add_atom(self, atom_id, atom_type):
        assume(atom_id not in self.g.atoms)
        self.g.add_atom(atom_id, atom_type)

    @precondition(lambda self: len(self.g.atoms) != 0)
    @rule(data=st.data())
    def add_bond(self, data):
        atom_id1 = data.draw(st.sampled_from(list(self.g.atoms)))
        atom_id2 = data.draw(st.sampled_from(list(self.g.atoms)))
        assume(atom_id1 != atom_id2)
        self.g.add_bond(atom_id1, atom_id2)


class StereoMolGraphModifier(MolGraphModifier):
    _TestClass = MolGraph


class HashCollisionTester(MolGraphModifier):
    @invariant()
    def do_nothing(self):

        other_g = self._TestClass()
        id_type_list = list(zip(self.g.atoms, self.g.atom_types))
        bond_list = list(self.g.bonds)
        for _ in range(10):
            random.shuffle(id_type_list)
            random.shuffle(bond_list)
            for a, a_type in id_type_list:
                other_g.add_atom(a, a_type)
            for bond in bond_list:
                other_g.add_bond(*bond)

            assert hash(self.g) == hash(other_g)


HashTest = HashCollisionTester.TestCase