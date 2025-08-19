from pathlib import Path
import shelve

from hypothesis import assume, given
import hypothesis.strategies as st
from hypothesis.stateful import Bundle, RuleBasedStateMachine, initialize, rule, invariant, precondition

from stereomolgraph.periodic_table import PERIODIC_TABLE
from stereomolgraph import MolGraph


class GraphModifier(RuleBasedStateMachine):
    def __init__(self):
        super().__init__()
        self.g = MolGraph()
        #self._visited_path = Path("test_shelve")
        #self.visited = shelve.open(str(self._visited_path))

    #def teardown(self):
    #    self.visited.close()
    #    if self._visited_path.exists():
    #        self._visited_path.unlink()


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

    @invariant()
    def do_nothing(self):
        assert True

GraphTest = GraphModifier.TestCase