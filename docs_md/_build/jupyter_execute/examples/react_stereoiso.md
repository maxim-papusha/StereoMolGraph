# Stereochemistry of Phosgenation Transition States

```ipython3
import rdkit.Chem as Chem

from stereomolgraph import StereoCondensedReactionGraph
from stereomolgraph.stereodescriptors import Tetrahedral
from stereomolgraph.experimental import generate_fleeting_stereoisomers
from stereomolgraph.ipython import View2D

from IPython.core.interactiveshell import InteractiveShell

InteractiveShell.ast_node_interactivity = "all"

View2D.show_atom_numbers = True
View2D.generate_bond_orders = False
```

```ipython3
tmda_rdmol = Chem.AddHs(Chem.MolFromSmiles("NCCCCN"))
phosgene_rdmol = Chem.AddHs(Chem.MolFromSmiles("ClC(Cl)=O"))

scrg1 = StereoCondensedReactionGraph.from_rdmol(
    Chem.CombineMols(tmda_rdmol, phosgene_rdmol)
)
scrg1
```

```ipython3
scrg1.add_formed_bond(0, 19)
scrg1.add_formed_bond(18, 7)
scrg1.add_broken_bond(18, 19)
scrg1.add_broken_bond(0, 7)
scrg1
```

```ipython3
n_neighbors = scrg1.bonded_to(0)
c_neighbors = scrg1.bonded_to(19)
scrg1.set_atom_stereo_change(fleeting=Tetrahedral(atoms=(0, *n_neighbors), parity=None))
scrg1.set_atom_stereo_change(
    fleeting=Tetrahedral(atoms=(19, *c_neighbors), parity=None)
)
```

```ipython3
isomers = list(generate_fleeting_stereoisomers(scrg1))
len(isomers)
for iso in isomers:
    iso
```

```ipython3
scrg2 = StereoCondensedReactionGraph.from_rdmol(
    Chem.CombineMols(tmda_rdmol, phosgene_rdmol)
)
scrg2
```

```ipython3
scrg2.add_formed_bond(0, 19)
scrg2.add_formed_bond(18, 17)
scrg2.add_formed_bond(5, 7)

scrg2.add_broken_bond(18, 19)
scrg2.add_broken_bond(5, 17)
scrg2.add_broken_bond(0, 7)
scrg2
```

```ipython3
scrg2.set_atom_stereo_change(fleeting=Tetrahedral(atoms=(5, *scrg2.bonded_to(5))))
scrg2.set_atom_stereo_change(fleeting=Tetrahedral(atoms=(0, *scrg2.bonded_to(0))))
scrg2.set_atom_stereo_change(fleeting=Tetrahedral(atoms=(19, *scrg2.bonded_to(19))))
```

```ipython3
isomers2 = list(generate_fleeting_stereoisomers(scrg2))
len(isomers2)

for iso in isomers2:
    iso
```
