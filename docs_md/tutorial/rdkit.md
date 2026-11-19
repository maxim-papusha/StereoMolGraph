# RDKit Interoperability

StereoMolGraphs(SMGs) can be used with `RDKit <https://www.rdkit.org/>`\_, a popular cheminformatics toolkit and we recommend its usage for parsing of InChIs or SMILES.
:doc:`graph2rdmol </reference/graph2rdmol>` and :doc:`rdmol2graph </reference/rdmol2graph>` functions allow for easy conversion between StereoMolGraph’s graph representations and RDKit’s `Mol` objects.

## Relevant details

- SMGs use explicit hydrogen atoms.
  -> use rdkit.Chem.AddHs bevore conversion
- Stereochemistry is converted.
- StereoMolGraph does not have bond orders nor charges and they are ignored
