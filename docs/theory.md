# Theory


Molgraph: Connectivity

CondensedReactionGraph: Connectivity, ConnectivityChange

StereoMolGraph: Connectivity, Stereochemistry

StereoCondensedReactionGraph: Connectivity, Stereochemistry, ConnectivityChange, StereochemistryChange



| Graph Type                  | Connectivity | Connectivity/Change | Stereochemistry |Stereochemistry/Change  |
|-----------------------------|--------------|-------------|-----------------|-------------|
| `Molgraph`                  | ✓            |             |                 |             |
| `CondensedReactionGraph`    | ✓            | ✓          |                 |             |
| `StereoMolGraph`            | ✓            |             | ✓               |             |
| `StereoCondensedReactionGraph` | ✓         | ✓          | ✓               | ✓           |