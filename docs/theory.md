# Theory


Molgraph: Connectivity

CondensedReactionGraph: Connectivity, ConnectivityChange

StereoMolGraph: Connectivity, Stereochemistry

StereoCondensedReactionGraph: Connectivity, Stereochemistry, ConnectivityChange, StereochemistryChange



| Graph Type                  | Connectivity | Connectivity <br> Change | Stereochemistry |Stereochemistry <br> Change  |
|-----------------------------|--------------|-------------|-----------------|-------------|
| `Molgraph`                  | ✓            |             |                 |             |
| `CondensedReactionGraph`    | ✓            | ✓          |                 |             |
| `StereoMolGraph`            | ✓            |             | ✓               |             |
| `StereoCondensedReactionGraph` | ✓         | ✓          | ✓               | ✓           |