# Theory


Molgraph: Connectivity

CondensedReactionGraph: Connectivity, ConnectivityChange

StereoMolGraph: Connectivity, Stereochemistry

StereoCondensedReactionGraph: Connectivity, Stereochemistry, ConnectivityChange, StereochemistryChange



| Graph Type                  | Connectivity | Connectivity <br> Change | Stereo |Stereo <br> Change  |
|-----------------------------|--------------|-------------|-----------------|-------------|
| `Molgraph`                  | ✓            |             |                 |             |
| `CondensedReactionGraph`    | ✓            | ✓          |                 |             |
| `StereoMolGraph`            | ✓            |             | ✓               |             |
| `StereoCondensedReactionGraph` | ✓         | ✓          | ✓               | ✓           |