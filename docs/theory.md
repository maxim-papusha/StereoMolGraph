# Theory


Molgraph: Connectivity

CondensedReactionGraph: Connectivity, ConnectivityChange

StereoMolGraph: Connectivity, Stereochemistry

StereoCondensedReactionGraph: Connectivity, Stereochemistry, ConnectivityChange, StereochemistryChange



| Graph Type     <br>         | Connectivity <br>| Connectivity <br> Change | Stereo <br> |Stereo <br> Change  |
|-----------------------------|------------------|--------------------------|-----------------|-------------|
| `Molgraph`                  | ✓                |                          |                 |             |
| `CondensedReactionGraph`    | ✓                | ✓                        |                 |             |
| `StereoMolGraph`            | ✓                |                          | ✓               |             |
| `StereoCondensedReactionGraph` | ✓             | ✓                       | ✓               | ✓           |