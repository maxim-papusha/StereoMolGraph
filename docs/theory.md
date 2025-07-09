# Theory


Molgraph: Connectivity

CondensedReactionGraph: Connectivity, ConnectivityChange

StereoMolGraph: Connectivity, Stereochemistry

StereoCondensedReactionGraph: Connectivity, Stereochemistry, ConnectivityChange, StereochemistryChange




| Graph Type                     | Connectivity | ConnectivityChange | Stereochemistry | StereochemistryChange |
|--------------------------------|--------------|--------------------|-----------------|-----------------------|
| Molgraph                       | ✓           |                    |                 |                       |
| CondensedReactionGraph         | ✓           | ✓                  |                 |                       |
| StereoMolGraph                 | ✓           |                    | ✓               |                       |
| StereoCondensedReactionGraph   | ✓           | ✓                  | ✓               | ✓                     |