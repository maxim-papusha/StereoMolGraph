"""Embed a StereoMolGraph into 3D Cartesian coordinates using RDKit."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem
from rdkit import Geometry as RDGeometry
from rdkit import __version__ as rdkit_version
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers

from stereomolgraph import StereoMolGraph
from stereomolgraph.coords import Geometry
from stereomolgraph.graph2rdmol import stereo_mol_graph_to_rdmol


@dataclass
class EmbedParameters:
    """Parameters for RDKit distance-geometry embedding.

    All fields mirror the options available in RDKit's
    :class:`rdkit.Chem.rdDistGeom.EmbedParameters`.
    """

    seed: int = 42
    """Random seed for reproducible embeddings."""

    max_attempts: int = 1_000_000_000
    """Maximum number of attempts for the distance-geometry embedder."""

    n_conformers: int = 1
    """Number of conformers to generate per call."""

    optimize: bool = True
    """If ``True``, each conformer is refined with a UFF force-field
    minimization."""

    ETversion: int = 2
    """ETKDG version (1 or 2)."""

    useRandomCoords: bool = False
    """Use random coordinates as starting point."""

    useBasicKnowledge: bool = True
    """Use basic knowledge (bond lengths, angles) during embedding."""

    useExpTorsionAnglePrefs: bool = True
    """Use experimental torsion-angle preferences."""

    enforceChirality: bool = True
    """Enforce chirality constraints during embedding."""

    embedFragmentsSeparately: bool = False
    """Embed disconnected fragments separately."""

    useSmallRingTorsions: bool = False
    """Use small-ring torsion corrections."""

    useMacrocycleTorsions: bool = False
    """Use macrocycle torsion corrections."""

    def to_rdkit(self) -> rdDistGeom.EmbedParameters:
        """Build an RDKit :class:`~rdkit.Chem.rdDistGeom.EmbedParameters`
        from this dataclass."""
        ps = rdDistGeom.EmbedParameters()
        ps.ETversion = self.ETversion
        ps.randomSeed = self.seed
        ps.useRandomCoords = self.useRandomCoords
        ps.useBasicKnowledge = self.useBasicKnowledge
        ps.useExpTorsionAnglePrefs = self.useExpTorsionAnglePrefs
        ps.enforceChirality = self.enforceChirality
        ps.embedFragmentsSeparately = self.embedFragmentsSeparately
        ps.useSmallRingTorsions = self.useSmallRingTorsions
        ps.useMacrocycleTorsions = self.useMacrocycleTorsions
        # maxAttempts was renamed to maxIterations in RDKit >= 2025
        if int(rdkit_version.split(".")[0]) >= 2025:
            ps.maxIterations = self.max_attempts
        else:
            ps.maxAttempts = self.max_attempts  # type: ignore[attr-defined]
        return ps


@dataclass
class SMG2Geo:
    """Generate 3D geometries from :class:`StereoMolGraph` objects using RDKit.

    This is a callable dataclass: the :meth:`__call__` method is a
    **generator** that yields :class:`~stereomolgraph.coords.Geometry`
    instances, one per conformer.

    Usage::

        embed = SMG2Geo(seed=123)
        for geo in embed(smg):
            print(geo.xyz_str())

        # single geometry with custom parameters
        geo = next(SMG2Geo(params=EmbedParameters(optimize=False))(smg))

        # with fixed atom coordinates
        geo = next(SMG2Geo()(smg, fixed_coords={1: (0.0, 0.0, 0.0)}))
    """

    params: EmbedParameters = field(default_factory=EmbedParameters)
    """Embedding parameters passed to RDKit."""

    def __call__(
        self,
        smg: StereoMolGraph,
        *,
        fixed_coords: dict[int, tuple[float, float, float]] | None = None,
    ) -> Iterator[Geometry]:
        """Yield 3D geometries generated from *smg*.

        :param smg: input molecular graph (all hydrogens must be explicit)
        :param fixed_coords: optional mapping of SMG atom IDs to fixed
            Cartesian coordinates ``(x, y, z)`` in Angstrom.  Atoms listed
            here are constrained to the given positions during embedding;
            all other atoms are embedded freely by the distance-geometry
            engine.
        :yields: :class:`Geometry` objects, one per conformer
        :raises ValueError: if RDKit cannot embed the molecule
        """
        # --- Convert SMG → RDKit molecule -----------------------------------
        rw_mol, idx_map = stereo_mol_graph_to_rdmol(smg)
        mol: Chem.Mol = rw_mol.GetMol()
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)

        # --- Embed conformers ------------------------------------------------
        ps = self.params.to_rdkit()
        n_confs = self.params.n_conformers

        # Apply fixed-coordinate constraints (SMG atom IDs → RDKit indices)
        if fixed_coords is not None:
            smg_to_rdkit = {v: k for k, v in idx_map.items()}
            coord_map: dict[int, RDGeometry.Point3D] = {}
            for smg_id, xyz in fixed_coords.items():
                if smg_id not in smg_to_rdkit:
                    raise ValueError(
                        f"Fixed atom ID {smg_id} not found in the molecule."
                    )
                coord_map[smg_to_rdkit[smg_id]] = RDGeometry.Point3D(*xyz)
            ps.SetCoordMap(coord_map)

        if n_confs <= 1:
            result = rdDistGeom.EmbedMolecule(mol, ps)
            if result < 0:
                raise ValueError(
                    f"RDKit could not embed the molecule (error code {result})."
                )
            conf_ids = [mol.GetConformer().GetId()]
        else:
            conf_ids = list(
                rdDistGeom.EmbedMultipleConfs(mol, numConfs=n_confs, params=ps)
            )
            if not conf_ids:
                raise ValueError(
                    "RDKit could not embed any conformer for the molecule."
                )

        # --- Yield each conformer --------------------------------------------
        for cid in conf_ids:
            conf = mol.GetConformer(cid)

            if self.params.optimize:
                Chem.SetHybridization(mol)
                # rdForceFieldHelpers.UFFHasAllMoleculeParams
                ff = rdForceFieldHelpers.UFFGetMoleculeForceField(mol, confId=cid)
                if ff:
                    ff.Minimize(energyTol=1e-6, maxIts=100)
                    # Re-fetch positions after optimization
                    conf = mol.GetConformer(cid)

            coords = np.array(conf.GetPositions(), dtype=np.float64)
            yield Geometry(atom_types=smg.atom_types, coords=coords)
