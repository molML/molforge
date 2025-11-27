from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from io import BytesIO

from ..constants import HAS_OPENEYE, oechem

def oe_to_rdkit(oe_mol: 'oechem.OEMol') -> Chem.Mol:
    """
    Convert OpenEye molecule with conformers to RDKit molecule.
    Preserves atom ordering by using SDF format.

    Args:
        oe_mol: OpenEye molecule with conformers

    Returns:
        RDKit molecule with all conformers

    Raises:
        ImportError: If OpenEye is not available
    """
    if not HAS_OPENEYE:
        raise ImportError(
            "OpenEye toolkit is required for oe_to_rdkit conversion. "
            "Install openeye-toolkits to use this function."
        )

    # Write first conformer to SDF to get the structure with correct atom ordering
    ofs = oechem.oemolostream()
    ofs.SetFormat(oechem.OEFormat_SDF)
    ofs.openstring()
    
    # Write just the first conformer to establish structure
    if oe_mol.NumConfs() > 0:
        first_conf = oe_mol.GetConf(oechem.OEHasConfIdx(0))
        oechem.OEWriteConstMolecule(ofs, first_conf)
    else:
        oechem.OEWriteMolecule(ofs, oe_mol)
    
    sdf_string = ofs.GetString()
    ofs.close()
    
    # Read into RDKit
    rd_mol = Chem.MolFromMolBlock(sdf_string, removeHs=False)
    if rd_mol is None:
        raise ValueError("Failed to create RDKit molecule from SDF")
    
    n_atoms = rd_mol.GetNumAtoms()
    n_confs = oe_mol.NumConfs()
    
    if n_atoms != oe_mol.NumAtoms():
        raise ValueError(f"Atom count mismatch: OE={oe_mol.NumAtoms()}, RDKit={n_atoms}")
    
    # Remove all conformers and add them back with correct coordinates
    rd_mol.RemoveAllConformers()
    
    coords_array = oechem.OEFloatArray(n_atoms * 3)
    
    for i, conf in enumerate(oe_mol.GetConfs()):
        conf.GetCoords(coords_array)
        coords = np.array(coords_array, dtype=np.float64).reshape(n_atoms, 3)
        
        # Create RDKit conformer
        rd_conf = Chem.Conformer(n_atoms)
        for atom_idx in range(n_atoms):
            rd_conf.SetAtomPosition(atom_idx, coords[atom_idx].tolist())
        
        rd_mol.AddConformer(rd_conf, assignId=True)
    
    return rd_mol

