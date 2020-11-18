"""
Extract nanobody (2Rs15d is real structure of VHH needed for analysis)
from crystal structure of a HER2-Nb complex (5my6)
and compare it with another structures
"""

from Bio.PDB import PDBParser, Superimposer

from PDBA import get_structure_slice_by_residues, get_gyration_radius, get_synchronized_atoms


# build structures
pdb_parser = PDBParser()
real_structure = pdb_parser.get_structure("5my6", "pdb5my6.ent")
nanobody = get_structure_slice_by_residues(real_structure, 1, 1, 115, "2Rs15d")
nanobody_residues = list(nanobody.get_residues())

modeled_structure = pdb_parser.get_structure("5", "2Rs15d_1tag_model.pdb")
domain1 = get_structure_slice_by_residues(modeled_structure, 0, 1, 115, "domain1")
domain1_residues = list(domain1.get_residues())

domain2 = get_structure_slice_by_residues(modeled_structure, 0, 120, 234, "domain2")
domain2_residues = list(domain2.get_residues())

# calculation parameters
rg_nanobody = get_gyration_radius(nanobody)

sup = Superimposer()
sup.set_atoms(*get_synchronized_atoms(nanobody_residues, domain1_residues))
domain1_rms = sup.rms

sup.set_atoms(*get_synchronized_atoms(nanobody_residues, domain2_residues))
domain2_rms = sup.rms

# printing results
print(f"nanobody gyration radius: {rg_nanobody}")
print(f"domain1 RMSD: {domain1_rms}")
print(f"domain2 RMSD: {domain2_rms}")
