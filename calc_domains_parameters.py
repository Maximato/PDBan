"""
Calculate some parameters of two specified domains in a structure
"""

from Bio.PDB import PDBParser

from PDBA import *


# build structures
pdb_parser = PDBParser()
structure = pdb_parser.get_structure("2Rs15d", "2Rs15d_1tag_model.pdb")
domain1 = get_structure_slice_by_residues(structure, 0, 1, 115, "domain1")
domain2 = get_structure_slice_by_residues(structure, 0, 120, 234, "domain2")

# calculation parameters
center_domain1 = get_geometric_center(domain1)
center_domain2 = get_geometric_center(domain2)
rg_domain1 = get_gyration_radius(domain1)
rg_domain2 = get_gyration_radius(domain2)
domain_dist = math.dist(center_domain1, center_domain2)

# printing results
print(" == Domain1 == ")
print(f"mass center: {center_domain1}")
print(f"gyration radius: {rg_domain1}\n")

print(" == Domain2 == ")
print(f"mass center: {center_domain2}")
print(f"gyration radius: {rg_domain2}\n")

print("== Domain distance ==")
print(f"distance (A): {domain_dist}")
