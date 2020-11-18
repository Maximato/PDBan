# PDB analysis script

import math
from Bio.PDB import PDBParser, Structure, Model, Chain


MOL_MASSES = {"O": 16, "C": 12, "N": 14, "S": 32}


def get_molecular_mass(struct: Structure) -> float:
    """
    Calculate molecular mass of protein structure
    :param struct: structure of protein
    :return: molecular mass
    """
    mass = 0
    for atom in struct.get_atoms():
        mass += MOL_MASSES[atom.get_name()[0]]
    return mass


def get_geometric_center(struct: Structure) -> list:
    """
    Calculate center of protein structure by arithmetic mean of all atoms coordinates
    :param struct: structure of protein
    :return: coordinates of center
    """
    x_sum, y_sum, z_sum = 0, 0, 0
    n = 0
    for atom in struct.get_atoms():
        n += 1
        x_sum += atom.get_coord()[0]
        y_sum += atom.get_coord()[1]
        z_sum += atom.get_coord()[2]
    return [x_sum/n, y_sum/n, z_sum/n]


def get_mass_center(struct: Structure) -> list:
    """
    Calculate mass center of protein structure
    :param struct: structure of protein
    :return: coordinates of mass center
    """
    x_sum, y_sum, z_sum = 0, 0, 0
    mass = get_molecular_mass(struct)
    for atom in struct.get_atoms():
        atom_mass = MOL_MASSES[atom.get_name()[0]]
        x_sum += atom.get_coord()[0]*atom_mass
        y_sum += atom.get_coord()[1]*atom_mass
        z_sum += atom.get_coord()[2]*atom_mass
    return [x_sum/mass, y_sum/mass, z_sum/mass]


def get_gyration_radius(struct: Structure) -> float:
    """
    Calculate radius of gyration protein
    https://en.wikipedia.org/wiki/Radius_of_gyration
    :param struct: structure of protein
    :return: radius of gyration
    """
    mass_center = get_mass_center(struct)

    coords = [a.get_coord() for a in struct.get_atoms()]
    masses = [MOL_MASSES[a.get_name()[0]] for a in struct.get_atoms()]

    mr2_summ = 0
    for r, m in zip(coords, masses):
        mr2_summ += m*(math.dist(r, mass_center))**2
    return math.sqrt(mr2_summ/sum(masses))


def get_structure_slice_by_residues(struct: Structure, start: int, finish: int, domain_name: str) -> Structure:
    """
    Return new structure that contains new model (id=1), new chain (id=1) with residues from 'start' to 'finish' of
    first chain input structure
    :param struct: input structure to slice
    :param start: start residue
    :param finish: finish residues
    :param domain_name: new structure name
    :return: new structure
    """
    new_chain = Chain.Chain(1)
    chain = list(struct.get_chains())[0]
    for i in range(start, finish+1):
        new_chain.add(chain[i])

    model = Model.Model(1)
    model.add(new_chain)
    domain = Structure.Structure(domain_name)
    domain.add(model)
    return domain


# build structures
pdb_parser = PDBParser()
structure = pdb_parser.get_structure("ts", "test.pdb")
domain1 = get_structure_slice_by_residues(structure, 1, 114, "domain1")
domain2 = get_structure_slice_by_residues(structure, 1, 114, "domain2")

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
