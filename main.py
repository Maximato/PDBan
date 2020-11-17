# PDB analysis script

import math
from Bio.PDB import PDBParser, Structure


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


# build structure
pdb_parser = PDBParser()
structure = pdb_parser.get_structure("ts", "test.pdb")

print(get_geometric_center(structure))
print(get_mass_center(structure))
print(get_gyration_radius(structure))
