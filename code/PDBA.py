# PDB analysis (PDBA)

import math
from Bio.PDB import Structure, Model, Chain, Atom
from Bio.PDB.PDBExceptions import PDBException


MOL_MASSES = {"O": 16, "C": 12, "N": 14, "S": 32, "H": 1}


def get_molecular_mass(struct: Structure) -> float:
    """
    Calculate molecular mass of protein structure
    :param struct: structure of protein
    :return: molecular mass
    """
    mass = 0
    for atom in struct.get_atoms():
        mass += _get_atom_mass(atom)
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
        atom_mass = _get_atom_mass(atom)
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
    masses = [_get_atom_mass(a) for a in struct.get_atoms()]

    mr2_sum = 0
    for r, m in zip(coords, masses):
        mr2_sum += m*(math.dist(r, mass_center))**2
    return math.sqrt(mr2_sum/sum(masses))


def get_structure_slice_by_residues(struct: Structure, domain_name: str, chain_order: int, start: int, finish: int) -> Structure:
    """
    Return new structure that contains new model (id=1), new chain (id=1) with residues from 'start' to 'finish' of
    specified chain of input structure
    :param struct: input structure to slice
    :param chain_order: order of chain to extract residues
    :param start: start residue
    :param finish: finish residues
    :param domain_name: new structure name
    :return: new structure
    """
    new_chain = Chain.Chain(1)
    chain = list(struct.get_chains())[chain_order]
    for i in range(start, finish+1):
        new_chain.add(chain[i])

    model = Model.Model(1)
    model.add(new_chain)
    domain = Structure.Structure(domain_name)
    domain.add(model)
    return domain


def get_synchronized_atoms(residues1: list, residues2: list) -> (list, list):
    """
    Synchronize atoms in residues1 and residues2.
    Needs for superimposing structures (class Superimposer, set_atoms function)
    :param residues1: list of residues
    :param residues2: list of residues
    :return: two lists of synchronized_atoms
    """

    # different size checking
    if len(residues1) != len(residues2):
        raise PDBException("Residues differ in size")

    atoms1, atoms2 = [], []
    for r1, r2 in zip(residues1, residues2):
        # different residues checking
        if r1.get_resname() != r2.get_resname():
            raise PDBException(f"Different residues: {r1} and {r2}")

        r1_atoms = r1.get_list()
        r2_atoms = r2.get_list()

        if len(r1_atoms) == len(r2_atoms):
            atoms1.extend(r1_atoms)
            atoms2.extend(r2_atoms)
        else:
            # atom synchronization
            for a1 in r1_atoms:
                for a2 in r2_atoms:
                    if a1.get_name() == a2.get_name():
                        atoms1.append(a1)
                        atoms2.append(a2)
                        break
    return atoms1, atoms2


def _get_atom_mass(atom: Atom):
    i = 0
    atom_name = atom.get_name()
    while atom_name[i].isdigit():
        i += 1
    return MOL_MASSES[atom_name[i]]
