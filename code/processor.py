from Bio.PDB import PDBParser, Superimposer

from code.PDBA import *


def calc_domains_params(id: str, pdb_filename: str, d1_info: tuple, d2_info: tuple):
    """
    Calculate some parameters of two specified domains in a structure
    :param id: structure id
    :param pdb_filename: filename to PDB structure
    :param d1_info: domain 1 info (turple: chain_number, start_residue, finish_residue)
    :param d2_info: domain 2 info (turple: chain_number, start_residue, finish_residue)
    """

    # build structures
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(id, pdb_filename)
    domain1 = get_structure_slice_by_residues(structure, "domain1", *d1_info)
    domain2 = get_structure_slice_by_residues(structure, "domain2", *d2_info)

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


def compare_domains(s1_id: str, s1_filename: str, ds1_info: tuple, s2_id: str, s2_filename: str, ds2_info: tuple):
    """
    Extract nanobody (2Rs15d is real structure of VHH needed for analysis)
    from crystal structure of a HER2-Nb complex (5my6)
    and compare it with another structures
    Compare two PDB structures by calculating their gyration radius and RMSD
    PDB structures should be equal primary protein structures: equal size, equal sequences
    :param s1_id: structure 1 id
    :param s1_filename: filename to PDB structure 1
    :param ds1_info: domain info in structure 1 (turple: chain_number, start_residue, finish_residue)
    :param s2_id: structure 2 id
    :param s2_filename: filename to PDB structure 2
    :param ds2_info: domain info in structure 2 (turple: chain_number, start_residue, finish_residue)
    """

    # build structures
    pdb_parser = PDBParser()

    structure1 = pdb_parser.get_structure(s1_id, s1_filename)
    domain1 = get_structure_slice_by_residues(structure1, "real_structure", *ds1_info)
    domain1_residues = list(domain1.get_residues())

    structure2 = pdb_parser.get_structure(s2_id, s2_filename)
    domain2 = get_structure_slice_by_residues(structure2, "domain1", *ds2_info)
    domain2_residues = list(domain2.get_residues())

    # calculation parameters
    rg1 = get_gyration_radius(domain1)
    rg2 = get_gyration_radius(domain2)

    sup = Superimposer()
    sup.set_atoms(*get_synchronized_atoms(domain1_residues, domain2_residues))
    domain_rms = sup.rms

    # printing results
    print(f"domain 1 gyration radius: {rg1}")
    print(f"domain 2 gyration radius: {rg2}")
    print(f"domain RMSD: {domain_rms}")
