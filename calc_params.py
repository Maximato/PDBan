from code.processor import calc_domains_params

# configuration
structure_id = "test"
pdb_filename = "test_data/2Rs15d-link-2Rs15d.pdb"
d1_info = (0, 1, 115)
d2_info = (0, 120, 234)


if __name__ == '__main__':
    # model parameters calculation
    calc_domains_params(structure_id, pdb_filename, d1_info, d2_info)
