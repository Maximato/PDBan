from code.processor import compare_domains


# configuration
structure1_id = "test"
s1_filename = "test_data/2Rs15d-link-2Rs15d.pdb"
d1_info = (0, 1, 115)

structure2_id = "test"
s2_filename = "test_data/2Rs15d-link-2Rs15d.pdb"
d2_info = (0, 120, 234)

if __name__ == '__main__':
    # compare two structures
    compare_domains(structure1_id, s1_filename, d1_info, structure2_id, s2_filename, d2_info)
