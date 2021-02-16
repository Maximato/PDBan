### About
This project is designed to help in work with PDB files, which represent the three-dimensional structure of proteins. It allows to calculate some parameters of the protein structure: the geometric center, center of mass, gyration radius of domains, distance between domains. It is also possible to compare two structures by calculating the RMSD (root-mean-square deviation of atomic positions).

### Scripts
#### `calc_params.py`:
Script `calc_params.py` is designed to calculate parameters of two specified domains in a structure. To run this script set parameters:
- `structure_id` --- id of structure
- `pdb_filename` --- filename to PDB structure
- `d1_info` --- domain 1 info (turple: chain_number, start_residue, finish_residue)
- `d2_info` --- domain 2 info (turple: chain_number, start_residue, finish_residue)

After running this script the geometric center, center of mass, gyration radius of domains, distance between domains will be calculated.
Geometric center calculated as:

![formula](https://render.githubusercontent.com/render/math?math=r_c=\frac{\sum_{i=1}^{N}r_i}{N})

Mass center calculated as:

![formula](https://render.githubusercontent.com/render/math?math=r_c=\frac{\sum_{i=1}^{N}m_i r_i}{\sum_{i=1}^{N}m_i})

Gyration radius calculated as:

![formula](https://render.githubusercontent.com/render/math?math=R_g=\sqrt{\frac{\sum_{i=1}^{N}m_i (r_i-r_c)^2}{\sum_{i=1}^{N}m_i}})

Where: *N* - number of atoms, m_i - mass of *i* atom, r_i - coordinates of *i* atom.

Distance calculated as the Euclidean distance between coordinates of mass centers of domains.

#### `compare_domains.py`:
Script `compare_domains.py` is designed to compare the two structures by calculating the RMSD (root-mean-square deviation of atomic positions). To run this script set parameters:
- `structure1_id`, `structure2_id` --- id of structure
- `s1_filename`, `s2_filename` --- filename to PDB structure
- `d1_info`, `d2_info` --- domain info (turple: chain_number, start_residue, finish_residue)

After running this script the gyration radius of domains and the RMSD will be calculated. RMSD calculation is implemented using the **set_atoms** function from the **Biopython** library.

### Example 1
To test scripts you can run it for modeled structure of nanobody **2Rs15d-GGGS-2Rs15d** where 2Rs15d is sequences of VHH antibody *Camelus dromedarius*, GGGS - linker. Coordinates of first subunit is 1-115 residues, coordinates of second subunit is 120-234 residues.

__`calc_params.py` configuration:__
```
# configuration
structure_id = "test"
pdb_filename = "test_data/2Rs15d-link-2Rs15d.pdb"
d1_info = (0, 1, 115)
d2_info = (0, 120, 234)
```

__`compare_domains.py` configuration:__
```
# configuration
structure1_id = "test"
s1_filename = "test_data/2Rs15d-link-2Rs15d.pdb"
d1_info = (0, 1, 115)

structure2_id = "test"
s2_filename = "test_data/2Rs15d-link-2Rs15d.pdb"
d2_info = (0, 120, 234)
```

### Example 2
There is real structure of VHH 2Rs15d in protein data bank (B chain of HER2-Nb complex [5my6](https://www.rcsb.org/structure/5MY6)). This structure was download and saved in **5my6.ent** file. You can compare crystal structure of 2Rs15d (1 chain, 1-115 residues of 5my6 complex) with modeled structure 2Rs15d (0 chain, 1-115 or 120-234 residues of **2Rs15d-GGGS-2Rs15d** model) by using `compare_domains.py` script.

__`compare_domains.py` configuration:__
```
# configuration
structure1_id = "modeled"
s1_filename = "test_data/2Rs15d-link-2Rs15d.pdb"
d1_info = (0, 1, 115)

structure2_id = "crystal_structure"
s2_filename = "test_data/5my6.ent"
d2_info = (1, 1, 115)
```
