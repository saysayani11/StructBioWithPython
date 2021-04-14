# Bio.PDB is a Biopython module that focuses on working with crystal structures of biological macromolecules
# The glob module finds all the pathnames matching a specified pattern according to the rules used by the Unix shell, although results are returned in arbitrary order. 
from Bio.PDB import MMCIFParser
from glob import glob
parser = MMCIFParser()

# \U in "C:\Users... starts an eight-character Unicode escape, such as \U00014321
# You either need to duplicate all backslashes or prefix the string with r (to produce a raw string)

# Assuming the downloaded files are saed by their default filenames
# pdb_files = glob(r'required_file_location\*')
pdb_files = glob(r'C:\Users\saySa\OneDrive\Desktop\pdbfiles\*')

#slice the 4 letter structure name (for example, '5p21') and the full file name (for example, '5p21.cif')
for x in pdb_files:
     structure_name=x[-8:-4]
     file_name=x[-8:]
     # structure = parser.get_structure('5p21,'5p21.cif')
     structure = parser.get_structure(structure_name,file_name)
     # Rest of your code
