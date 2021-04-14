# The glob module finds all the pathnames matching a specified pattern according to the rules used by the Unix shell, although results are returned in arbitrary order. 
# Pandas provide a fast and efficient DataFrame object for data manipulation with integrated indexing
# Bio.PDB is a Biopython module that focuses on working with crystal structures of biological macromolecules
import matplotlib.pyplot as plt
from glob import glob
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIFParser
parser = MMCIFParser()

pdb_files = glob(r'C:\Users\saySa\OneDrive\Desktop\pdbfiles\*')
no_of_files=len(pdb_files)

# minB calculates the minimum B-factor value in a PDB structure and returns a value
def minB(structure):
    b_factor=[]
    for model in structure:
          for chain in model:
              for residue in chain:
                  for atom in residue:
                      b_factor.append(atom.get_bfactor())
    return(min(b_factor))
    
# maxB calculates the maximum B-factor value in a PDB structure and returns a value
def maxB(structure): 
    b_factor=[]
    for model in structure:
          for chain in model:
              for residue in chain:
                  for atom in residue:
                      b_factor.append(atom.get_bfactor())
    return(max(b_factor))
    
min_b_list=[]
max_b_list=[]
index_names=[]

# Slice the 4 letter structure name (for example, '5p21') and the full file name (for example, '5p21.cif')
# Call the functions minB and maxB to calculate he minima and maxima for a directory of PDB files
for x in pdb_files:
    structure_name=x[-8:-4]
    file_name=x[-8:]
    index_names.append(structure_name)
    structure = parser.get_structure(structure_name,file_name)
    min_b_list.append(minB(structure))
    max_b_list.append(maxB(structure))
      
temp = {'MinValue':min_b_list,'MaxValue':max_b_list}

# Create a Pandas dataframe to store thedata in the following order: Structure name-MinB-maxB
# Write the dataframe into a .csv file
df = pd.DataFrame(temp,columns = ['MinValue','MaxValue'], index=index_names)
df.to_csv('B_Factors_Ranges.csv')
ordered_df = df.sort_values(by='MinValue')

# Plot the dataframe into a stacked dumb-bell plot
my_range=range(1,len(df.index)+1)
plt.hlines(y=my_range, xmin=ordered_df['MinValue'], xmax=ordered_df['MaxValue'],color='grey', alpha=0.4)
plt.scatter(ordered_df['MinValue'],my_range, s=300 , color='skyblue', alpha=1, label='MinBFactor')
plt.scatter(ordered_df['MaxValue'], my_range, s=300, color='green', alpha=0.4 , label='MaxBFactor')
plt.rcParams["figure.figsize"] = (50,40)

plt.xticks(fontsize=50)

ticks=list(range(1,(no_of_files+1)))
plt.yticks(ticks,[ '2pvr', '3h30','3q04','3q9w','3war', '5cu4','5csp','3nsz','4kwp', '5cvg', '3r0t', '5clp', '3pw1','3owk', '5cu6', '5csv', '3mb7', '3bqc', '4rll','5cs6'],fontsize=50)

plt.title("B-Factor Ranges",size=50)
plt.xlabel('B-Factors',size=50)
plt.show()
