import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
from Bio.PDB.MMCIFParser import MMCIFParser
parser = MMCIFParser()

pdb_files = glob(r'C:\Users\saySa\OneDrive\Desktop\pdbfiles\*')
no_of_files=len(pdb_files)

def minB(structure):
    b_factor=[]
    for model in structure:
          for chain in model:
              for residue in chain:
                  for atom in residue:
                      b_factor.append(atom.get_bfactor())
    return(min(b_factor))
    
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

for x in pdb_files:
    structure_name=x[-8:-4]
    file_name=x[-8:]
    index_names.append(structure_name)
    structure = parser.get_structure(structure_name,file_name)
    min_b_list.append(minB(structure))
    max_b_list.append(maxB(structure))
      
temp = {'MinValue':min_b_list,'MaxValue':max_b_list}
df = pd.DataFrame(temp,columns = ['MinValue','MaxValue'], index=index_names)
df.to_csv('B_Factors_Ranges.csv')

ordered_df = df.sort_values(by='MinValue')
my_range=range(1,len(df.index)+1)
plt.hlines(y=my_range, xmin=ordered_df['MinValue'], xmax=ordered_df['MaxValue'],color='grey', alpha=0.4)
plt.scatter(ordered_df['MinValue'],my_range, s=800 , color='skyblue', alpha=1, label='MinBFactor')
plt.scatter(ordered_df['MaxValue'], my_range, s=800, color='green', alpha=0.4 , label='MaxBFactor')
plt.rcParams["figure.figsize"] = (40,30)
plt.xticks(fontsize=40)
ticks=list(range(no_of_files))
plt.yticks(ticks,ordered_df.index, fontsize=30)
plt.title("B-Factor Ranges",size=90)
plt.xlabel('B-Factors',size=40)
plt.show()
