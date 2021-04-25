#Itertools iterates over data structures by stepping over a for-loop usage
from itertools import chain 
from itertools import groupby
from itertools import islice

#Numpy provides a multidimensional array object, 
#various derived objects such as masked arrays and matrices), 
#and an assortment of routines for fast operations on arrays.

import numpy as np
from statistics import stdev

#MatplotLib creates static, animated, and interactive visualizations
import matplotlib.pyplot as plt

#Biopython makes Python for bioinformatics easy
#by creating high-quality, reusable modules and classes.
from Bio.PDB.MMCIFParser import MMCIFParser

#-------------------------------------------------------
# Create an object of the mmCIF file(s)
parser = MMCIFParser()
structure = parser.get_structure("3bqc", "3bqc.cif")

# Declare lists to store the following:
# 1. bfac:  Extract the list of atomic B-factors(n=3095)
# 2. segid: Extract the list of SEGIDs
# 3. indices: Returns the list of indices for atom ID="N"
# 4. modindices: Returns the list of corrected indices, starting at 1
# 5. new: Returns the list of corrected indices pruning the start value
# 6. output: Splits the list 'bfac' into sublists at breaking points taken from the list 'new'
# 7. r: generate the average values of each sublist of the list 'output'
# 8. a_mean: Assumened mean value of the list 'bfac'
# 9. std_bfac: Standard Deviation of the list 'bfac'
# 10. n_bfac: Normalized B-Factors

segid=[]
bfac=[]
modindices=[]
std_bfac=[]
n_bfac=[]

#SMCRA loop over the mmCIF file to fetch atomic B-factors and SEGIDs
for model in structure:
    for chain in model:
          for residue in chain:
               for atom in residue:
                   bfac.append(atom.get_bfactor())
                   segid.append(atom.get_name())

# Generate a list of indices for the occurence of SEGID="N"                             
indices = [i for i, x in enumerate(segid) if x == "N"]

# Modified indices
for i in indices:
    i=i+1
    modindices.append(i)
new= modindices[1:]

#Split the list 'bfac' according to the indices generated in the list 'new'
output=np.split(bfac,new) 

# Calculate average B-factor of atoms in each Residues of the crystal structure
r=list(map(np.mean,output))

# Generate Plots of Average B-factors

plt.figure(figsize=(60,30))
plt.title("Average B-Factors of Atoms in Residues", fontsize=90)
plt.plot(r,color='black',marker='o',markerfacecolor='red',markersize=2,label="Average   B-Factors")
plt.xlabel("Residues",fontsize=60)
plt.ylabel("Average B-Factors",fontsize=60)
plt.xticks(np.arange(0, 350,5), rotation='vertical',fontsize=30)
plt.yticks(np.arange(0, 80,5),fontsize=30)
plt.legend(prop={"size":50})
plt.show()
plt.savefig("Bfac_figure.tiff",dpi=600)

# Normalization of B-factors
a_mean=sum(bfac)/len(bfac)

# Standard Deviation
std_bfac=stdev(bfac)

# Normalized B-factors
for i in r:
    temp=(i-a_mean)/std_bfac
    n_bfac.append(temp)
