import numpy as np

#-- Multiple Sequence Alignment
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
align = AlignIO.read("sequences_pdb20.aln", "clustal")
#print(align[:,30:330:])
align[:,30:330:]
#--------------------------------------------------------------------

# Define ranges
start_id = 30
end_id   = 300
atoms_to_be_aligned = range(start_id, end_id + 1)
#--------------------------------------------------------------------

#-- Parse necessary files with MMCIFParser
from Bio.PDB import *
pdb1 = PDBList()
parser = MMCIFParser(QUIET = True)
pdb1.retrieve_pdb_file('3war')
pdb1.retrieve_pdb_file('3BQC')
ref_structure=parser.get_structure('3WAR','3war.cif')
sample_structure=parser.get_structure('3bqc','3bqc.cif')
#--------------------------------------------------------------------

#-- Select the first models 
ref_model    = ref_structure[0]  
sample_model = sample_structure[0]

#-- Extract backbones
ref_atoms = []
ref_ca=[]
sample_atoms = []
sample_ca=[]

ref_chain=ref_model['A']
sample_chain=sample_model['A']

#-- Access residues of the Reference Structure: 3war
for ref_res in ref_chain:
    if ref_res.get_id()[1] in atoms_to_be_aligned:
       ref_atoms.append(ref_res['CA'])
       ref_ca.append(ref_res["CA"].get_coord())
         
#-- Access residues of the Variant Structure: x    

for sample_res in sample_chain:
    if sample_res.get_id()[1] in atoms_to_be_aligned:
       sample_atoms.append(sample_res['CA'])
       sample_ca.append(sample_res["CA"].get_coord())
#-----------------------------------------------------------------

from Bio.SVDSuperimposer import SVDSuperimposer
from numpy import array, dot 
      
#-- Transform lists to numpy arrays
x=np.array(ref_ca)
x.shape
y=np.array(sample_ca)
y.shape

#--Initiate SVDSuperimposer and perform a LSQ Fitting
sup = SVDSuperimposer()   

#-- y on x 
sup.set(x,y) 
sup.run()
sup.get_rms()
#--Matrix transformations: rotate and translate y on x
rot, tran = sup.get_rotran()
y_on_x1 = dot(y, rot)+tran
y_on_x2 = sup.get_transformed()


#-- NumPy operations followed by calculation of residue-wise RMSD 
x_ci=x[:,0]
y_ci=x[:,1]
z_ci=x[:,2]

x_di=y_on_x1[:,0]
y_di=y_on_x1[:,1]
z_di=y_on_x1[:,2]

dx = np.subtract(x_ci, x_di) 
dy = np.subtract(y_ci, x_di) 
dz = np.subtract(z_ci, x_di) 

dx2=np.square(dx)
dy2=np.square(dy)
dz2=np.square(dz)

d_sum=np.add(dx2,dy2,dz2)
length=end_id-start_id
d_i= np.divide(d_sum,length)
#--------------------------------------------------------------

#-- Per residue RMSD
per_residue_rmsd=np.sqrt(d_i)

#-- Average RMSD
avg_rms=np.sum(per_residue_rmsd)/length
#--------------------------------------------------------------

#--Plot
import matplotlib.pyplot as plt
plt.figure(figsize=(70,30))
plt.plot(per_residue_rmsd,color='blue',marker='o',markerfacecolor='red',markersize=5)
plt.savefig("test.tiff",dpi=600)

#--3D Plot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize = (100, 100))
ax = fig.add_subplot(444, projection='3d')
ax.scatter(x_ci,y_ci,z_ci,s=100, color='orange')
ax.scatter(x_di,y_di,z_di,s=100)

plt.show()































