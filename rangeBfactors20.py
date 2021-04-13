import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB.MMCIFParser import MMCIFParser
parser = MMCIFParser()
structure1 = parser.get_structure("3bqc", "3bqc.cif")
structure2 = parser.get_structure("3h30", "3h30.cif")
structure3 = parser.get_structure("3mb7", "3mb7.cif")
structure4 = parser.get_structure("3nsz", "3nsz.cif")
structure5 = parser.get_structure("3pe1", "3pe1.cif")
structure6 = parser.get_structure("3war", "3war.cif")
structure7 = parser.get_structure("4kwp", "4kwp.cif")
structure8 = parser.get_structure("5csv", "5csv.cif")
structure9 = parser.get_structure("5cu4", "5cu4.cif")
structure10 = parser.get_structure("5cu6", "5cu6.cif")
structure11 = parser.get_structure("2pvr", "2pvr.cif")
structure12 = parser.get_structure("5clp", "5clp.cif")
structure13 = parser.get_structure("5csp", "5csp.cif")
structure14 = parser.get_structure("5cvg", "5cvg.cif")
structure15 = parser.get_structure("5cvg", "5cvg.cif")
structure16 = parser.get_structure("3q9w", "3q9w.cif")
structure17 = parser.get_structure("3q04", "3q04.cif")
structure18 = parser.get_structure("5cs6", "5cs6.cif")
structure19 = parser.get_structure("3owk", "3owk.cif")
structure20 = parser.get_structure("4rll", "4rll.cif")
structureData1=[]
structureData2=[]
structureData3=[]
structureData4=[]
structureData5=[]
structureData6=[]
structureData7=[]
structureData8=[]
structureData9=[]
structureData10=[]
structureData11=[]
structureData12=[]
structureData13=[]
structureData14=[]
structureData15=[]
structureData16=[]
structureData17=[]
structureData18=[]
structureData19=[]
structureData20=[]

for model in structure1:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData1.append(atom.get_bfactor())
min_Bfactor1=min(structureData1)
max_Bfactor1=max(structureData1)

for model in structure2:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData2.append(atom.get_bfactor())
min_Bfactor2=min(structureData2)
max_Bfactor2=max(structureData2)

for model in structure3:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData3.append(atom.get_bfactor())
min_Bfactor3=min(structureData3)
max_Bfactor3=max(structureData3)

for model in structure4:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData4.append(atom.get_bfactor())
min_Bfactor4=min(structureData4)
max_Bfactor4=max(structureData4)

for model in structure5:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData5.append(atom.get_bfactor())
min_Bfactor5=min(structureData5)
max_Bfactor5=max(structureData5)

for model in structure6:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData6.append(atom.get_bfactor())
min_Bfactor6=min(structureData6)
max_Bfactor6=max(structureData6)

for model in structure7:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData7.append(atom.get_bfactor())
min_Bfactor7=min(structureData7)
max_Bfactor7=max(structureData7)
                
for model in structure8:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData8.append(atom.get_bfactor()) 
min_Bfactor8=min(structureData8)
max_Bfactor8=max(structureData8)
                
for model in structure9:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData9.append(atom.get_bfactor())
min_Bfactor9=min(structureData9)
max_Bfactor9=max(structureData9)
                
for model in structure10:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData10.append(atom.get_bfactor())
min_Bfactor10=min(structureData10)
max_Bfactor10=max(structureData10)

for model in structure11:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData11.append(atom.get_bfactor())
min_Bfactor11=min(structureData11)
max_Bfactor11=max(structureData11)

for model in structure12:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData12.append(atom.get_bfactor())
min_Bfactor12=min(structureData12)
max_Bfactor12=max(structureData12)

for model in structure13:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData13.append(atom.get_bfactor())
min_Bfactor13=min(structureData13)
max_Bfactor13=max(structureData13)

for model in structure14:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData14.append(atom.get_bfactor())
min_Bfactor14=min(structureData14)
max_Bfactor14=max(structureData14)

for model in structure15:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData15.append(atom.get_bfactor())
min_Bfactor15=min(structureData15)
max_Bfactor15=max(structureData15)

for model in structure16:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData16.append(atom.get_bfactor())
min_Bfactor16=min(structureData16)
max_Bfactor16=max(structureData16)

for model in structure17:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData17.append(atom.get_bfactor())
min_Bfactor17=min(structureData17)
max_Bfactor17=max(structureData17)
                
for model in structure18:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData18.append(atom.get_bfactor()) 
min_Bfactor18=min(structureData18)
max_Bfactor18=max(structureData18)
                
for model in structure19:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData19.append(atom.get_bfactor())
min_Bfactor19=min(structureData19)
max_Bfactor19=max(structureData19)
                
for model in structure20:
    for chain in model:
         for residue in chain:
            for atom in residue:
                structureData20.append(atom.get_bfactor())
min_Bfactor20=min(structureData20)
max_Bfactor20=max(structureData20)

temp={'MinValue':[min_Bfactor1,min_Bfactor2,min_Bfactor3,min_Bfactor4,min_Bfactor5,min_Bfactor6,min_Bfactor7,min_Bfactor8,min_Bfactor9,min_Bfactor10, min_Bfactor11,min_Bfactor12,min_Bfactor13,min_Bfactor14,min_Bfactor15,min_Bfactor16,min_Bfactor17,min_Bfactor18,min_Bfactor19,min_Bfactor20],'MaxValue':[max_Bfactor1,max_Bfactor2,max_Bfactor3,max_Bfactor4,max_Bfactor5,max_Bfactor6,max_Bfactor7,max_Bfactor8,max_Bfactor9,max_Bfactor10, max_Bfactor11,max_Bfactor12,max_Bfactor13,max_Bfactor14,max_Bfactor15,max_Bfactor16,max_Bfactor17,max_Bfactor18,max_Bfactor19,max_Bfactor20]}
df = pd.DataFrame(temp,columns = ['MinValue','MaxValue'], index=['3bqc','3h30','3mb7','3nsz','3pe1','3war','4kwp','5csv','5cu4','5cu6','2pvr','5clp', '5csp', '5cvg', '3r0t', '3q9w', '3q04', '5cs6', '3owk', '4rll'])

print (df)

ordered_df = df.sort_values(by='MinValue')
my_range=range(1,len(df.index)+1)
plt.hlines(y=my_range, xmin=ordered_df['MinValue'], xmax=ordered_df['MaxValue'],color='grey', alpha=0.4)
plt.scatter(ordered_df['MinValue'],my_range, s=800 , color='skyblue', alpha=1, label='MinBFactor')
plt.scatter(ordered_df['MaxValue'], my_range, s=800, color='green', alpha=0.4 , label='MaxBFactor')
plt.rcParams["figure.figsize"] = (40,30)
plt.xticks(fontsize=40)
ticks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
plt.yticks(ticks,ordered_df.index, fontsize=30)
plt.title("B-Factor Ranges",size=90)
plt.xlabel('B-Factors',size=40)
plt.show()
plt.savefig('BFactorRange.png')  