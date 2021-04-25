# Read More: https://jakevdp.github.io/PythonDataScienceHandbook/05.09-principal-component-analysis.html
# Principal component analysis is a fast and flexible unsupervised method for dimensionality reduction in data
import pandas as pd
import numpy as np
from pandas import DataFrame
from sklearn.decomposition import PCA
from numpy import linalg as LA
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

# The dataset used is of size (270 X 19) and comprises the relative C-alpha backbone distances of different PDB structures of the human CK2 protein to a reference CK2 protein
# Load the prepared dataset
df = pd.read_csv(r"C:\Users\saySa\OneDrive\Desktop\RMSD_data.csv")

#Standardise the dataset so that mean for each Feature=0 and Standard Deviation =1
df_std= (df-df.mean()) / (df.std()) 

#Coviariance matrix of the dataset (19 X 19)
k=df_std.cov()

# The numpy.linalg.eig function returns a tuple consisting of a vector and an array. 
# The vector (here w) contains the eigenvalues. The array (here v) contains the corresponding eigenvectors, one eigenvector per column. 
# The eigenvectors are normalized so their Euclidean norms are 1. 
w,v = LA.eig(k)

# Pick k eigenvalues from a matrix of eigenvectors
u = v[:,0:2]

# [Features_matrix] (270*19)  X  [Top_Eigenvector_Matrix] (19*2)  =  [Transformed_matrix] (270*2)
transformed=-(df_std.dot(u)) 
transformed2=transformed.to_numpy()
# Target set
y = df_std['3bqc (30-300)']

# We can now plot the first two principal components of each point to learn about the data
plt.scatter(transformed2[:, 0], transformed2[:, 1],  c=y, edgecolor='none', alpha=0.7,cmap=plt.cm.get_cmap('nipy_spectral', 10))
plt.xlabel('Component 1')
plt.ylabel('Component 2')
plt.colorbar();
