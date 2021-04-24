import pandas as pd
import numpy as np
from pandas import DataFrame
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

dist_data = pd.read_csv(r"C:\Users\saySa\OneDrive\Desktop\RMSD_data.csv")
dist_data.shape

pca = PCA(2)  # project from 64 to 2 dimensions
projected = pca.fit_transform(dist_data)
print(dist_data.shape)
print(projected.shape)

y = dist_data['3bqc (30-300)']


plt.scatter(projected[:, 0], projected[:, 1],  c=y, edgecolor='none', alpha=0.7,cmap=plt.cm.get_cmap('nipy_spectral', 10))
plt.xlabel('component 1')
plt.ylabel('component 2')
plt.colorbar();

