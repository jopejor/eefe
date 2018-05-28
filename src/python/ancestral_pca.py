# Libraries

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy import interpolate
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import sys
from pylab import *
#import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
#sns.set_style('darkgrid')
from sklearn.cluster import KMeans
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
from sklearn import cross_validation
from sklearn.linear_model import LinearRegression
from sklearn import linear_model, decomposition, datasets
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
import matplotlib as mpl


#Read data

name="ancestral.csv"
df=pd.read_csv(name)
X = df.iloc[:,4:].values
y = df.Condition.values#Or any other factor
treatment = df.Linage.values#Or any other factor

########################### PCA 2D  ##################################

pca = PCA()
X_r = pca.fit_transform(scale(X,with_std=False))


# Create 3 different clusters
kmeans = KMeans(n_clusters=5)
kmeans.fit(X_r)
# Get the cluster centroids
centers = kmeans.cluster_centers_
centers

# 2D PCA
plt.figure()
ig, ax = plt.subplots()
colors = ['royalblue', 'sandybrown', 'firebrick']
lw = 2
conditions=unique(y)


for color, i, conditions in zip(colors, [15, 37, 43], conditions):
    plt.scatter(X_r[y == i, 0], X_r[y == i, 1], color=color, alpha=.8, lw=lw,
                label=conditions)
plt.legend(loc='best', shadow=False, scatterpoints=1)
plt.grid(linestyle='dotted')
# Overlay the centroids on the scatter plot
plt.scatter(centers[:, 0], centers[:, 1], c='black', s=100,marker='x')
ax.set_xlabel('PC1',fontsize=12, fontweight='bold')
ax.set_ylabel('PC2',fontsize=12, fontweight='bold')
fontsize = 12
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
     tick.label1.set_fontsize(fontsize)
     tick.label1.set_fontweight('bold')


labels=df.Linage.values
for label, x, y in zip(labels, X_r[:, 0], X_r[:, 1]):
	plt.annotate(
		label,
		xy=(x, y),
        textcoords='data',fontsize=7)


pp = PdfPages('ancestral_pca.pdf')
plt.savefig(pp,format='pdf')
pp.close()


########################### PCA Variance  ##################################

# Select all 64 principal components
pca = PCA(92)  # project from 64 to 2 dimensions
X_r = pca.fit_transform(scale(X,with_std=False))

# Obtain the explained variance for each principal component
varianceExp= pca.explained_variance_ratio_
# Compute the total sum of variance
totVarExp=np.cumsum(np.round(pca.explained_variance_ratio_, decimals=4)*100)

# Plot the variance explained as a function of the number of principal components
plt.plot(totVarExp)
plt.xlabel('No of principal components')
plt.ylabel('% variance explained')
plt.title('No of Principal Components vs Total Variance explained')

plt.show()

################# t-SNE


X = df.iloc[:,4:].values
y = df.Condition.values#Or any other factor

# Fit and transform with a TSNE
from sklearn.manifold import TSNE
tsne = TSNE(n_components=2, perplexity=20, n_iter=3000)
X_2d = tsne.fit_transform(X)

# Visualize the data
plt.scatter(X_2d[:, 0], X_2d[:, 1], c=y) 
plt.show()


################ Isomap
from sklearn import (manifold, datasets, decomposition, ensemble,
                     discriminant_analysis, random_projection)
n_neighbors = 32

print("Computing Isomap embedding")
X_iso = manifold.Isomap(n_neighbors, n_components=2).fit_transform(X)
print("Done.")
plt.scatter(X_iso[:, 0], X_iso[:, 1], c=y) 
plt.show()

####### LLE
from sklearn import manifold, datasets
X = df.iloc[:,4:].values
y = df.Condition.values#Or any other factor


print("Computing LLE embedding")
X_r, err = manifold.locally_linear_embedding(X, n_neighbors=32,
                                             n_components=3)
print("Done. Reconstruction error: %g" % err)

#----------------------------------------------------------------------
# Plot result

fig = plt.figure()

ax = fig.add_subplot(211, projection='3d')
ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap=plt.cm.Spectral)

ax.set_title("Original data")
ax = fig.add_subplot(212)
ax.scatter(X_r[:, 0], X_r[:, 1], c=y, cmap=plt.cm.Spectral)
plt.axis('tight')
plt.xticks([]), plt.yticks([])
plt.title('Projected data')
plt.show()