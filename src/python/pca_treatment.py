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



########################### PCA 2D  ##################################

name="all.csv"
df=pd.read_csv(name)
X = df.iloc[:,4:].values
y = df.Treatment.values#Or any other factor

pca = PCA()
X_r = pca.fit_transform(scale(X,with_std=False))
X_ancestor=X_r[1:18,:]



# Create 3 different clusters
kmeans = KMeans(n_clusters=3)
kmeans.fit(X_ancestor)
# Get the cluster centroids
centers = kmeans.cluster_centers_
centers

# 2D PCA
plt.figure()
ig, ax = plt.subplots()
colors = ['sandybrown', 'blueviolet', 'khaki','lightcoral']
lw = 2
conditions=unique(y)


for color, i, conditions in zip(colors, ['ancestor', 'fast', 'slow','random'], conditions):
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

pp = PdfPages('treat_pca_all.pdf')
plt.savefig(pp,format='pdf')
pp.close()



########################### PCA 2D  ##################################

name="all.csv"
df=pd.read_csv(name)
X = df.iloc[:,4:].values
y = df.Linage.values#Or any other factor

pca = PCA()
X_r = pca.fit_transform(scale(X,with_std=False))
X_ancestor=X_r[1:35,:]



# Create 3 different clusters
kmeans = KMeans(n_clusters=3)
kmeans.fit(X_ancestor)
# Get the cluster centroids
centers = kmeans.cluster_centers_
centers

# 2D PCA
plt.figure()
ig, ax = plt.subplots()
colors = ['green', 'red']
lw = 2
conditions=unique(y)


for color, i, conditions in zip(colors, [606, 607], conditions):
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

pp = PdfPages('linage_pca_all.pdf')
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
plt.xlabel('No of principal components',fontsize=12, fontweight='bold')
plt.ylabel('% variance explained',fontsize=12, fontweight='bold')
fontsize = 12
ax = gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    tick.label1.set_fontweight('bold')
for tick in ax.yaxis.get_major_ticks():
     tick.label1.set_fontsize(fontsize)
     tick.label1.set_fontweight('bold')

pp = PdfPages('var_606.pdf')
plt.savefig(pp,format='pdf')
pp.close()


################# t-SNE


X = df.iloc[:,4:].values


# Fit and transform with a TSNE
from sklearn.manifold import TSNE
tsne = TSNE(n_components=3, perplexity=15, n_iter=6000)
X_2d = tsne.fit_transform(X)

conditions=unique(y)
fig = plt.figure(1, figsize=(8, 6))
ax = Axes3D(fig, elev=-150, azim=110)
for color, i, conditions in zip(colors, [606, 607], conditions):
    ax.scatter(X_2d[y == i, 0], X_2d[y == i, 1],  X_2d[y == i, 2], color=color, alpha=.8,
        label=conditions)
plt.show()

pp = PdfPages('tSNE_606.pdf')
plt.savefig(pp,format='pdf')
pp.close()


################ Isomap
from sklearn import (manifold, datasets, decomposition, ensemble,
                     discriminant_analysis, random_projection)
n_neighbors = 12

print("Computing Isomap embedding")
X_iso = manifold.Isomap(n_neighbors, n_components=3).fit_transform(X)
print("Done.")

fig = plt.figure(1, figsize=(8, 6))
ax = Axes3D(fig, elev=-150, azim=110)
conditions=unique(y)
for color, i, conditions in zip(colors, [606,607], conditions):
    ax.scatter(X_iso[y == i, 0], X_iso[y == i, 1],  X_iso[y == i, 2], color=color, alpha=.8, lw=lw,
        label=conditions)


plt.show()

pp = PdfPages('iso_606.pdf')
plt.savefig(pp,format='pdf')
pp.close()

####### LLE
from sklearn import manifold, datasets
X = df.iloc[:,4:].values



print("Computing LLE embedding")
X_r, err = manifold.locally_linear_embedding(X, n_neighbors=32,
                                             n_components=3)
print("Done. Reconstruction error: %g" % err)


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

pp = PdfPages('LLE_606.pdf')
plt.savefig(pp,format='pdf')
pp.close()



def my_kernel(X, Y):
    """
    We create a custom kernel:

                 (2  0)
    k(X, Y) = X  (    ) Y.T
                 (0  1)
    """
    M = np.array([[2, 0], [0, 1.0]])
    return np.dot(np.dot(X, M), Y.T)


h = .02  # step size in the mesh

# we create an instance of SVM and fit out data.
clf = svm.SVC(kernel=my_kernel)
clf.fit(X, Y)

# Plot the decision boundary. For that, we will assign a color to each
# point in the mesh [x_min, x_max]x[y_min, y_max].
x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])

# Put the result into a color plot
Z = Z.reshape(xx.shape)
plt.pcolormesh(xx, yy, Z, cmap=plt.cm.Paired)

# Plot also the training points
plt.scatter(X[:, 0], X[:, 1], c=Y, cmap=plt.cm.Paired, edgecolors='k')
plt.title('3-Class classification using Support Vector Machine with custom'
          ' kernel')
plt.axis('tight')
plt.show()