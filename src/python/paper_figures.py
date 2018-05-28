"""
plot.py
"""
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


########################### Bar plots Figure 1 ##################################
import seaborn as sns
sns.set_style("ticks")
import matplotlib as mpl
label_size = 30
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
mpl.rcParams['lines.linewidth'] = 30
mpl.rcParams['figure.figsize'] = 20,15



name="606_clones.csv"
df=pd.read_csv(name,index_col=0)
df.plot.bar(stacked=True,color=sns.hls_palette(8, l=.3, s=.8))
sns.despine()


axes = plt.gca()
axes.set_ylim([0,6])
plt.ylabel('Number of Mutations',fontsize=30, fontweight='bold')
plt.xlabel('',fontsize=30, fontweight='bold')
plt.grid(True, which='both')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12),
          fancybox=True, shadow=True, ncol=5,fontsize=25)


pp = PdfPages('population_barplot_606.pdf')
plt.savefig(pp,format='pdf')
pp.close()

########################### Bar plots Figure 2 ##################################

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
import seaborn as sns

sns.set_style("ticks")
import matplotlib as mpl
label_size = 20

mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['figure.figsize'] = 20,15

name="607_marga_strong_norm.csv"
df=pd.read_csv(name)
axes = plt.gca()
axes.set_ylim([0,2.2])
sns.barplot(x="strain", y="growth_rate", hue="Pass", data=df)
plt.ylabel('Relative Growth Rate',fontsize=22, fontweight='bold')
plt.xlabel('',fontsize=30, fontweight='bold')

pp = PdfPages('temperature_dependence_2.pdf')
plt.savefig(pp,format='pdf')
pp.close()
########################### PCA 3D  ##################################

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
import seaborn as sns
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

name="all.csv"
df=pd.read_csv(name)
X = df.iloc[:,4:].values
y = df.Linage.values#Or any other factor

pca = PCA()
X_reduced = pca.fit_transform(scale(X,with_std=False))
np.cumsum(np.round(pca.explained_variance_ratio_, decimals=4)*100)

# To getter a better understanding of interaction of the dimensions
# plot the first three PCA dimensions
fig = plt.figure(1, figsize=(8, 6))
ax = Axes3D(fig, elev=-150, azim=110)
X_reduced = PCA(n_components=3).fit_transform(scale(X))
ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=y,
           cmap=plt.cm.Set1, edgecolor='k', s=40)
ax.set_title("First three PCA directions")
ax.set_xlabel("1st eigenvector")
ax.w_xaxis.set_ticklabels([])
ax.set_ylabel("2nd eigenvector")
ax.w_yaxis.set_ticklabels([])
ax.set_zlabel("3rd eigenvector")
ax.w_zaxis.set_ticklabels([])
plt.show()

pca = PCA(n_components=2)
X_r = pca.fit_transform(scale(X))

########################### PCA 2D  ##################################

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

name="607_linage.csv"
df=pd.read_csv(name)
X = df.iloc[:,4:].values
y = df.Condition.values#Or any other factor
treatment = df.Treatment.values#Or any other factor

pca = PCA()
X_r = pca.fit_transform(scale(X,with_std=False))
X_ancestor=X_r[1:15,:]

from sklearn.decomposition import TruncatedSVD

svd = TruncatedSVD(n_components=2,algorithm='randomized', n_iter=5, random_state=None, tol=0.0)
X_r = svd.fit(X)
X_ancestor=X_r[1:15,:]


# Create 3 different clusters
kmeans = KMeans(n_clusters=3)
kmeans.fit(X_ancestor)
# Get the cluster centroids
centers = kmeans.cluster_centers_
centers

# 2D PCA
plt.figure()
ig, ax = plt.subplots()
colors = ['navy', 'sandybrown', 'firebrick']
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
     tick.label1.set_fontweight('bold'):

pp = PdfPages('pca_606.pdf')
plt.savefig(pp,format='pdf')
pp.close()

labels=df.Strain.values
for label, x, y in zip(labels, X_r[:, 0], X_r[:, 1]):
	plt.annotate(
		label,
		xy=(x, y), xytext=(-20, 20),
        textcoords='data', ha='right', va='bottom')

#Save score
dfscore=pd.DataFrame(X_r)
dfscore.to_csv("scores_607.csv")


###################################### Hierarchical Clustering ##############################

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
from scipy.cluster.hierarchy import dendrogram, linkage

# Generate the linkage matrix using the average method
Z = linkage(X_r, 'ward')
#Plot the dendrogram
dendrogram(Z)
plt.xlabel('Data')
plt.ylabel('Distance')
plt.suptitle('Samples clustering', fontweight='bold', fontsize=14);


g = sns.clustermap(X_r)





















# Percentage of variance explained for each components
print('explained variance ratio (first two components): %s'
      % str(pca.explained_variance_ratio_))

logistic = linear_model.LogisticRegression()
pca = decomposition.PCA()
pipe = Pipeline(steps=[('pca', pca), ('logistic', logistic)])
pca.fit_transform(scale(X))
plt.figure(1, figsize=(4, 3))
plt.clf()
plt.axes([.2, .2, .7, .7])
plt.plot(pca.explained_variance_, linewidth=2)
plt.axis('tight')
plt.xlabel('n_components')
plt.ylabel('explained_variance_')

# Prediction
n_components = [30, 50, 92]
Cs = np.logspace(-4, 4, 3)

# Parameters of pipelines can be set using ‘__’ separated parameter names:
estimator = GridSearchCV(pipe,
                         dict(pca__n_components=n_components,
                              logistic__C=Cs))
estimator.fit(X, y)

plt.axvline(estimator.best_estimator_.named_steps['pca'].n_components,
            linestyle=':', label='n_components chosen')
plt.legend(prop=dict(size=12))
plt.show()

# 2D PCA
plt.figure()
colors = ['navy', 'sandybrown', 'firebrick']
lw = 2
conditions=unique(y)


for color, i, conditions in zip(colors, [15, 37, 43], conditions):
    plt.scatter(X_r[y == i, 0], X_r[y == i, 1], color=color, alpha=.8, lw=lw,
                label=conditions)
plt.legend(loc='best', shadow=False, scatterpoints=1)
plt.grid(linestyle='dotted')


labels=df.Strain.values
for label, x, y in zip(labels, X_r[:, 0], X_r[:, 1]):
    plt.annotate(
        label,
        xy=(x, y), xytext=(-20, 20),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

plt.show()



