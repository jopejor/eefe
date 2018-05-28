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


if __name__ == "__main__":

	strain=sys.argv[1]

	df=pd.read_csv(strain,index_col=0)
	trans_df=df.transpose()
	cmap = cmx.get_cmap('Spectral')


	trans_df.plot(linewidth=2,cmap=cmap);
	plt.xlabel('Time',fontsize=18)
	plt.ylabel('Mutation Frequency', fontsize=18)



	filname=strain+".pdf"

	pp = PdfPages(filname)
	plt.savefig(pp,format='pdf')
	pp.close()
	


