"""
Model for the temperature dependent growth rate
"""
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as nprandom
import scipy.stats as sps
import copy
import random as random
import csv
from scipy import integrate
import math



# Constants

m=0.3
l=200
c=30

start=0
end=1000
numsteps=1000
time=np.linspace(start,end,numsteps)
y0=np.array([1])

#Functions

def deriv(y,t,m,l,c):
	yprime = np.array([m*math.exp((m*(math.exp(1)*l - math.exp(1)*t))/c-math.exp((m*(math.exp(1)*l - math.exp(1)*t))/c + 1) + 2)+0*y[0]])
	return yprime

# Calculations

y=integrate.odeint(deriv,y0,time,(m,l,c))

plt.plot(time,y[:])
plt.show()
