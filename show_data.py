# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
from mdtools import StatQ, Isomorph, RDFAnalysis

plt.style.use('default')
os.chdir("/home/gn/Code/MD-simulation/examples/examplebin")


N = [6, 7, 8, 9, 10, 12]
RHO = [0.3, 0.5, 1.0, 1.5]
T = [0.5, 1.0, 2.0]
A = [0, 0.25, 0.50, 0.75, 0.8, 0.90]

run = StatQ(10000, 6**3*4)
run.sf(0.5, 100, 8, 0.5)

plt.show()


#%%
