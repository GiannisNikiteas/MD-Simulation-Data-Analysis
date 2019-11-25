# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
from mdtools import StatQ, Isomorph, RDFAnalysis

plt.style.use('default')
os.chdir("/home/gn/Code/MD-simulation/examples/examplebin")

# Lennard Jones pair potential
run = StatQ(5000, 10**3)
run.rdf_plot("lj_simple_run_", 0.5, 0.5)
run.msd("lj_simple_run_", 0.5, 0.5)
run.vaf("lj_simple_run_", 0.5, 0.5)

# Exponential pair potential
run = StatQ(5000, 10**3)
run.rdf_plot("exp_simple_run_", 0.5, 0.5, 8, 0.5)
run.msd("exp_simple_run_", 0.5, 0.5, 8, 0.5)
run.vaf("exp_simple_run_", 0.5, 0.5, 8, 0.5)

# GCM pair potential
run = StatQ(5000, 10**3)
run.rdf_plot("gcm_simple_run_", 0.5, 0.5)
run.msd("gcm_simple_run_", 0.5, 0.5)
run.vaf("gcm_simple_run_", 0.5, 0.5)

# BIP pair potential
run = StatQ(5000, 10**3)
run.rdf_plot("bip_simple_run_", 0.5, 0.5, 8, 0.5)
run.msd("bip_simple_run_", 0.5, 0.5, 8, 0.5)
run.vaf("bip_simple_run_", 0.5, 0.5, 8, 0.5)

plt.show()


#%%
