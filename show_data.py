# %%
import os
import matplotlib.pyplot as plt
from mdtools import StatQ, Isomorph, RDFAnalysis, ParticleVisualisation

plt.style.use('default')
# %%
# NOTE: change to your directory where the data is
os.chdir("/home/gn/Code/md-sim/examples/examplebin")

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

# plt.show()
# %%

plt.style.use('default')
os.chdir("/home/gn/Code/md-sim/examples/examplebin")

# Plot the fluid visualisation methods

# Lennard Jones pair potential
run = ParticleVisualisation(5000, 10**3)
run.particle_plot("lj_simple_run_", 0.5, 0.5)
run.vector_field("lj_simple_run_", 0.5, 0.5)
run.vector_field_3d("lj_simple_run_", 0.5, 0.5)
plt.show()

# %%

plt.style.use('default')
os.chdir("/home/gn/Code/md-sim/examples/examplebin")

# Plot the 3D data
run = ParticleVisualisation(1000, 5**3)
run.animation3D("3D_view_", 0.5, 0.5, 8, 0.5)
# %%
