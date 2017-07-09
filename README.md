# MD-Simulation-Data-Analysis

Meant to be used in combination with the ![MD-Simulation](https://github.com/GiannisNikiteas/MD-simulation) in order to assist in the analysis of the properties od the MD fluid.

## Getting Started

* Set ```__pathWIN``` and ```__pathLINUX``` to the location of where the data files are stored from the MD-simulation.

## Run

* Run the ***Run.py*** file with the methods from the ***FileLoadingPlotting.py***.


## Components
The ***FileLoadingPlotting.py*** file contains the following methods to analyse the fluid:
* **Energy_plots**: Plots a the Kinetic, Potential and Total Energy against time.
* **Particle_plot**: Creates a representation of the particles' last positions in 3D.
* **Vector_field**: Creates a snapshot of the particles' last positions and velocities.
* **Radial_dist_func**: Plots the Radial Distribution Fucntion (RDF) of the fluid.
* **Vel_autocor_func**: Plots the Velocity Autocorrelation Function (VAF).
* **Mean_sqr_disp**: Plots the Mean Square Displacement (MSD).
* **Pc**: Plots the Configurational (virial) pressure against parameter ***A***.
* **Avg_pressure**: Plots the Average Configurational pressure against parameter ***A***.
* **Avg_Kin**: Plots the Average Kinetic Energy against parameter ***A***.
* **Avg_Pot**: Plots the Average Potential Energy against parameter ***A***.
* **Avg_En**: Plots the Total Energy of the fluid against parameter ***A***.
* **Diffusion_plot**: Creates a plot with the diffusion coefficients of multiple runs.
* **Potential**: Plots the Pair-Potential of the fluid.
* **Force**: Plots the Force between the particles of the fluid.
* **RDF2**: Plots the Ideal Gas RDF, without the expansion terms.
