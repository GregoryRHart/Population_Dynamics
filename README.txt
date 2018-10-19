-> The code requires openMP and trng; you must have these libraries installed
-> You may need to edit Makefile to direct the compiler to the correct locations to source the libraries (the -L commands)
--> sometimes the system is smart enough to know where the libs are and you can remove the -L commands altogether, try that first
--> you can probably get openmp from Synaptics package manager; trng files are included in the repository
--> I'm afraid that getting the code linked to libraries can be a pain.

-> There are two main files (main.cpp.ode and main.cpp.gilps). For the ODE version of the simulator copy main.cpp.ode to main.cpp before compiling. For the stochastic version copy main.cpp.gilps to main.cpp before compiling.
-> To compile:
make clean
make

---

-> For a parallel run: You need to set the enviroment variable OMP_NUM_THREADS
./main

-> In MATLAB, to plot trajectories for individual runs:
> MCSampler_Potts_plotter()
> Tcell_traj() 
> Pop_dy()
all of these functions have optional arguements 

-> Rest of the *.m files are for analysising groups of runs.
---

NOTES

The code takes as input files: inputs.dat, P1_target.dat, P2_target.dat, resIdx.dat, h.dat, J.dat, epitopes.dat (optional) and the files listed in epitopes.dat.
