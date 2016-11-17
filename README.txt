-> The code requires openMP and trng; you must have these libraries installed
-> You may need to edit Makefile to direct the compiler to the correct locations to source the libraries (the -L commands)
--> sometimes the system is smart enough to know where the libs are and you can remove the -L commands altogether, try that first
--> you can probably get openmp from Synaptics package manager; trng files are included in the repository
--> I'm afraid that getting the code linked to libraries can be a pain, let me know if you have big problems
-> This branch of the code uses the Wright-Fisher model with variable population size

-> To compile:
make clean
make

---

-> sample io is in:
./NS4A

-> For a parallel run:
You need to set the enviroment variable OMP_NUM_THREADS
./main

-> In MATLAB, to plot trajectories:
> MCSampler_Potts_plotter()
> Tcell_traj() 
> Pop_dy()
all of these functions have optional arguements 

---

NOTES

The code takes as input files: inputs.dat, P1_target.dat, P2_target.dat, resIdx.dat, h.dat, J.dat, epitopes.dat (optional) and the files listed in epitopes.dat.

There is a lot of work that still needs to be done on this code. Both in terms of functionality and optimization.
