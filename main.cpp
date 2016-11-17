/*
 *  main.cpp
 *  Population_Dynamics
 *
 *  Created by Greg Hart on May 2015
 *  Copyright 2015 ALF. All rights reserved.
 *
 */
/*
 * IN
 * resIdx.dat 		- ordered list of residue types - as integer codes - present at each site from most to least probable
 * h.dat 		- h parameters
 * J.dat 		- J parameters
 * init_seq.dat		- (optional) a starting sequence, if file is not found then the WT is used.
 *
 * OUT
 * P1_model_traj.dat 	- running trajectory of P1 values of the population
 * P2_model_traj.dat 	- running trajectory of P2 values of the population
 * P1_model.dat 	- terminal P1 values computed by the population path
 * P2_model.dat	 	- terminal P2 values computed by the population path
 * MC_seqs.dat 		- sequences sampled from population trajectory
 * pop_stats.dat	- data on the population at different time points
 * Tcell_traj.dat	- Number of Tcells and targets at different time points
 */
 
#include "main.h"
#include "functions.h"
#include "d_functions.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <limits>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
//#include <stdint.h>
//#include <thrust/sort.h>

#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/poisson_dist.hpp>
#include <boost/numeric/odeint.hpp>

#define PARTPARALLEL 1

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::vector;
using std::sort;
using namespace boost::numeric::odeint;


typedef std::vector<double> state_type;
const double a = 1e-7;
const double ap = 1e-6;
const double b = 1e-4;
const double d = 0.2;//0.5;
const double dp = 3;
const double re = 6;//4;
const double e = 3*1e-4;
const double g = 0.03;
const double hr = 0.03;
const double w = .1;//7.5;
const double aff_rate = 1e-2;
double Ichi = 0;
int rep_lim = 0;

// integrator if there is only one generation of effector cells
void int_Tcell_1(const state_type &x, state_type &dxdt, double t){
    dxdt[0] = -a*x[0]*Ichi + b*x[0] - e*x[0] + w;
    dxdt[1] = (a*x[0]+ap*x[2])*Ichi - dp*x[1] - g*x[1];
    dxdt[2] = -ap*x[2]*Ichi - hr*x[2] + g*x[1];
}

// integrator if there are only two generations of effector cells
void int_Tcell_2(const state_type &x, state_type &dxdt, double t){
    dxdt[0] = -a*x[0]*Ichi + b*x[0] - e*x[0] + w;
    dxdt[1] = (a*x[0]+ap*x[3])*Ichi - re*x[1] - d*x[1] - g*x[1];
    dxdt[2] =  2*re*x[1] - dp*x[2] - g*x[2];
    dxdt[3] = -ap*x[3]*Ichi - hr*x[3] + g*(x[1]+x[2]);
}

// integrator if there three or more generations of effector cells
void int_Tcell_N(const state_type &x, state_type &dxdt, double t){
    double Etotal = 0.0;
    dxdt[0] = -a*x[0]*Ichi + b*x[0] - e*x[0] + w;
    dxdt[1] = (a*x[0]+ap*x[rep_lim+1])*Ichi - re*x[1] - d*x[1] - g*x[1];
    Etotal += x[1];
    for(int i=2; i<rep_lim; i++){
        dxdt[i] = -re*x[i] + 2*re*x[i-1] - d*x[i] - g*x[i];
        Etotal += x[i];
    }
    dxdt[rep_lim] =  2*re*x[rep_lim-1] - dp*x[rep_lim] - g*x[rep_lim];
    Etotal += x[rep_lim];
    dxdt[rep_lim+1] = -ap*x[rep_lim+1]*Ichi - hr*x[rep_lim+1] + g*(Etotal);
}


// main function preforming Fisher-wright dynamics on the viral side and intergrating ODEs for the Tcell side
int main (int argc, char *argv[]) {
	
	cout.setf(std::ios_base::scientific);
	cout.precision(3);	
	
	#ifdef UNIT_TESTING
	
	#endif  // #ifdef UNIT_TESTING

	
	// start timer
    double wall_start = get_wall_time();
    double cpu_start  = get_cpu_time();
	
	
	// inputs with default values
    long n_inputs=15;
	
    long seed = -240164;
    long m=3;
    long N=10000;
    double rate=0.0001;
    long progeny = 10;
    long n_cycles = 20000;
    long burnin = 2000;
    long sample_mod = 10000;
    double T = 1;
    long print_mod = 1;
    long write_mod = 5;
    long n_epitope = 1;
    double n_T = 1000;
    rep_lim = 9;
    double T_penalty = 1;
    
    std::ifstream fin_inputs;
    fin_inputs.open("./inputs.dat");
    if (!fin_inputs) {
	cerr << "Cannot open input file inputs.dat in the current directory; aborting." << endl;
	exit(-1);
    }
    long input_cntr=0;
    while (!fin_inputs.eof()) {
			
	std::string input_name, tmp_str;
	fin_inputs >> input_name;
			
	if (strcmp(input_name.c_str(),"seed")==0) {
	    fin_inputs >> seed;
	    input_cntr++;
	} else if (strcmp(input_name.c_str(),"m")==0) {
	    fin_inputs >> m;
	    input_cntr++;
	} else if (strcmp(input_name.c_str(),"N")==0) {
	    fin_inputs >> N;
	    input_cntr++;
	} else if (strcmp(input_name.c_str(),"rate")==0) {
	    fin_inputs >> rate;
	    input_cntr++;
        } else if (strcmp(input_name.c_str(),"progeny")==0) {
            fin_inputs >> progeny;
            input_cntr++;
	} else if (strcmp(input_name.c_str(),"n_cycles")==0) {
	    fin_inputs >> n_cycles;
	    input_cntr++;
	} else if (strcmp(input_name.c_str(),"burnin")==0) {
	    fin_inputs >> burnin;
	    input_cntr++;
	} else if (strcmp(input_name.c_str(),"sample_mod")==0) {
	    fin_inputs >> sample_mod;
	    input_cntr++;
	} else if (strcmp(input_name.c_str(),"T")==0) {
	    fin_inputs >> T;
	    input_cntr++;
	} else if (strcmp(input_name.c_str(),"print_mod")==0) {
	    fin_inputs >> print_mod;
	    input_cntr++;
	} else if (strcmp(input_name.c_str(),"write_mod")==0) {
	    fin_inputs >> write_mod;
	    input_cntr++;
        } else if (strcmp(input_name.c_str(),"n_epitope")==0) {
            fin_inputs >> n_epitope;
            input_cntr++;
        } else if (strcmp(input_name.c_str(),"n_T")==0) {
            fin_inputs >> n_T;
            input_cntr++;
        } else if (strcmp(input_name.c_str(),"rep_lim")==0) {
            fin_inputs >> rep_lim;
            input_cntr++;
        } else if (strcmp(input_name.c_str(),"T_penalty")==0) {
            fin_inputs >> T_penalty;
            input_cntr++;
	}
	getline(fin_inputs,tmp_str);
    }
	
    if (input_cntr!=n_inputs) {
	cerr << "ERROR - did not read expected number of inputs from inputs.dat; aborting." << endl;
	exit(-1);
    }
	
    fin_inputs.close();
	
    // loading parameters
    //- loading ordered residue indices at each site (resIdx.dat), computing total number of h and J parameters, 
    //  and candidateSites list of sites containing more than one residue type as only sites at which we attempt mutations
    std::vector<long> nRes(m,-1);
    std::vector< std::vector<long> > resIdx(m);
    load_resIdx("./resIdx.dat",nRes,resIdx,m);
    long n_h=0;
    long n_J=0;
    for(long i=0; i<m; i++) {
	n_h += nRes[i];
        for(long j=i+1; j<m; j++) {
		n_J += nRes[i]*nRes[j];
	}
    }
   
    // 6 hours to replicate, 7 day half life, 40 replication complexes per cell.
//    long progeny = 4*7*40*(n_h-m)/(19*m);
    cout << "progeny = " << progeny << endl;
 
    std::vector<long> epi_start(n_epitope,-1);
    std::vector<long> epi_end(n_epitope,-1);
    std::vector< std::vector<double> > chi(n_epitope);
    if(n_epitope > 0){
        load_epitopes("./epitopes.dat", epi_start, epi_end, chi, n_epitope, resIdx, nRes);
    }
    
    cout << "\n";
    cout << "$$ # parameters: $$" << endl;
    cout << "  # h parameters = " << n_h << endl;
    cout << "  # J parameters = " << n_J << endl;
    cout << "  # epitopes = " << n_epitope << endl;
    cout << endl;
	
    std::vector <long> candidateSites(0);
    for(long i=0; i<m; i++) {
	if (nRes[i]>1) {
	    candidateSites.push_back(i);
	}
    }
    long n_candidateSites = candidateSites.size();
	
		
    //- loading h parameters
    std::vector< std::vector<double> > h(m);
    for(long i=0; i<m; i++) {
	(h[i]).resize(nRes[i],(double)0);
    }
    load_X1("./h.dat",nRes,h,m);
	
	
    //- loading J parameters
    // \-> loading only upper triangle 
    // \-> loading residuewise elements in column major order (i.e., over i residues first)
    std::vector< std::vector< std::vector< std::vector<double> > > > J(m);
    for(long i=0; i<m; i++) {
	(J[i]).resize(m);
        for(long j=i+1; j<m; j++) {
	    (J[i][j]).resize(nRes[i]);
            for(long p=0; p<nRes[i]; p++) {
		(J[i][j][p]).resize(nRes[j],(double)0);
	    }
	}
    }
    load_X2("./J.dat",nRes,J,m);
	
	
      // initializing trajectory files
    FILE* fout_P1_traj;
    {
        std::string fstr_P1_traj="./P1_model_traj.dat";
        fout_P1_traj = fopen(fstr_P1_traj.c_str(),"wb");
        if (fout_P1_traj==NULL) {
	    cerr << "Cannot open output file " << fstr_P1_traj << " in the current directory; aborting." << endl;
	    exit(-1);
        }
        int32_t sizes [2];
        sizes[0] = sizeof(long);
        sizes[1] = sizeof(double);
        fwrite(sizes,4,2,fout_P1_traj);    
    } 
	
    FILE* fout_P2_traj;
    {
        std::string fstr_P2_traj="./P2_model_traj.dat";
        fout_P2_traj = fopen(fstr_P2_traj.c_str(),"wb");
        if (fout_P2_traj==NULL) {
	    cerr << "Cannot open output file " << fstr_P2_traj << " in the current directory; aborting." << endl;
	    exit(-1);
	}
        int32_t sizes [2];
        sizes[0] = sizeof(long);
        sizes[1] = sizeof(double);
        fwrite(sizes,4,2,fout_P2_traj);    
   }

   FILE* fout_MC_seqs;
   {
        std::string fstr_MC_seqs="./MC_seqs.dat";
        fout_MC_seqs = fopen(fstr_MC_seqs.c_str(),"wb");
        if (fout_MC_seqs==NULL) {
	    cerr << "Cannot open output file " << fstr_MC_seqs << " in the current directory; aborting." << endl;
	    exit(-1);
        }
        int32_t sizes [2];
        sizes[0] = sizeof(long);
        sizes[1] = sizeof(int8_t);
        fwrite(sizes,4,2,fout_MC_seqs);
        fwrite(&N,sizeof(long),1,fout_MC_seqs); 
    }
    
    
    FILE* fout_pop_stats;
    {
        std::string fstr_pop_stats="./pop_stats.dat";
        fout_pop_stats = fopen(fstr_pop_stats.c_str(),"wb");
        if (fout_pop_stats==NULL) {
            cerr << "Cannot open output file " << fstr_pop_stats << " in the current directory; aborting." << endl;
            exit(-1);
        }
        int32_t sizes [2];
        sizes[0] = sizeof(long);
        sizes[1] = sizeof(double);
        fwrite(sizes,4,2,fout_pop_stats);    
    }
    

    FILE* fout_Tcell_traj;
    {
        std::string fstr_Tcell_traj="./Tcell_traj.dat";
        fout_Tcell_traj = fopen(fstr_Tcell_traj.c_str(),"wb");
        if (fout_Tcell_traj==NULL) {
            cerr << "Cannot open output file " << fstr_Tcell_traj << " in the current directory; aborting." << endl;
            exit(-1);
        }
        int32_t sizes [3];
        sizes[0] = sizeof(long);
        sizes[1] = sizeof(double);
        sizes[2] = n_epitope;
        fwrite(sizes,4,3,fout_Tcell_traj);    
    }

	
    // initializing model probability arrays
    std::vector< std::vector<double> > n1(h);
    std::vector< std::vector< std::vector< std::vector<double> > > > n2(J);
    for(long i=0; i<m; i++) {
        for(long p=0; p<nRes[i]; p++) {
	    n1[i][p] = 0;
            for(long j=i+1; j<m; j++) {
                for(long q=0; q<nRes[j]; q++) {
		    n2[i][j][p][q] = 0;
		}
	    }
	}
    }
    
	printf("PARTPARALLEL = %d\n", PARTPARALLEL);
	
    // Population dynamics over empirical fitness landscape
    cout << "Commencing evolutionary dynamics on the Potts model..." << endl;
	
    // initialize population arrays
    std::vector< std::vector<int8_t> > population(N, std::vector<int8_t>(m));
    std::vector< std::vector<int8_t> > temp_pop(progeny*N, std::vector<int8_t>(m));
    double *boltzfac = new double[N];
    vector<double> temp_boltz(N*progeny);
    double *part_sum = new double[N*progeny];
	bool *part_sum_used = new bool[N*progeny];
    double Z = 1.0;
    double Ztemp = 1.0;
    double *boltzfac_eff = new double[N];
    double Z_eff = 1.0;
    int *pop_init = new int[m];
    bool *occupied = new bool[N*progeny];
    int *idx = new int[N*progeny];
	
    // Initialize every sequence in the population to the sequence found in the file or the WT if no file.
    load_seq("init_seq.dat", pop_init, nRes, resIdx, m);
    for(long i=0; i<N; i++){
        for(long j=0; j<m; j++){
            population[i][j]=pop_init[j];
	}
    }
	
	// initialize occupied arrays
	for(long i = 0; i < N*progeny; i++)
		occupied[i] = false;
	
    //setup T cell populations
    std::vector<state_type> Tcells(n_epitope, state_type(rep_lim+2));
    double chiI[n_epitope];
    std::vector< std::vector<int> > place_value(n_epitope);
    std::vector< std::vector<int> > n_WTepitopes(n_epitope);
    for (long i=0; i<n_epitope; i++){
        (n_WTepitopes[i]).resize((chi[i]).size(), 0);
	// create a way to index all mutations in an epitope
        (place_value[i]).resize(epi_end[i]-epi_start[i], 0);
        int size = 1;
        int count = 0;
	for(int j=epi_start[i]; j<epi_end[i]; j++){
	    place_value[i][count] = size;
	    size *= nRes[j];
	    count++;
	}
        Tcells[i][0] = 200.0;//500;
        Tcells[i][rep_lim+1] = 0.0;
        for(long j=1; j<=rep_lim; j++){
            Tcells[i][j] = 0.0;
        }
    }

    double temp = 0;
    for(long i=0; i<m; i++) {
        for(long p=0; p<nRes[i]; p++) {
            temp += h[i][p];
        }
    }
    T_penalty = T_penalty*temp/n_h;
    cout << "T penalty: " << T_penalty << endl;

    // Calulate the boltzmann factor for each sequence and the partition function as if only the sequences in the population are possible.
    boltzfac[0] = exp(-energy(population[0],m,h,J)/T);
    boltzfac_eff[0] = boltzfac[0];
    for(long i=1; i<N; i++){
        boltzfac[i]=boltzfac[0];
        boltzfac_eff[i] = boltzfac[0];
    }
    Z = boltzfac[0];
    Ztemp = Z;
    Z_eff = boltzfac_eff[0];
    
    // calculate initial refrencies of amino acids
    for(long k=0; k<N; k++){
        for(long i=0; i<m; i++) {
             n1[i][population[k][i]] += 1.0;
             for(long j=i+1; j<m; j++) {
                n2[i][j][population[k][i]][population[k][j]] += 1.0;
             }
        }
    }
    
    // write initial state to file
    long zero = 0;
    fwrite(&zero, sizeof(long), 1, fout_pop_stats);
    fwrite(&N, sizeof(long), 1, fout_pop_stats);
    fwrite(&Z, sizeof(double), 1, fout_pop_stats);
    fwrite(&Ztemp, sizeof(double), 1, fout_pop_stats);
    fwrite(&Z_eff, sizeof(double), 1, fout_pop_stats);
    fwrite(&zero, sizeof(long), 1, fout_P1_traj);
    for(long i=0; i<m; i++) {
        for(long p=0; p<nRes[i]; p++) {
            n1[i][p] = n1[i][p]/N;
            fwrite(&n1[i][p], sizeof(double), 1, fout_P1_traj);
            n1[i][p] = n1[i][p]*N;
        }
    }
                 
    fwrite(&zero, sizeof(long), 1, fout_P2_traj);
    for(long i=0; i<m; i++) {
        for(long j=i+1; j<m; j++) {
            for(long q=0; q<nRes[j]; q++) {
                for(long p=0; p<nRes[i]; p++) {
                    n2[i][j][p][q] = n2[i][j][p][q]/N;
                    fwrite(&n2[i][j][p][q], sizeof(double), 1, fout_P2_traj);
                    n2[i][j][p][q] = n2[i][j][p][q]*N;
                }
            }
        }
    }

    // reset amino acid frequencies
    for(long i=0; i<m; i++) {
        for(long p=0; p<nRes[i]; p++) {
            n1[i][p] = 0.0;
            for(long j=i+1; j<m; j++) {
                for(long q=0; q<nRes[j]; q++) {
                    n2[i][j][p][q] = 0.0;
                }
            }
        }
    }

    fwrite(&zero, sizeof(long), 1, fout_Tcell_traj);
    for(long i=0; i<n_epitope; i++){
        double Etotal = 0.0;
        fwrite(&Tcells[i][0], sizeof(double), 1, fout_Tcell_traj);
        for(long j=1; j<=rep_lim; j++){
            Etotal += Tcells[i][j];
        }
        fwrite(&Etotal, sizeof(double), 1, fout_Tcell_traj);
        fwrite(&Tcells[i][rep_lim+1], sizeof(double), 1, fout_Tcell_traj);
        int index = 0;
        int count = 0;
        for(long position=epi_start[i]; position<epi_end[i]; position++){
            index += place_value[i][count]*population[0][position];
            count++;
        }
        chi[i][index] = chi[i][index]*N;
        fwrite(&chi[i][index], sizeof(double), 1, fout_Tcell_traj);
        chi[i][index] = chi[i][index]/N;
    }
                
    fwrite(&zero, sizeof(long), 1, fout_MC_seqs);
    for(long k=0; k<N; k++){
        for(long i=0; i<m; i++) {
            fwrite(&resIdx[i][population[k][i]], sizeof(int8_t), 1, fout_MC_seqs);
        }
    }

	
	// N.B. For execution acceleration, this subroutine operates under fake residue numbering scheme (fakeres) whereby 
	//      most probable residue at site i is coded as 0, next most probable as 1, ..., least probable as nRes[i]-1
	//      to obviate frequent (slow) lookup of residue types in resIdx
	//
	//		Under this scheme, the fakeres code is identical to its position in resIdx, facilitating rapid lookup of 
	//		corresponding n1,n2,h,J matrix elements
	//
		
    long pos = 0;
    long n_samples = 0;
    double Etemp[rep_lim];
    // \-> safeUnity guards against blue moon segmentation faults for ran2 returning precisely unity due to interaction with floor call (ran2 precisely zero is not problematic)
    static double safeUnity = 1.0 - std::numeric_limits<double>::epsilon();
	double E = 0;
	double penalty_S;
	int thread_num;
	#pragma omp parallel shared(thread_num)
	{
		if(omp_get_thread_num()==0)
			thread_num = omp_get_num_threads();
	}
	int seq_S;
#if PARTPARALLEL == 2 || PARTPARALLEL == 1 // unused in 1,0
	int *seq_arr = new int[thread_num];
#elif PARTPARALLEL == 3 || PARTPARALLEL == 0
	int *seq_arr = new int[N];
#endif
	vector<vector <int>> energy_idx; // This matrix of index is for parallelizing
										 // the energy calculation in PART 2
	init_energy_index(m, thread_num, energy_idx);
	

    // If omp is in use parallelize all operations on the population for speed.
    #pragma omp parallel default(none) shared(stdout, part_sum_used, seq_arr, energy_idx, E, seq_S, penalty_S, nRes, resIdx, h, J, T, population, temp_pop, temp_boltz, m, N, boltzfac, boltzfac_eff, Z, Z_eff, Ztemp, fout_P1_traj, fout_P2_traj, fout_MC_seqs, fout_pop_stats, fout_Tcell_traj, write_mod, print_mod, sample_mod, n_samples, cout, cerr, n_cycles, n1, n2, burnin, progeny, rate, seed, safeUnity, pos, Tcells, Etemp, n_T, chi, chiI, Ichi, T_penalty, n_WTepitopes, n_epitope, epi_start, epi_end, rep_lim, occupied, idx, part_sum, place_value)
    {
    
    // Parallel safe PRNG
    trng::yarn2 r((long)fabs(seed));        // PRNG
    trng::uniform01_dist<> u;               // uniform distribution for probabilities generates number between [0, 1]
    trng::poisson_dist IntGenerator(3*m*rate);  // poisson distribution for number of mutants in a sequence

    #ifdef _OPENMP
    const int size=omp_get_num_threads();    // get total number of processes
    const int rank=omp_get_thread_num();     // get rank of current process
	r.split(size, rank);               // split PRN sequences by leapfrog method - choose sub-stream no. rank out of size streams
	if(rank == 0)
		printf("OPENMP: %d threads\n", size);
    #else
    int size = 1;
    int rank = 0;
    #endif
	
    // break up the population into a chunk for each processor.
    const long len = N/size;
    const long begin = rank*len;
    const long end = rank + 1 == size ? N : begin+len;
	
	// break up DNA in later computation
	const long len_m = m/size;
    const long begin_m = rank*len_m;
    const long end_m = rank + 1 == size ? m : begin_m+len_m;
	
    // start main loop
    for (long cycle=1; cycle<=n_cycles; cycle++) 
	{
        
	// reset values used in every iteration of the loop
        #pragma omp single
        {
        pos = 0;
        Ztemp = 0;
        for (long i=0; i<n_epitope; i++)
		{
            for(long j=0; j<(n_WTepitopes[i]).size(); j++)
			{
                n_WTepitopes[i][j] = 0;
            }
       }
        }
       #pragma omp barrier

        int ranInt1 = 0;   // IntGenerator(r)=number of mution for the parent sequence;
        double ranNum1 = 0;     //=u(r)*safeUnity;   // range [0,1)
        long parent = 0;
        long site = 0;
        long res = 0;
        double Z_local = 0;
        double Ztemp_local = 0;
        double Z_local_eff = 0;
      
      
        // Have every sequence in the population produce offspring equal to the value of progeny. Mutating each offspring sequence as it is copied
        for(long k=0; k<progeny; k++)
      {       // loop over the number of progeny
            for(long i=begin; i<end; i++)
         {   // loop over the processor's part of the population
                ranInt1 = IntGenerator(r);   // number of mution for the parent sequence;
                for(long j=0; j<m; j++)
            {     // copy parent sequence
                    temp_pop[i+(N*k)][j] = population[i][j];
                }
                for(long j=0; j<ranInt1; j++)
            { // introduce mutations
                    ranNum1 = u(r)*safeUnity;
                    site = (long)floor((double)ranNum1 * (double)m);
                    ranNum1 = u(r)*safeUnity;
                    res = (long)floor((double)ranNum1 * (double)(nRes[site]-1));
                    if (temp_pop[i+(N*k)][site] <= res) 
               {
                  res++;
               }
                    temp_pop[i+(N*k)][site] = res;
                }
            
                idx[i+(N*k)] = i+(N*k);
                double penalty = 0;
                for(long j=0; j<n_epitope; j++)
            { // calculate T-cell susceptibility
                    int index = 0;
                    int count = 0;
                    for(long position=epi_start[j]; position<epi_end[j]; position++)
               {
                        index += place_value[j][count]*temp_pop[i+(N*k)][position];
                        count++;
                    }
                    double num_Ecells = 0;
                    for(long t=1; t<=rep_lim; t++)
               {
                        num_Ecells += Tcells[j][t];
                    }
                    penalty += T_penalty*chi[j][index]*num_Ecells/n_T;
                }
                // calculate boltzmann factor of daughter sequence and add it to the partition function
                temp_boltz[i+(N*k)] = exp(-(energy(temp_pop[i+(N*k)],m,h,J))/T-penalty/T);
                Ztemp_local += temp_boltz[i+(N*k)];
                occupied[i+(N*k)] = true;
            }
      }

        // combine partition function from each processor
        #pragma omp atomic
            Ztemp += Ztemp_local;
   
        #pragma omp barrier
      

#if PARTPARALLEL == 0      
        #pragma omp single
        {
            sort(idx, idx+N*progeny, [&temp_boltz](size_t i1, size_t i2){return temp_boltz[i1] < temp_boltz[i2];});
            part_sum[0] = temp_boltz[idx[0]];
            for(long i=1; i<N*progeny; i++){
                part_sum[i] = part_sum[i-1] + temp_boltz[idx[i]];
            }
      
            // pick sequences that survive to the final population based on their fitness
            Z = 0;
            Z_eff = 0;
            while(pos < N){
                ranNum1 = u(r)*safeUnity*Ztemp;
                long seq = std::distance(part_sum, std::lower_bound(part_sum, part_sum + N*progeny,ranNum1));
      seq_arr[pos] = seq;
                while(!occupied[idx[seq]]){
                    seq++;
                }
                if(seq >= N*progeny){
                    cout << "Error in find, seq = " << seq << ", ranNum1 = " << ranNum1 << ", Ztemp = " << Ztemp << endl;
                }
                for(long j=0; j<m; j++){
                    population[pos][j]=temp_pop[idx[seq]][j];
                }
                // prime T cells
                double penalty = 0;
                for(long j=0; j<n_epitope; j++){
                    int index = 0;
                    int count = 0;
                    for(long position=epi_start[j]; position<epi_end[j]; position++){
                        index += place_value[j][count]*population[pos][position];
                        count++;
                        }
                        double num_Ecells = 0;
                        for(long t=1; t<=rep_lim; t++){
                            num_Ecells += Tcells[j][t];
                        } 
                    penalty += T_penalty*chi[j][index]*num_Ecells/n_T;
                    n_WTepitopes[j][index]+=1.0;
                }
                boltzfac[pos] = exp(-(energy(population[pos],m,h,J))/T);
                boltzfac_eff[pos] = boltzfac[pos]*exp(-penalty/T);
                Z += boltzfac[pos];
                Z_eff += boltzfac_eff[pos];
                pos++;
                occupied[idx[seq]] = false;
                Ztemp -= temp_boltz[idx[seq]];
                for(long i=seq; i<N*progeny; i++){
                    part_sum[i] -= temp_boltz[idx[seq]];
                }
            }
      
        } // end #pragma omp single
        #pragma omp barrier
#endif

#if PARTPARALLEL == 1
///////////////////////////////
// PART 2 Introducing Mutation
///////////////////////////////

   #pragma omp single
   {
      // make sure smaller numbers don't get rounded off
      sort(idx, idx+N*progeny, [&temp_boltz](size_t i1, size_t i2){return temp_boltz[i1] < temp_boltz[i2];});
        part_sum[0] = temp_boltz[idx[0]];
        for(long i=1; i<N*progeny; i++)
      {
            part_sum[i] = part_sum[i-1] + temp_boltz[idx[i]];
        }
      // pick sequences that survive to the final population based on their fitness
        Z = 0;
        Z_eff = 0;
   }
   
        while(pos < N)
      {
          if(rank==0)
         //#pragma omp single
         {// only main thread start
         double ranNum1 = u(r)*safeUnity*Ztemp;
         long seq = std::distance(part_sum, std::lower_bound(part_sum, part_sum + N*progeny,ranNum1));
         seq_S = seq;
         while(!occupied[idx[seq]])
         {
            seq++;
         }
         if(seq >= N*progeny)
         {
            cout << "Error in find, seq = " << seq << ", ranNum1 = " << ranNum1 << ", Ztemp = " << Ztemp << endl;
         }
         for(long j=0; j<m; j++)
         {
            population[pos][j]=temp_pop[idx[seq]][j];
         }
         // prime T cells
         double penalty = 0;
         for(long j=0; j<n_epitope; j++)
         {
            int index = 0;
            int count = 0;
            for(long position=epi_start[j]; position<epi_end[j]; position++)
            {
               index += place_value[j][count]*population[pos][position];
               count++;
            }
            double num_Ecells = 0;
            for(long t=1; t<=rep_lim; t++)
            {
               num_Ecells += Tcells[j][t];
            } 
            penalty += T_penalty*chi[j][index]*num_Ecells/n_T;
            n_WTepitopes[j][index]+=1.0;
         }
         penalty_S = penalty;
         }// only main threadthread end 
         else
         {
            
         }
         
         #pragma omp barrier 
                  
         double E_P = 0;
         auto sequence = population[pos];
         for   (auto i=energy_idx[rank].begin(); i!=energy_idx[rank].end(); i++) 
         {
            E_P += h[*i][sequence[*i]];
            for (long j=(*i)+1; j<m; j++)
            {
               E_P += J[*i][j][sequence[*i]][sequence[j]];
            }
         }
         #pragma omp atomic
            E += E_P;
            
         #pragma omp barrier   
            
         if(rank==0)
         //#pragma omp single
         {   // only main thread start
            boltzfac[pos] = exp(-E/T);
            boltzfac_eff[pos] = boltzfac[pos]*exp(-penalty_S/T);
            Z += boltzfac[pos];
            Z_eff += boltzfac_eff[pos];
            pos++;
            occupied[idx[seq_S]] = false;
            Ztemp -= temp_boltz[idx[seq_S]];
         }// only main thread end

         if(size == 1 || rank != 0)
         {
             #pragma omp for schedule(static)
            for(long i=seq_S; i<N*progeny; i++)
            {
               part_sum[i] -= temp_boltz[idx[seq_S]];
            } 
         }
         #pragma omp barrier
        }// end while
   
//////////////////////////////
// END Part 2
/////////////////////////////


#endif 

#if PARTPARALLEL == 2
///////////////////////////////
// PART 2 Introducing Mutation
///////////////////////////////
   #pragma omp single
   {
      // make sure smaller numbers don't get rounded off
      sort(idx, idx+N*progeny, [&temp_boltz](size_t i1, size_t i2){return temp_boltz[i1] < temp_boltz[i2];});
        part_sum[0] = temp_boltz[idx[0]];
        for(long i=1; i<N*progeny; i++)
      {
            part_sum[i] = part_sum[i-1] + temp_boltz[idx[i]];
        }
      // pick sequences that survive to the final population based on their fitness
        Z = 0;
        Z_eff = 0;
   }

   //for(int k = 0; k < N; k+=size)
   while(pos < N)
   {
      int pos_L = pos + rank;
      double ranNum1 = u(r)*safeUnity*Ztemp;
      long seq = std::distance(part_sum, std::lower_bound(part_sum, part_sum + N*progeny,ranNum1));
      int active_threads = N > pos + (size - 1) ? size : N - pos;
      seq_arr[rank] = seq;
        #pragma omp barrier
      #pragma omp single
      {
         int dup_idx;
         while( (dup_idx = check_duplicate(seq_arr, active_threads)) != 0)
         {
            ranNum1 = u(r)*safeUnity*Ztemp;
            seq_arr[dup_idx] = std::distance(part_sum, std::lower_bound(part_sum, part_sum + N*progeny,ranNum1));
         }
      } 
      seq = seq_arr[rank];
      while(!occupied[idx[seq]])
      {
         seq++;
      }
      if(seq >= N*progeny)
      {
         cout << "Error in find, seq = " << seq << ", ranNum1 = " << ranNum1 << ", Ztemp = " << Ztemp << endl;
      }
      if(pos_L < N)
      {
         for(long j=0; j<m; j++)
         {
            population[pos_L][j]=temp_pop[idx[seq]][j];
         }
         // prime T cells
         double penalty = 0;
         for(long j=0; j<n_epitope; j++)
         {
            int index = 0;
            int count = 0;
            for(long position=epi_start[j]; position<epi_end[j]; position++)
            {
               index += place_value[j][count]*population[pos_L][position];
               count++;
            }
            double num_Ecells = 0;
            for(long t=1; t<=rep_lim; t++)
            {
               num_Ecells += Tcells[j][t];
            } 
            penalty += T_penalty*chi[j][index]*num_Ecells/n_T;
            #pragma omp atomic
               n_WTepitopes[j][index]+=1.0;
         }
         boltzfac[pos_L] = exp(-(energy(population[pos_L],m,h,J))/T); 
         boltzfac_eff[pos_L] = boltzfac[pos_L]*exp(-penalty/T); 
         
         #pragma omp critical
         {
            Z += boltzfac[pos_L]; 
            Z_eff += boltzfac_eff[pos_L]; 
            Ztemp -= temp_boltz[idx[seq]];
         }
         
         occupied[idx[seq]] = false;
      }
      
       for(int j = 0; j < active_threads; j++)
      {
         #pragma omp barrier
         #pragma omp for schedule(static)
         for(long i=seq_arr[j]; i<N*progeny; i++)
         {
            part_sum[i] -= temp_boltz[idx[seq_arr[j]]];
         }
      }

      #pragma omp single
      {
         pos += size;
      }
   }
   
//////////////////////////////
// END Part 2
/////////////////////////////

#endif

#if PARTPARALLEL == 3
///////////////////////////////
// PART 2 Introducing Mutation
///////////////////////////////


   // occupied here will denote the sequences that have already been removed from the array
   #pragma omp for schedule(static)
   for(int i = 0; i < N * progeny; i++)
   {
      part_sum_used[i] = false;
      occupied[i] = false;
   }
      
   
   #pragma omp single
   {
      // make sure smaller numbers don't get rounded off
      sort(idx, idx+N*progeny, [&temp_boltz](size_t i1, size_t i2){return temp_boltz[i1] < temp_boltz[i2];});
        part_sum[0] = temp_boltz[idx[0]];
        for(long i=1; i<N*progeny; i++)
      {
            part_sum[i] = part_sum[i-1] + temp_boltz[idx[i]];
        }
      // pick sequences that survive to the final population based on their fitness
        Z = 0;
        Z_eff = 0;
      
      // picking all sequences in parallel
       for(int i = 0; i < N; i++)
      {
         double ranNum1 = u(r)*safeUnity*Ztemp;
         int seq = std::distance(part_sum, std::lower_bound(part_sum, part_sum + N*progeny,ranNum1));
         int part_count = 0;
         
         // search for a different seq if it's used
         while(part_sum_used[seq] == true)
         {
            // escape condition: remove all used seq from array
            if(part_count > 1000)
            {
               for(int j = 0; j < N*progeny; j++)
               {
                  if(part_sum_used[j] == true && occupied[idx[j]] == false)
                  {
                     for(int k = j; k < N*progeny; k++)
                     {
                        part_sum[k] -= temp_boltz[idx[j]];
                     }
                     Ztemp -= temp_boltz[idx[j]];
                     occupied[idx[j]] = true;
                  }
               }
               part_count = 0;
            }
            //printf("seq: %d\n", seq);
            ranNum1 = u(r)*safeUnity*Ztemp;
            seq = std::distance(part_sum, std::lower_bound(part_sum, part_sum + N*progeny,ranNum1));            
            part_count++;
         }
         seq_arr[i] = seq;
         part_sum_used[seq] = true;
         if(seq >= N*progeny)
         {
            printf("error: seq = %d\n", seq);
         }
      }
      
/*       int im = 0;
      for(int in = N * progeny - N; in < N * progeny; in++)
      {
         int seq = rand() % in;
         if(part_sum_used[seq] == true)
            seq = in;
         part_sum_used[seq] = true;
         seq_arr[im] = seq;
         im++;
      }  */
   }
   
   
   
     #pragma omp for schedule(static) reduction(+: Z, Z_eff, Ztemp)
   for(int pos_L = 0; pos_L < N; pos_L++)
   {
      int seq = seq_arr[pos_L];
      for(long j=0; j<m; j++)
      {
         population[pos_L][j]=temp_pop[idx[seq]][j];
      }
      // prime T cells
      double penalty = 0;
      for(long j=0; j<n_epitope; j++)
      {
         int index = 0;
         int count = 0;
         for(long position=epi_start[j]; position<epi_end[j]; position++)
         {
            index += place_value[j][count]*population[pos_L][position];
            count++;
         }
         double num_Ecells = 0;
         for(long t=1; t<=rep_lim; t++)
         {
            num_Ecells += Tcells[j][t];
         } 
         penalty += T_penalty*chi[j][index]*num_Ecells/n_T;
         #pragma omp atomic
            n_WTepitopes[j][index]+=1.0;
      }
      boltzfac[pos_L] = exp(-(energy(population[pos_L],m,h,J))/T); 
      boltzfac_eff[pos_L] = boltzfac[pos_L]*exp(-penalty/T); 
      
      Z += boltzfac[pos_L]; 
      Z_eff += boltzfac_eff[pos_L];
      // only subtract from Ztemp if it hasn't been subtracted yet
      if(occupied[idx[seq]] == false)
         Ztemp += -temp_boltz[idx[seq]]; 
      
      occupied[idx[seq]] = true;
   }  
   
   #pragma omp barrier


//////////////////////////////
// END Part 2
/////////////////////////////
#endif

    #pragma omp for
    for(int i=0; i<n_epitope; i++){
        chiI[i] = 0.0;
        for(int j=0; j<(chi[i]).size(); j++){
            chiI[i] += chi[i][j]*n_WTepitopes[i][j];  
        }      
     }

#pragma omp single
{
     double totalT = 0;
     for(long i=0; i<n_epitope; i++)
    {
        Ichi = chiI[i];
      if(rep_lim < 2)
      { // effector generations equal 1
            integrate(int_Tcell_1, Tcells[i], 0.0, 1.0, 0.05);
        } 
      else if(rep_lim < 3)
      { // effector generations equal 2
            integrate(int_Tcell_2, Tcells[i], 0.0, 1.0, 0.05);
        } 
      else 
      {  // effector generations 3 or more
            integrate(int_Tcell_N, Tcells[i], 0.0, 1.0, 0.05);
        }
   // add up effector cells for T-cell cap
        for(long j=1; j<=rep_lim; j++){
                totalT +=  Tcells[i][j];
        }
     }

    // apply T-cell cap
    if(totalT > n_T){
        for(long i=0; i<n_epitope; i++){
            for(long j=1; j<rep_lim+1; j++){
                Tcells[i][j] = Tcells[i][j]/totalT*n_T;
            }
        }
    }
}
       
    #pragma omp barrier
    #pragma omp single
    {
        Z /= N;
        Z_eff /= N;
    }
        
    // sample block
    if (!(cycle % sample_mod) && cycle > burnin) 
   {
    for(long k=0; k<N; k++)
   {
        for(long i=begin_m; i<end_m; i++) 
      {
            n1[i][population[k][i]] += 1.0;
            for(long j=i+1; j<m; j++) 
         {
                n2[i][j][population[k][i]][population[k][j]] += 1.0;
           }
       }
   }
    #pragma omp single
    {
            n_samples += N;
   }
    }
        
    #pragma omp barrier
    #pragma omp single
    {
   // print block
   if (cycle % print_mod == 0 && cycle > 0) {
       cout << "  Completed cycle " << std::setw(5) << cycle << " of " << std::setw(5) << n_cycles << "\n";
   }
      
      
   // write block
        if (cycle % write_mod == 0 && cycle > burnin) 
      {
            fwrite(&cycle, sizeof(long), 1, fout_pop_stats);
            fwrite(&pos, sizeof(long), 1, fout_pop_stats);
            fwrite(&Z, sizeof(double), 1, fout_pop_stats);
            fwrite(&Ztemp, sizeof(double), 1, fout_pop_stats);
            fwrite(&Z_eff, sizeof(double), 1, fout_pop_stats);
            fwrite(&cycle, sizeof(long), 1, fout_P1_traj);
            for(long i=0; i<m; i++) {
                for(long p=0; p<nRes[i]; p++) {
                     n1[i][p] = n1[i][p]/n_samples;
                     fwrite(&n1[i][p], sizeof(double), 1, fout_P1_traj);
                     n1[i][p] = n1[i][p]*n_samples;
           }
       }
         
            fwrite(&cycle, sizeof(long), 1, fout_P2_traj);
            for(long i=0; i<m; i++) {
                for(long j=i+1; j<m; j++) {
                    for(long q=0; q<nRes[j]; q++) {
                        for(long p=0; p<nRes[i]; p++) {
                            n2[i][j][p][q] = n2[i][j][p][q]/n_samples;
                            fwrite(&n2[i][j][p][q], sizeof(double), 1, fout_P2_traj);
                            n2[i][j][p][q] = n2[i][j][p][q]*n_samples;
         }
          }
      }
       }

            fwrite(&cycle, sizeof(long), 1, fout_Tcell_traj);
            for(long i=0; i<n_epitope; i++){
                double Etotal = 0;
                fwrite(&Tcells[i][0], sizeof(double), 1, fout_Tcell_traj);
                for(long j=1; j<=rep_lim; j++){
                    Etotal += Tcells[i][j];
                }
                fwrite(&Etotal, sizeof(double), 1, fout_Tcell_traj);
                fwrite(&Tcells[i][rep_lim+1], sizeof(double), 1, fout_Tcell_traj);
                fwrite(&chiI[i], sizeof(double), 1, fout_Tcell_traj);
            }
         
            fwrite(&cycle, sizeof(long), 1, fout_MC_seqs);
       for(long k=0; k<N; k++){
                for(long i=0; i<m; i++) {
                    fwrite(&resIdx[i][population[k][i]], sizeof(int8_t), 1, fout_MC_seqs);
           }
       }
   }
        } // end #pragma omp single
        #pragma omp barrier
    } // end main loop
    } // end #pragma omp parallel
    cout << "DONE!" << endl << endl;

   // free memory
   delete []idx;
   delete []part_sum;
   delete []boltzfac_eff;
   delete []boltzfac;
   delete []pop_init;
   delete []occupied;
//   delete []seq_arr;
   delete []part_sum_used;
   
    // closing traj files
    fclose(fout_P1_traj);
    fclose(fout_P2_traj);
    fclose(fout_MC_seqs);
    fclose(fout_pop_stats);
    fclose(fout_Tcell_traj);
   
   
    // reporting final probabilities
    FILE* fout_P1;
    {
   std::string fstr_P1="./P1_model.dat";
        fout_P1 = fopen(fstr_P1.c_str(),"wb");
        if (fout_P1==NULL) {
       cerr << "Cannot open output file " << fstr_P1 << " in the current directory; aborting." << endl;
       exit(-1);
   }
        int32_t sizes [1];
        sizes[0] = sizeof(double);
        fwrite(sizes,4,1,fout_P1);    
    } 
   

    FILE* fout_P2;
    {
   std::string fstr_P2="./P2_model.dat";
        fout_P2 = fopen(fstr_P2.c_str(),"wb");
        if (fout_P2==NULL) {
       cerr << "Cannot open output file " << fstr_P2 << " in the current directory; aborting." << endl;
       exit(-1);
   }
        int32_t sizes [1];
        sizes[0] = sizeof(double);
        fwrite(sizes,4,1,fout_P2);    
    } 

    for(long k=0; k<N; k++){
        for(long i=0; i<m; i++) {
            n1[i][population[k][i]] += 1.0;
            for(long j=i+1; j<m; j++) {
                n2[i][j][population[k][i]][population[k][j]] += 1.0;
       }
   }
    }
    n_samples += N;
   
    for(long i=0; i<m; i++) {
        for(long p=0; p<nRes[i]; p++) {
              n1[i][p] = n1[i][p]/n_samples;
              fwrite(&n1[i][p], sizeof(double), 1, fout_P1);
   }
    }
      
    for(long i=0; i<m; i++) {
        for(long j=i+1; j<m; j++) {
            for(long q=0; q<nRes[j]; q++) {
                for(long p=0; p<nRes[i]; p++) {
                    n2[i][j][p][q] = n2[i][j][p][q]/n_samples;
                    fwrite(&n2[i][j][p][q], sizeof(double), 1, fout_P2);
      }
       }
   }
    }
   
    fclose(fout_P1);
    fclose(fout_P2);

#if PARTPARALLEL == 0 || PARTPARALLEL == 3
    FILE* seq_file;
    // write unsorted
    seq_file = fopen("seq_unsorted.txt", "w");
    for(long i =0; i < N; i++)
       fprintf(seq_file, "%d\n", seq_arr[i]);
    fclose(seq_file);

    // write sorted
    sort(seq_arr, seq_arr + N);
    seq_file = fopen("seq_sorted.txt", "w");
    for(long i =0; i < N; i++)
       fprintf(seq_file, "%d\n", seq_arr[i]);
    fclose(seq_file);
#endif // check sequence
    delete []seq_arr;


   
    // stop timer
    double wall_stop = get_wall_time();
    double cpu_stop  = get_cpu_time();
   
    printf("Wall Time = %e s\n",wall_stop - wall_start);
    printf("CPU Time  = %e s\n",cpu_stop - cpu_start);
   
    return 0;
}
