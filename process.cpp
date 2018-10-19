/*
 *  process.cpp
 *  process Population_Dynamics data
 *
 *  Created by Greg Hart on Sep 2016
 *  Copyright 2015 ALF. All rights reserved.
 *
 */
/*
 * IN
 * resIdx.dat         - ordered list of residue types - as integer codes - present at each site from most to least probable
 * MC_seqs.dat         - sequences sampled from population trajectory
 * pop_stats.dat    - data on the population at different time points
 * epitopes.dat
 * inputs.dat
 */
 
#include "process.h"
//#include "functions.h"

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

//#include <boost/numeric/odeint.hpp>

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::vector;
using std::sort;
//using namespace boost::numeric::odeint;


// main function preforming Fisher-wright dynamics on the viral side and intergrating ODEs for the Tcell side
int main (int argc, char *argv[]) {
    
    cout.setf(std::ios_base::scientific);
    cout.precision(3);    
    
    #ifdef UNIT_TESTING
    
    #endif  // #ifdef UNIT_TESTING

    
      // start timer
//    double wall_start = get_wall_time();
    
    
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
    long rep_lim = 9;
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

   // Load in MC seps
   long n_steps;
   long N_read;
   int32_t sizes[2];
   FILE* fout_MC_seqs;
   {
        std::string fstr_MC_seqs="./MC_seqs.dat";
        std::ifstream resume;
        resume.open(fstr_MC_seqs.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
        long size = resume.tellg();
        resume.close();
        fout_MC_seqs = fopen(fstr_MC_seqs.c_str(),"rb");
        if (fout_MC_seqs==NULL) {
        cerr << "Cannot open output file " << fstr_MC_seqs << " in the current directory; aborting." << endl;
        exit(-1);
        }
        fread(sizes, 4, 2, fout_MC_seqs);
        fread(&N_read, sizes[0], 1, fout_MC_seqs);
        if(N != N_read){
            cerr << "***WARNING*** Read population size differs from that in the inputs file" << endl;
            N = N_read;
        }
        n_steps = (size-2)/(sizes[0] + sizes[1]*m*N);
    }    

    int8_t *population = (int8_t *)malloc(n_steps*N*m * sizeof(int8_t));
    long time[n_steps];
    for(long i=0; i<n_steps; i++){
        fread(&time[i], sizes[0], 1, fout_MC_seqs);
        fread(&population[i*N*m], sizes[1],N* m, fout_MC_seqs);
    }
    fclose(fout_MC_seqs);

    // find epitope section of MC seqs and write to file
    int epitope_start[n_epitope];
    int epitope_end[n_epitope];
    std::string fstr = "./epitopes.dat";
    std::ifstream fin;
    fin.open(fstr.c_str());
    if (!fin) {
        std::cerr << "Cannot open input file " << fstr << "; aborting." << std::endl;
        exit(-1);
    }
    
    std::string epi_file;
    for(long i=0; i<n_epitope; i++) {
        if (fin.eof()) {
           std::cerr << "Prematurely encountered eof in " << fstr << "; aborting." << std::endl;
           exit(-1);
        }
        
        std::getline(fin, epi_file);
        size_t pos = 0;
        pos = epi_file.find_last_of('/');
        epi_file = epi_file.substr(pos+1);
        std::string token;
        pos = epi_file.find('_');
        epi_file.erase(0, pos + 1);
        pos = epi_file.find('_');
        token = epi_file.substr(0, pos);
        pos = token.find('-');
        epitope_start[i] = atoi(token.substr(0,pos).c_str())-1;
        epitope_end[i] = atoi(token.substr(pos+1).c_str());
        FILE* fout_epi_traj;
        {
            fout_epi_traj = fopen(epi_file.c_str(),"wb");
            if (fout_epi_traj==NULL) {
                cerr << "Cannot open output file " << epi_file << " in the current directory; aborting." << endl;
                exit(-1);
            }
            int32_t sizes [2];
            sizes[0] = sizeof(long);
            sizes[1] = sizeof(int8_t);
            fwrite(sizes,4,2,fout_epi_traj);    
            fwrite(&N,sizeof(long),1,fout_epi_traj);
            for(long j=0; j<n_steps; j++){
                for(long k=0; k<N; k++){
                    fwrite(&time[j],sizes[0],1,fout_epi_traj);
                    fwrite(&population[j*N*m+k*m+epitope_start[i]],sizes[1],epitope_end[i]-epitope_start[i],fout_epi_traj);
                }
            }
        } 
    }
    
     while(!fin.eof()) {
          std::getline(fin,epi_file);
          if(!epi_file.empty()){
             std::cerr << "Elements remain in " << fstr << " after populating array; aborting." << std::endl;
             exit(-1);
        }
    }
    
    fin.close();

    //entropy calculation
    long idx[N];
    bool used[N];
    for(long i=0;i<N;i++){
        idx[i] = i;
        used[i] = false;
    }
    std::random_shuffle(&idx[0],&idx[N]);
    double entropy_traj[n_steps][10];
    for(long i=0; i<n_steps; i++){
        for(long j=0; j<10; j++){
            for(long k=0; k<floor(N*((j+1)/10)); k++){
                double P = 0;
                if(!used[k]){
                    for(long l=0; l<floor(N*((j+1)/10)); l++){
                        bool same = true;
                        if(!used[l]){
                            for(long p=0; p<m; p++){
                                if(population[i*N*m + idx[k]*m+p]!=population[i*N*m + idx[l]*m+p]){
                                    same = false;
                                    break;
                                }
                            }
                            if(same){
                                P++;
                                used[l] = true;
                            }
                        }
                    }
                    P /= floor(N*((j+1)/10));
                    entropy_traj[i][j] += P*log(P);
                }
            }
        }
    }

    FILE* fout_S_traj;
    {
        fout_S_traj = fopen("entropy_traj.dat","wb");
        if (fout_S_traj==NULL) {
           cerr << "Cannot open output file " << "entorpy_traj.dat" << " in the current directory; aborting." << endl;
           exit(-1);
        }
        int32_t sizes [2];
        sizes[0] = sizeof(long);
        sizes[1] = sizeof(double);
        fwrite(sizes,4,2,fout_S_traj);    
        for(long j=0; j<n_steps; j++){
            fwrite(&time[j],sizes[0],1,fout_S_traj);
            fwrite(entropy_traj[j], sizes[1], 10, fout_S_traj);
        }
    }
 
    // stop timer
//    double wall_stop = get_wall_time();
    
//    printf("Wall Time = %e s\n",wall_stop - wall_start);

    free(population);
    cout << "DONE" << endl;
    
    return 0;
}
