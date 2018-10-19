/*
 *  functions.cpp
 *  viralMC
 *
 *  Created by Andrew Ferguson on 4/15/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include "functions.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <sys/time.h>


double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

int load_seq(std::string fstr, int seq[], std::vector<long>& nRes, std::vector< std::vector<long> >& resIdx, long m) {

	std::ifstream fin;
	fin.open(fstr.c_str());
	if (!fin) {
		std::cout << "Cannot open input file " << fstr << "; Using WT." << std::endl;
		for(long i=0; i<m; i++){
			seq[i] = 0;
		}
		return 0;
	}
	
	std::string readStr;
	std::getline(fin,readStr);
	std::stringstream readStrStrm(readStr);
	
	std::string word;
	for(long i=0; i<m; i++) {
		if(readStrStrm.eof()){
			std::cerr << "Encountered end of line read from " << fstr << "before fully populationg arry; aborting." << std::endl;
			exit(-1);
			}
		readStrStrm >> word;
		seq[i] = atof(word.c_str());
		for(long j=0; j<nRes[i]; j++){
			if(seq[i] == resIdx[i][j]){
				seq[i] = j;
				break;
			}
		}
	}

	if (!readStrStrm.eof()) {
		std::cerr << "Elements remain in " << fstr << " after populating array; aborting." << std::endl;
		exit(-1);
	}
	
	fin.close();
		
	return 0;
}

int load_epitopes(std::string fstr, std::vector<long>& epitope_start, std::vector<long>& epitope_end, std::vector< std::vector<double> >& chi, long n, std::vector< std::vector<long> >& resIdx, std::vector<long>& nRes, std::vector<state_type>& Tcells, int rep_lim) {

     std::ifstream fin;
	fin.open(fstr.c_str());
	if (!fin) {
		std::cerr << "Cannot open input file " << fstr << "; aborting." << std::endl;
		exit(-1);
	}
	
	std::string epi_file;
	for(long i=0; i<n; i++) {
	    if (fin.eof()) {
	        std::cerr << "Prematurely encountered eof in " << fstr << "; aborting." << std::endl;
		exit(-1);
	    }
		
	    std::ifstream epi_in;
	    {
            std::getline(fin, epi_file);
	    std::stringstream readStrStrm(epi_file);
	    std::string word;
            if(readStrStrm >> word){
	        epi_in.open(word.c_str());
                if (!epi_in) {
	           std::cerr << "Cannot open input file " << epi_file << "; aborting." << std::endl;
	           exit(-1);
	        }
            } else {
	        std::cerr << "Not enough lines in " << fstr << "; aborting." << std::endl;
	        exit(-1);
            }
            Tcells[i][0] = 200.0;//500;
            Tcells[i][rep_lim+1] = 0.0;
            for(long j=1; j<=rep_lim; j++){
                Tcells[i][j] = 0.0;
            }
	    if(readStrStrm >> word){
	        Tcells[i][0] = atof(word.c_str());
	    } else {
	        std::cerr << "Number of cells for epitope " << i << " not given. Using 200.0 Naive, 0.0 Effector, and 0.0 Memory." << std::endl;
	    }
	    if(readStrStrm >> word){
	        Tcells[i][1] = atof(word.c_str());
	    }
	    if(readStrStrm >> word){
	        Tcells[i][rep_lim+1] = atof(word.c_str());
	    } 
	    if(readStrStrm >> word) {
	        std::cerr << "Extra entries for epitope " << i << "; aborting." << std::endl;
	        exit(-1);
	    }
            }

	    std::string readStr;
	    std::getline(epi_in,readStr);
            std::stringstream readStrStrm(readStr);
	    std::string word;
	    double average = 0;
	    if(readStrStrm >> word){
	        average = atof(word.c_str());
	    } else {
	        std::cerr << "Empty line in " << epi_file << "; aborting." << std::endl;
	        exit(-1);
            }
	    word = "";
	    readStrStrm >> word;
	    if(!word.empty()){
	        std::cerr << "Extra characters on the first line of " << epi_file << "; aborting." << std::endl;
	        exit(-1);
	    }
		
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
	     int place_value[epitope_end[i]-epitope_start[i]];
	     int size = 1;
	     int count = 0;
	     for(int j=epitope_start[i]; j<epitope_end[i]; j++){
	          place_value[count] = size;
	          size *= nRes[j];
	          count ++;
	     }
	     (chi[i]).resize(size, average);
	     
	     int j=0;
	     while(std::getline(epi_in,readStr)){
	          int index = 0;
	          readStrStrm.clear();
     		readStrStrm.str(readStr);
     		for(int k=0; k<epitope_end[i]-epitope_start[i]; k++){
     		     if(!(readStrStrm >> word)){
     		          std::cerr << "Too few residues on line " << j+2 << " of " << epi_file << "; aborting." << std::endl;
     		          exit(-1);
     		     }
     		     for(int p=0; p<resIdx[epitope_start[i]+k].size(); p++){
	                    if(atoi(word.c_str()) == resIdx[epitope_start[i]+k][p]){
	                         index += place_value[k]*p;
	                         break;
	                    }
	               }
     		}
     		if(readStrStrm >> word){
     		     chi[i][index] = atof(word.c_str());
     		} else {
     		     std::cerr << "No chi on line " << j+2 << " of " << epi_file << "; aborting." << std::endl;
     		     exit(-1);
     		}
     		word = "";
		     readStrStrm >> word;
     		if(!word.empty()){
		          std::cerr << "Too many residues on line " << j+2 << " of " << epi_file << "; aborting." << std::endl;
		          exit(-1);
		     }
		     j++;
	     }
	     std::cout << i+1 << "th epitope starts at " << epitope_start[i] << " and ends at " << epitope_end[i] << std::endl;
	     while(!epi_in.eof()) {
               std::getline(epi_in,readStr);
               if(!readStr.empty()){
          		std::cerr << "Elements remain in " << epi_file << " after populating array; aborting." << std::endl;
	          	exit(-1);
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
	
     return 0;
}

int load_epitopes(std::string fstr, std::vector<long>& epitope_start, std::vector<long>& epitope_end, std::vector< std::vector<int> >& epitopes, long n, std::vector<double>& chi, std::vector< std::vector<long> >& resIdx) {

     std::ifstream fin;
	fin.open(fstr.c_str());
	if (!fin) {
		std::cerr << "Cannot open input file " << fstr << "; aborting." << std::endl;
		exit(-1);
	}
	
	std::string readStr;
	for(long i=0; i<n; i++) {
		if (fin.eof()) {
			std::cerr << "Prematurely encountered eof in " << fstr << "; aborting." << std::endl;
			exit(-1);
		}
		
		std::getline(fin,readStr);
		std::stringstream readStrStrm(readStr);
		
		std::string word;
		if(readStrStrm >> word){
		    chi[i] = atof(word.c_str());
		} else {
		     std::cerr << "Empty line in " << fstr << "; aborting." << std::endl;
		}
		if(readStrStrm >> word){
		     epitope_start[i] = atol(word.c_str());
		} else {
		     std::cerr << "Line ended early in " << fstr << "; aborting." << std::endl;
		}
		if(readStrStrm >> word){
		     epitope_end[i] = atol(word.c_str());
		} else {
		     std::cerr << "Missing the second position on line " << i+1 << " of " << fstr << "; aborting." << std::endl;
		}
		int counter = 0;
		while (readStrStrm >> word) {
			(epitopes[i]).push_back(atoi(word.c_str()));
			counter++;
		}
		if(counter > epitope_end[i]-epitope_start[i]+1){
		     std::cerr << "Too many residues on line " << i+1 << " of " << fstr << "; aborting." << std::endl;
		} else if(counter <epitope_end[i]-epitope_start[i]+1){
		     std::cerr << "Too few residues on line " << i+1 << " of " << fstr << "; aborting." << std::endl;
		}
	}

	if (!fin.eof()) {
		std::cerr << "Elements remain in " << fstr << " after populating array; aborting." << std::endl;
		exit(-1);
	}
	
	fin.close();
	
	for(long i=0; i<n; i++){
	     for(long j=0; j<epitope_end[i]-epitope_start[i]+1; j++){
	          for(long k=0; k<resIdx[epitope_start[i]+j-1].size(); k++){
	               if(epitopes[i][j] == resIdx[epitope_start[i]+j-1][k]){
	                    epitopes[i][j] = k;
	                    break;
	               }
	          }
	     }
	}
	
     return 0;
}

int load_resIdx(std::string fstr, std::vector<long>& nRes, std::vector< std::vector<long> >& resIdx, long m) {

	std::ifstream fin;
	fin.open(fstr.c_str());
	if (!fin) {
		std::cerr << "Cannot open input file " << fstr << "; aborting." << std::endl;
		exit(-1);
	}
	
	std::string readStr;
	for(long i=0; i!=m; i++) {
		if (fin.eof()) {
			std::cerr << "Prematurely encountered eof in " << fstr << "; aborting." << std::endl;
			exit(-1);
		}
		
		std::getline(fin,readStr);
		std::stringstream readStrStrm(readStr);
		
		std::string word;
		readStrStrm >> word;	// dumping position index
		long cntr=0;
		while (readStrStrm >> word) {
			(resIdx[i]).push_back(atof(word.c_str()));
			cntr++;
		}
		nRes[i]=cntr;
	}

	if (!fin.eof()) {
		std::cerr << "Elements remain in " << fstr << " after populating array; aborting." << std::endl;
		exit(-1);
	}
	
	fin.close();
	
	return 0;
}

int load_X1(std::string fstr, std::vector<long>& nRes, std::vector< std::vector<double> >& X1, long m) {

	std::ifstream fin;
	fin.open(fstr.c_str());
	if (!fin) {
		std::cerr << "Cannot open input file " << fstr << "; aborting." << std::endl;
		exit(-1);
	}
	
	std::string readStr;
	std::getline(fin,readStr);
	std::stringstream readStrStrm(readStr);
	
	std::string word;
	for(long i=0; i!=m; i++) {
		for(long p=0; p!=nRes[i]; p++) {
			if (readStrStrm.eof()) {
				std::cerr << "Encountered end of line read from " << fstr << " before fully populating array; aborting." << std::endl;
				exit(-1);
			}
			readStrStrm >> word;
			(X1[i][p]) = atof(word.c_str());
		}
	}
	
	if (!readStrStrm.eof()) {
		std::cerr << "Elements remain in line read from " << fstr << " after populating array; aborting." << std::endl;
		exit(-1);
	}

	if (!fin.eof()) {
		std::cerr << "More than one line detected in " << fstr << "; aborting." << std::endl;
		exit(-1);
	}
	
	fin.close();
	
	return 0;
}

int load_X2(std::string fstr, std::vector<long>& nRes, std::vector< std::vector< std::vector< std::vector<double> > > >& X2, long m) {
	
	std::ifstream fin;
	fin.open(fstr.c_str());
	if (!fin) {
		std::cerr << "Cannot open input file " << fstr << "; aborting." << std::endl;
		exit(-1);
	}
	
	std::string readStr;
	std::getline(fin,readStr);
	std::stringstream readStrStrm(readStr);
	
	std::string word;
	for(long i=0; i!=m; i++) {
		for(long j=i+1; j!=m; j++) {
			for(long q=0; q!=nRes[j]; q++) {
				for(long p=0; p!=nRes[i]; p++) {
					if (readStrStrm.eof()) {
						std::cerr << "Encountered end of line read from " << fstr << " before fully populating array; aborting." << std::endl;
						exit(-1);
					}
					readStrStrm >> word;
					X2[i][j][p][q] = atof(word.c_str());	// column major order file load
				}
			}
		}
	}
	
	if (!readStrStrm.eof()) {
		std::cerr << "Elements remain in line read from " << fstr << " after populating array; aborting." << std::endl;
		exit(-1);
	}

	if (!fin.eof()) {
		std::cerr << "More than one line detected in " << fstr << "; aborting." << std::endl;
		exit(-1);
	}
	
	fin.close();
	
	return 0;
}

double energy(std::vector<int8_t>& sequence, long m, std::vector< std::vector<double> >& h, std::vector< std::vector< std::vector< std::vector<double> > > >& J){
    	double E=0;
    	for	(long i=0; i<m; i++) {
          E = E + h[i][sequence[i]];
          for (long j=i+1; j<m; j++){
	     	E = E + J[i][j][sequence[i]][sequence[j]];
          }
		}
               
     return E;
}

double energy_P(std::vector<int8_t>& sequence, long m, std::vector< std::vector<double> >& h, std::vector< std::vector< std::vector< std::vector<double> > > >& J){
    	double E=0;
		
		#pragma omp parallel shared(E)
		{
			int size = omp_get_num_threads();
			int rank = omp_get_thread_num();
			
			int len = (m-1)/size + 1;
			int begin = rank * len;
			int end = begin + len;
			if(end > m)
				end = m;
			double E_P = 0;

				for	(long i=begin; i<end; i++) 
				{
					E_P = E_P + h[i][sequence[i]];
					for (long j=i+1; j<m; j++)
					{
						E = E + J[i][j][sequence[i]][sequence[j]];
					}
				}
			#pragma omp atomic
				E += E_P;
			#pragma omp barrier
     }
               
     return E;
}

double Tresponse(double affinity, double num_Ecells, double b, long N, double aff_rate){
     double k = .01;//3;
     double mid = .1*N;
     return b*affinity*1.0/(1.0 + exp(-k*(num_Ecells-mid)));
}

void init_energy_index(int m, int size, std::vector<std::vector<int>> & arr_idx){
    arr_idx.resize(size);
    int i = 0;
    bool rev = true;
    while(i < m) {
        for(int j = 0; j < size; j++) {
	    if(i < m)
	        arr_idx[j].push_back(i);
	    else
                break;
	    i++;
        }
    }
}

int check_duplicate(int* arr, int size){
	std::sort(arr, arr+size);
	int prev = arr[0];
	for(int i = 1; i < size; i++)
	{
		if(arr[i] == prev) return i;
		else prev = arr[i];
	}
	return 0;
}
