/*
 *  functions.h
 *  viralMC
 *
 *  Created by Andrew Ferguson on 4/15/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <string>
#include <vector>
#include <stdint.h>

double get_wall_time();
double get_cpu_time();
int load_seq(std::string fstr, int seq[], std::vector<long>& nRes, std::vector< std::vector<long> >& resIdx, long m);
int load_epitopes(std::string fstr, std::vector<long>& epitope_start, std::vector<long>& epitope_end, std::vector< std::vector<double> >& chi, long n, std::vector< std::vector<long> >& resIdx, std::vector<long>& nRes);
int load_epitopes(std::string fstr, std::vector<long>& epitope_start, std::vector<long>& epitope_end, std::vector< std::vector<int> >& epitopes, long n, std::vector<double>& chi, std::vector< std::vector<long> >& resIdx); 
int load_resIdx(std::string fstr, std::vector<long>& nRes, std::vector< std::vector<long> >& resIdx, long m);
int load_X1(std::string fstr, std::vector<long>& nRes, std::vector< std::vector<double> >& X1, long m);
int load_X2(std::string fstr, std::vector<long>& nRes, std::vector< std::vector< std::vector< std::vector<double> > > >& X2, long m);
double energy(std::vector<int8_t>& sequence, long m, std::vector< std::vector<double> >& h, std::vector< std::vector< std::vector< std::vector<double> > > >& J);
double Tresponse(double affinity, double num_Ecells, double b, long N, double aff_rate);

#endif
