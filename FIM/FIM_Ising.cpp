#include <armadillo>
#include <iostream>
#include <cstdio>

#define TOL 1E-10
#define LAMBDA_H 0
#define LAMBDA_J 0

double my_mean(arma::Col<int32_t> col, int32_t size) {
  double mean = 0.0;
  for (int32_t i = 0; i < size; i++) {
    mean += col[i];
  }
  return mean / size;
}

int main(int argc, char** argv) {

  int32_t N;
  int32_t M;
  int32_t K;
  int32_t NEIGS;

  arma::Mat<int32_t>* msa;
  if (argc-1 == 1) {
    FILE* input = fopen(argv[1], "r");

    fscanf(input, "%d", &N);
    fscanf(input, "%d", &M);
    fscanf(input, "%d", &NEIGS);
    NEIGS = NEIGS > M + M * (M - 1) / 2 ? M + M * (M - 1) / 2 : NEIGS;

    // initialize msa array
    msa = new arma::Mat<int32_t>(N, M);
    for (int32_t i = 0; i < N; i++) {
      for (int32_t j = 0; j < M; j++) {
        int32_t num;
        fscanf(input, "%d", &num);
        (*msa)(i, j) = num;
      }
    }
    fclose(input);
  } else {
    sscanf(argv[1], "%d", &N);
    sscanf(argv[2], "%d", &M);
    sscanf(argv[3], "%d", &NEIGS);
    NEIGS = NEIGS > M + M * (M - 1) / 2 ? M + M * (M - 1) / 2 : NEIGS;
    
    arma::arma_rng::set_seed_random();
    msa = new arma::Mat<int32_t>(arma::randi<arma::Mat<int32_t>>(N, M, arma::distr_param(0, 1))); 
  }

  K = N;

  int32_t fim_dim = M + M * (M - 1) / 2;
  arma::SpMat<double>* fim = new arma::SpMat<double>(fim_dim, fim_dim);
  int32_t omp_num_threads = omp_get_num_threads();

  std::cout << "threads: " << omp_num_threads << std::endl;

  #pragma omp parallel for schedule(dynamic) num_threads(M > omp_num_threads ? omp_num_threads : M) 
  for (int32_t i = 0; i < M; i++) {
    for (int32_t j = i; j < M; j++) {
      for (int32_t k = 0; k < M; k++) {
        for (int32_t l = k; l < M; l++) {
          int32_t row = i * M - i * (i - 1) / 2 + (j - i);
          int32_t col = k * M - k * (k - 1) / 2 + (l - k);
          if (col < row) {
            continue;
          }
          arma::Col<int32_t> col_i = (*msa).col(i);
          arma::Col<int32_t> col_j = (*msa).col(j);
          arma::Col<int32_t> col_k = (*msa).col(k);
          arma::Col<int32_t> col_l = (*msa).col(l);

          double p_ijkl = my_mean(col_i % col_j % col_k % col_l, N);
          double p_ij = my_mean(col_i % col_j, N);
          double p_kl = my_mean(col_k % col_l, N);

          double fim_row_col = 0.0;
          if (i == k && j == l) {
            if (i == j) {
              fim_row_col = -(p_ij * p_kl - p_ijkl - 2 * LAMBDA_H / K);
            } else {
              fim_row_col = -(p_ij * p_kl - p_ijkl - 2 * LAMBDA_J / K);
            }
          } else {
            fim_row_col = -(p_ij * p_kl - p_ijkl);
          }

          if (fim_row_col > TOL) {
            (*fim)(row, col) = fim_row_col;
            (*fim)(col, row) = fim_row_col;
          }
        }
      }
    }
  }

  delete msa; 

  std::cout << "Filled out FIM" << std::endl;
  
  std::cout << "dimension: " <<  fim_dim << " by " << fim_dim << std::endl;
  std::cout << "density: " << (float) nonzeros(*fim).size() * 100.0 / (*fim).size() << "\%" << std::endl;
 
  arma::Col<double> eigval;
  arma::Mat<double> eigvec;
  arma::eigs_sym(eigval, eigvec, *fim, NEIGS);
  eigval.print("eigval:");
  eigvec.print("eigvec:");
  delete fim;
}
