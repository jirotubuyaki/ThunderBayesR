// [[Rcpp::depends(RcppGSL)]]
//----------------------------------------------------------------------------------------Welcome
//                     Chinese Restaurant Process MCMC
//                                                                   2020/04/22
//                                                                 Masashi Okada
//                                                      okadaalgorithm@gmail.com
//
//
//
// Description: Chinese Restaurant Process is used in order to compose Dirichlet
// Process. The Clustering which uses Chinese Restaurant Process does not need to
// decide a number of clusters in advance. This algorithm automatically adjusts
// it. Main function is gibbs()

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <sys/time.h>
#include <RcppGSL.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_statistics_int.h>

using namespace Rcpp;

// [[Rcpp::export]]
List gibbs(const int data_length, const int dim, const RcppGSL::matrix<double> &data_raw, const RcppGSL::matrix<double> &data_raw_result, const double alpha, const int burn_in, const int iteration) {
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  struct timeval tv;
  gettimeofday(&tv,0);
  unsigned long mySeed = tv.tv_sec + tv.tv_usec;
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r, mySeed);

  RcppGSL::matrix<double> data(data_length, dim);
  RcppGSL::matrix<double> data_result(data_length, dim);

  for(int i = 0; i < data_length; i++){
    for(int j = 0; j < dim; j++){
      data(i, j) = data_raw(i, j);
      data_result(i, j) = data_raw_result(i, j);
    }
  }

  RcppGSL::vector<double> mu_0(dim);
  RcppGSL::matrix<double> sigma_new_init(dim, dim);
  for(int i = 0; i < dim; i++){
    double col[data_length];
    RcppGSL::vector_view<double> col_view = gsl_matrix_column(data, i);
    for(int j = 0; j < data_length; j++){
      col[j] = col_view[j];
    }
    double mean = gsl_stats_mean(col, 1, data_length);
    mu_0[i] = mean;
  }
  for (int i = 0; i < dim; i++) {
    for (int j = i; j < dim; j++) {
      RcppGSL::vector_view<double> col_view_i = gsl_matrix_column(data, i);
      RcppGSL::vector_view<double> col_view_j = gsl_matrix_column(data, j);
      double col_i[data_length];
      double col_j[data_length];
      for(int k = 0; k < data_length; k++){
        col_i[k] = col_view_i[k];
        col_j[k] = col_view_j[k];
      }
      double cov = gsl_stats_covariance(col_i, 1,col_j, 1, data_length);
      sigma_new_init(i, j) = cov;
      sigma_new_init(j, i) = cov;
    }
  }
  RcppGSL::vector<double> eval(dim);
  RcppGSL::Matrix evec(dim, dim);
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(dim);
  gsl_eigen_symmv(sigma_new_init, eval, evec, w);
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC); 
  int flag = 0;
  for (int k = 0; k < dim; k++) {
    if(eval[k] <= 0){
      flag = 1;
      Rcout << "Warning: Dataset's covariance-matrix is not positive definite." << "\n";
      for(int l = 0; l < dim; l++){
        for(int m = 0; m < dim; m++){
          sigma_new_init(l, m) = 1;
        }
      }
    }
  }
  eval.free();
  evec.free();
  RcppGSL::matrix<double> sigma_new(dim, dim);
  RcppGSL::matrix<double> sigma_new_mem(dim, dim);

  if(flag == 0){
    double s;
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < i; j++) {
      s = sigma_new_init(i,j);
        for (int k = 0; k < j; k++){
          s = s - sigma_new(i,k) * sigma_new(j,k);
        }
        sigma_new(i,j) = s / sigma_new(j,j);
      }
      s = sigma_new_init(i,i);
      for (int k = 0; k < i; k++){
        double sigma_new_in = sigma_new(i, k);
        s = s - pow(sigma_new_in, 2.0);
      }
      sigma_new(i,i) = sqrt(s);
    }
  }
  for(int l = 0; l < dim; l++){
    for(int m = 0; m < dim; m++){
      sigma_new_mem(l, m) = sigma_new(l, m);
    }
  }
  std::vector<std::vector<double>> mu_k;
  std::vector<std::vector<std::vector<double>>> sigma_k;
  mu_k.push_back(std::vector<double>(dim));
  mu_k.push_back(std::vector<double>(dim));
  sigma_k.push_back(std::vector<std::vector<double>>(dim, std::vector<double>(dim)));
  sigma_k.push_back(std::vector<std::vector<double>>(dim, std::vector<double>(dim)));

  RcppGSL::matrix<int> z(iteration, data_length + 1);


  std::vector<int> n_k(2);
  int k_count = 0;
  int k_count_used = 0;
  int count = 0;
 
  std::vector<double> prob_k(2);
  RcppGSL::vector<double> sample(dim);
  for(int i = 0;i < dim; i++){
    sample[i] = 0;
  }

  RcppGSL::vector<double> work(dim);
  for(int i = 0;i < dim; i++){
    work[i] = 0;
  }

  Rcout << "____ Simulation Start. ____\n";
  for(int i = 0; i < data_length; i++){
    if(i == 0){
      z(0, 0) = 0;
      gsl_ran_multivariate_gaussian(r, mu_0, sigma_new, sample);
      for(int k = 0; k < dim; k++){
        mu_k[0][k] = sample[k];
        for(int l = 0; l < dim; l++){
          sigma_k[0][k][l] = sigma_new(k, l);
        }
      }
      n_k[0] = 1;
      k_count = 0;
    }
    else{
      for(int j = 0; j <= k_count; j++){
        double prob_k_tmp = 0;
        RcppGSL::vector<double> data_v_in(dim);
        RcppGSL::vector<double> mu_k_v_in(dim);
        for(int k = 0; k < dim; k++){
          data_v_in[k] = data(i, k);
          mu_k_v_in[k] = mu_k[j][k];
        }
        gsl_ran_multivariate_gaussian_pdf(data_v_in, mu_k_v_in, sigma_new, &prob_k_tmp, work);
        data_v_in.free();
        mu_k_v_in.free();
        double prob_k_CRP = n_k[j] / (data_length - 1 + alpha);
        prob_k[j] = prob_k_tmp * prob_k_CRP;
        if(std::isnan(prob_k_tmp)){
          prob_k[j] = 0;
          Rcout << "Warning: "<< j << "th cluster's probability is zero.\n";
        }
      }
      gsl_ran_multivariate_gaussian(r, mu_0, sigma_new, sample);
      for(int j = 0; j < dim; j++){
        mu_k[k_count + 1][j] = sample[j];
      }
      double prob_k_tmp = 0;
      RcppGSL::vector<double> data_v_in(dim);
      RcppGSL::vector<double> mu_k_v_in(dim);
      for(int j = 0;j < dim; j++){
        data_v_in[j] = data(i, j);
      }
      for(int j = 0; j < dim; j++){
        mu_k_v_in[j] = mu_k[k_count + 1][j];
      }
      gsl_ran_multivariate_gaussian_pdf(data_v_in, mu_k_v_in, sigma_new, &prob_k_tmp, work);
      data_v_in.free();
      mu_k_v_in.free();
      double prob_k_CRP = alpha / (data_length - 1 + alpha);
      prob_k[k_count + 1] = prob_k_tmp * prob_k_CRP;
      if(std::isnan(prob_k_tmp)){
        prob_k[k_count + 1] = 0;
        Rcout << "Warning: "<< k_count + 1 << "th cluster's probability is zero.\n";
      }
      double ppp[k_count + 1];
      for(int j = 0; j<= k_count + 1; j++){
        ppp[j] = prob_k[j];
      }
      unsigned int mult_op[k_count + 1];
      gsl_ran_multinomial(r, k_count + 2, 1, ppp,mult_op);

      if(mult_op[k_count + 1] == 1){
        k_count = k_count + 1;
        mu_k.push_back(std::vector<double>(dim));
        sigma_k.push_back(std::vector<std::vector<double>>(dim, std::vector<double>(dim)));
        prob_k.push_back(0);
        n_k.push_back(0);
        if(k_count_used == 0){
          k_count_used = 1;
        }
        for(int k = 0; k < dim; k++){
          mu_k[k_count][k] = sample[k];
          for(int l = 0; l < dim; l++){
            sigma_k[k_count][k][l] = sigma_new(k, l) + sigma_new(k, l) * ((data_length / k_count_used) / data_length);
          }
        }
      }
      for(int j = 0; j <= k_count; j++){
        if(mult_op[j] == 1){
          z(0, i) = j;
          n_k[j] = n_k[j] + 1;
        }
      }
    }
  }
  int t = 0;
  int k_count_used_before = k_count_used;
  k_count_used = 0;
  for(int j = 0; j < dim; j++){
    for(int k = 0; k < dim; k++){
      sigma_new(j, k) = 0;
    }
  }
  for(int j = 0; j <=k_count; j++){
    if(n_k[j] > 3){
      count = 0;
      k_count_used = k_count_used + 1;
      RcppGSL::matrix<double> sigma_v_in(dim, dim);
      RcppGSL::matrix<double> data_k(n_k[j], dim);
      for(int k = 0; k < data_length; k++){
        if(z(t, k) == j){
          RcppGSL::vector_view<double> data_k_row = gsl_matrix_row(data, k);
          gsl_matrix_set_row(data_k, count, data_k_row);
          count = count + 1;
        }
      }
      for(int k = 0; k < dim; k++){
        double col[n_k[j]];
        RcppGSL::vector_view<double> col_view = gsl_matrix_column(data_k, k);
        for(int l = 0; l < n_k[j]; l++){
          col[l] = col_view[l];
        }
        double mean = gsl_stats_mean(col, 1, n_k[j]);
        mu_k[j][k] = mean;
      }
      for(int k = 0; k < dim; k++){
        for(int l = 0; l < dim; l++){
          RcppGSL::vector_view<double> col_view_k = gsl_matrix_column(data_k, k);
          RcppGSL::vector_view<double> col_view_l = gsl_matrix_column(data_k, l);
          double col_k[n_k[j]];
          double col_l[n_k[j]];
          for(int m = 0; m < n_k[j]; m++){
            col_k[m] = col_view_k[m];
            col_l[m] = col_view_l[m];
          }
          double cov = gsl_stats_covariance(col_k, 1,col_l, 1, n_k[j]);
          sigma_v_in(k, l) = cov;
          sigma_v_in(k, l) = sigma_v_in(k, l) +  sigma_v_in(k, l) * (n_k[j] / data_length);
          sigma_new(k, l) = sigma_new(k, l) + sigma_v_in(k, l);
        }
      }
      RcppGSL::vector<double> eval(dim);
      RcppGSL::Matrix evec(dim, dim);
      gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(dim);
      gsl_eigen_symmv(sigma_v_in, eval, evec, w);
      gsl_eigen_symmv_free(w);
      gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC); 
      int flag = 0;
      for (int k = 0; k < dim; k++) {
        if(eval[k] <= 0){
          flag = 1;
          Rcout << "Warning:" << j << "th cluster's covariance-matrix is not positive definite." << "\n";
          for(int l = 0; l < dim; l++){
            for(int m = 0; m < dim; m++){
              if(k_count_used_before == 0){
                k_count_used_before = 1;
              }
              sigma_k[j][l][m] = sigma_new_mem(l, m);
            }
          }
        }
      }
      data_k.free();
      eval.free();
      evec.free();
      if(flag == 0){
        RcppGSL::matrix<double> L(dim, dim);
        double s;
        for (int i_in = 0; i_in < dim; i_in++) {
          for (int j_in = 0; j_in < i_in; j_in++) {
            s = sigma_v_in(i_in,j_in);
            for (int k_in = 0; k_in < j_in; k_in++){
              s = s - L(i_in,k_in) * L(j_in,k_in);
            }
            L(i_in,j_in) = s / L(j_in,j_in);
          }
          s = sigma_v_in(i_in,i_in);
          for (int k_in = 0; k_in < i_in; k_in++){
            double L_in = L(i_in, k_in);
            s = s - pow(L_in, 2.0);
          }
          L(i_in,i_in) = sqrt(s);
        }
        for(int k = 0; k < dim; k++){
          for(int l = k + 1; l < dim; l++){
            L(k, l) = 0;
          }
        }
        for(int k = 0; k < dim; k++){
          for(int l = 0; l < dim; l++){
            sigma_k[j][k][l] = L(k, l);
          }
        }
        L.free();
      }
      sigma_v_in.free();
      data_k.free();
    }

  }
  if(k_count_used == 0){
    k_count_used = 1;
  }
  for(int j = 0; j < dim; j++){
    for(int k = 0; k < dim; k++){
      sigma_new(j, k) = sigma_new(j, k) / k_count_used;
    }
  }
  RcppGSL::vector<double> eval_sigma_new(dim);
  RcppGSL::Matrix evec_sigma_new(dim, dim);
  gsl_eigen_symmv_workspace * w_sigma_new = gsl_eigen_symmv_alloc(dim);
  gsl_eigen_symmv(sigma_new, eval_sigma_new, evec_sigma_new, w_sigma_new);
  gsl_eigen_symmv_free(w_sigma_new);
  gsl_eigen_symmv_sort(eval_sigma_new, evec_sigma_new, GSL_EIGEN_SORT_ABS_ASC); 
  flag = 0;
  for (int k = 0; k < dim; k++) {
    if(eval_sigma_new[k] <= 0){
      flag = 1;
      Rcout << "Warning: sigma_new's covariance-matrix is not positive definite." << "\n";
      for(int l = 0; l < dim; l++){
        for(int m = 0; m < dim; m++){
          sigma_new(l, m) = sigma_new_mem(l, m);
        }
      }
    }
  }
  eval_sigma_new.free();
  evec_sigma_new.free();
  if(flag == 0){
    RcppGSL::matrix<double> L(dim, dim);
    double s;
    for (int i_in = 0; i_in < dim; i_in++) {
      for (int j_in = 0; j_in < i_in; j_in++) {
        s = sigma_new(i_in,j_in);
        for (int k_in = 0; k_in < j_in; k_in++){
          s = s - L(i_in,k_in) * L(j_in,k_in);
        }
        L(i_in,j_in) = s / L(j_in,j_in);
      }
      s = sigma_new(i_in,i_in);
      for (int k_in = 0; k_in < i_in; k_in++){
        double L_in = L(i_in, k_in);
        s = s - pow(L_in, 2.0);
      }
      L(i_in,i_in) = sqrt(s);
    }
    for(int k = 0; k < dim; k++){
      for(int l = k + 1; l < dim; l++){
        L(k, l) = 0;
      }
    }
    for(int k = 0; k < dim; k++){
      for(int l = 0; l < dim; l++){
        sigma_new(k, l) = L(k, l);
      }
    }
    L.free();
  }
  Rcout << "____ Iteration Start. ____\n";
  for(int t = 1; t < iteration; t++){
    for(int i = 0; i < data_length; i++){
      n_k[z(t - 1, i)] = n_k[z(t - 1, i)] - 1;
      for(int j = 0; j <= k_count; j++){
        if(n_k[j] <= 0){
          prob_k[j] = 0;
        }
        else{
          double prob_k_tmp = 0;
          RcppGSL::vector<double> data_v_in(dim);
          RcppGSL::vector<double> mu_k_v_in(dim);
          RcppGSL::matrix<double> sigma_v_in(dim, dim);
          for(int k = 0; k < dim; k++){
            data_v_in[k] = data(i, k);
            mu_k_v_in[k] = mu_k[j][k];
            for(int l =0; l < dim; l++){
              sigma_v_in(k, l) = sigma_k[j][k][l];
            }
          }
          gsl_ran_multivariate_gaussian_pdf(data_v_in, mu_k_v_in, sigma_v_in, &prob_k_tmp, work);
          data_v_in.free();
          mu_k_v_in.free();
          sigma_v_in.free();
          double prob_k_CRP = n_k[j] / (data_length - 1 + alpha);
          prob_k[j] = prob_k_tmp * prob_k_CRP; 
          if(std::isnan(prob_k_tmp)){
            prob_k[j] = 0;
            Rcout << "Warning: "<< j << "th cluster's probability is zero.\n";
          }
        }
      }
      gsl_ran_multivariate_gaussian(r, mu_0, sigma_new, sample);
      double prob_k_tmp = 0;
      RcppGSL::vector<double> data_v_in(dim);
      for(int j = 0; j < dim; j++){
        data_v_in[j] = data(i, j);
      }
      gsl_ran_multivariate_gaussian_pdf(data_v_in, sample, sigma_new, &prob_k_tmp, work);
      data_v_in.free();
      double prob_k_CRP = alpha / (data_length - 1 + alpha);
      prob_k[k_count + 1] = prob_k_tmp * prob_k_CRP;
      if(std::isnan(prob_k_tmp)){
        prob_k[k_count + 1] = 0;
        Rcout << "Warning: "<< k_count + 1 << "th cluster's probability is zero.\n";
      }
      double ppp[k_count +1];
      for(int j = 0; j<= k_count + 1; j++){
        ppp[j] = prob_k[j];
      }
      unsigned int mult_op[k_count + 1];
      gsl_ran_multinomial(r, k_count + 2, 1, ppp,mult_op);
      if(mult_op[k_count + 1] == 1){
        k_count = k_count + 1;
        mu_k.push_back(std::vector<double>(dim));
        sigma_k.push_back(std::vector<std::vector<double>>(dim, std::vector<double>(dim)));
        prob_k.push_back(0);
        n_k.push_back(0);
        if(k_count_used == 0){
          k_count_used = 1;
        }
        for(int k = 0; k < dim; k++){
          mu_k[k_count][k] = sample[k];
          for(int l = 0; l < dim; l++){
            sigma_k[k_count][k][l] = sigma_new(k, l) + sigma_new(k, l) * ((data_length / k_count_used) / data_length);
          }
        }
      }
      for(int j = 0; j <= k_count; j++){
        if(mult_op[j] == 1){
          z(t, i)= j;
          n_k[j] = n_k[j] + 1;
        }
      }
    }
    int k_count_used_before = k_count_used;
    k_count_used = 0;
    for(int j = 0; j < dim; j++){
      for(int k = 0; k < dim; k++){
        sigma_new(j, k) = 0;
      }
    }
    for(int j = 0;j <= k_count; j++){
      if(n_k[j] > 3){
        k_count_used = k_count_used + 1;
        count = 0;
        RcppGSL::matrix<double> sigma_v_in(dim, dim);
        RcppGSL::matrix<double> data_k(n_k[j], dim);
        for(int k = 0; k < data_length; k++){
          if(z(t, k) == j){
            RcppGSL::vector_view<double> data_k_row = gsl_matrix_row(data, k);
            gsl_matrix_set_row(data_k, count, data_k_row);
            count = count + 1;
          }
        }
        for(int k = 0; k < dim; k++){
          double col[n_k[j]];
          RcppGSL::vector_view<double> col_view = gsl_matrix_column(data_k, k);
          for(int l = 0; l < n_k[j]; l++){
            col[l] = col_view[l];
          }
          double mean = gsl_stats_mean(col, 1, n_k[j]);
          mu_k[j][k] = mean;
        }
        for(int k = 0; k < dim; k++){
          for(int l = 0; l < dim; l++){
            RcppGSL::vector_view<double> col_view_k = gsl_matrix_column(data_k, k);
            RcppGSL::vector_view<double> col_view_l = gsl_matrix_column(data_k, l);
            double col_k[n_k[j]];
            double col_l[n_k[j]];
            for(int m = 0; m < n_k[j]; m++){
              col_k[m] = col_view_k[m];
              col_l[m] = col_view_l[m];
            }
            double cov = gsl_stats_covariance(col_k, 1,col_l, 1, n_k[j]);
            sigma_v_in(k, l) = cov;
            sigma_v_in(k, l) = sigma_v_in(k, l) +  sigma_v_in(k, l) * (n_k[j] / data_length);
            sigma_new(k, l) = sigma_new(k, l) + sigma_v_in(k, l);
          }
        }
        RcppGSL::vector<double> eval(dim);
        RcppGSL::Matrix evec(dim, dim);
        gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(dim);
        gsl_eigen_symmv(sigma_v_in, eval, evec, w);
        gsl_eigen_symmv_free(w);
        gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC); 
        int flag = 0;
        for (int k = 0; k < dim; k++) {
          if(eval[k] <= 0){
            flag = 1;
            Rcout << "Warning:" << j << "th cluster's covariance-matrix is not positive definite." << "\n";
            for(int l = 0; l < dim; l++){
              for(int m = 0; m < dim; m++){
                sigma_k[j][l][m] = sigma_new_mem(l, m);
              }
            }
          }
        }
        data_k.free();
        eval.free();
        evec.free();
        if(flag == 0){
          RcppGSL::matrix<double> L(dim, dim);
          double s;
          for (int i_in = 0; i_in < dim; i_in++) {
            for (int j_in = 0; j_in < i_in; j_in++) {
              s = sigma_v_in(i_in,j_in);
              for (int k_in = 0; k_in < j_in; k_in++){
                s = s - L(i_in,k_in) * L(j_in,k_in);
              }
              L(i_in,j_in) = s / L(j_in,j_in);
            }
            s = sigma_v_in(i_in,i_in);
            for (int k_in = 0; k_in < i_in; k_in++){
              double L_in = L(i_in, k_in);
              s = s - pow(L_in, 2.0);
            }
            L(i_in,i_in) = sqrt(s);
          }
          for(int k = 0; k < dim; k++){
            for(int l = k + 1; l < dim; l++){
              L(k, l) = 0;
            }
          }
          for(int k = 0; k < dim; k++){
            for(int l = 0; l < dim; l++){
              sigma_k[j][k][l] = L(k, l);
            }
          }
          L.free();          
        }
        sigma_v_in.free();
      }
    }
    if(k_count_used == 0){
      k_count_used = 1;
    }
    for(int j = 0; j < dim; j++){
      for(int k = 0; k < dim; k++){
        sigma_new(j, k) = sigma_new(j, k) / k_count_used;
      }
    }
    RcppGSL::vector<double> eval(dim);
    RcppGSL::Matrix evec(dim, dim);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(dim);
    gsl_eigen_symmv(sigma_new, eval, evec, w);
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC); 
    int flag = 0;
    for (int k = 0; k < dim; k++) {
      if(eval[k] <= 0){
        flag = 1;
        Rcout << "Warning: sigma_new's covariance-matrix is not positive definite." << "\n";
        for(int l = 0; l < dim; l++){
          for(int m = 0; m < dim; m++){
            sigma_new(l, m) = sigma_new_mem(l, m);
          }
        }
      }
    }
    eval.free();
    evec.free();
    if(flag == 0){
      RcppGSL::matrix<double> L(dim, dim);
      double s;
      for (int i_in = 0; i_in < dim; i_in++) {
        for (int j_in = 0; j_in < i_in; j_in++) {
          s = sigma_new(i_in,j_in);
          for (int k_in = 0; k_in < j_in; k_in++){
            s = s - L(i_in,k_in) * L(j_in,k_in);
          }
          L(i_in,j_in) = s / L(j_in,j_in);
        }
        s = sigma_new(i_in,i_in);
        for (int k_in = 0; k_in < i_in; k_in++){
          double L_in = L(i_in, k_in);
          s = s - pow(L_in, 2.0);
        }
        L(i_in,i_in) = sqrt(s);
      }
      for(int k = 0; k < dim; k++){
        for(int l = k + 1; l < dim; l++){
          L(k, l) = 0;
        }
      }
      for(int k = 0; k < dim; k++){
        for(int l = 0; l < dim; l++){
          sigma_new(k, l) = L(k, l);
        }
      }
      L.free();
    }
    Rcout << "iteration : " << t + 1 << "\n";
  }
  RcppGSL::matrix<int> z_count(data_length, k_count + 1);
  RcppGSL::vector<double> max(data_length);
  RcppGSL::vector<double> max_value(data_length);
  for(int i = 0; i < data_length; i++){
    max[i] = 0;
    max_value[i] = 0;
    for(int j = 0; j <= k_count; j++){
      z_count(i, j) = 0;
    }
    for(int t = burn_in; t < iteration; t++){
      z_count(i, z(t, i)) =z_count(i, z(t, i)) + 1;
    }
    for(int j = 0; j <= k_count; j++){
      if(n_k[j] != 0){
        if(max_value[i] <= z_count(i ,j)){
          max[i] = j;
          max_value[i] = z_count(i ,j);
        }
      }
    }
  }
  RcppGSL::vector<int> n_k_result(k_count + 1);
  for(int i = 0; i <= k_count; i++){
    n_k_result[i] = 0;
  }
  for(int i = 0; i < data_length; i++){
    n_k_result[max[i]] = n_k_result[max[i]] + 1;
  }
  count = 0;
  for(int i = 0; i <= k_count; i++){
    if(n_k_result[i] >= 1){
      count = count + 1;
    }
  }

  Rcpp::NumericMatrix result_arry(count, 2 + dim + dim * dim);
  std::vector<std::vector<double>> mu_result;
  std::vector<std::vector<std::vector<double>>> sigma_result;
  for(int i = 0; i <= k_count; i++){
    mu_result.push_back(std::vector<double>(dim));
    sigma_result.push_back(std::vector<std::vector<double>>(dim, std::vector<double>(dim)));
  }
  for(int i = 0; i <= k_count; i++){
    int n_k_result_in = n_k_result[i];
    if(n_k_result_in >= 1){
      count = 0;
      RcppGSL::matrix<double> data_k(n_k_result_in, dim);
      for(int j = 0; j < data_length; j++){
        if(max[j] == i){
          RcppGSL::vector_view<double> data_k_row = gsl_matrix_row(data, j);
          gsl_matrix_set_row(data_k, count, data_k_row);
          count = count + 1;
        }
      }
      for(int j = 0; j < dim; j++){
        double col[n_k_result_in];
        RcppGSL::vector_view<double> col_view = gsl_matrix_column(data_k, j);
        for(int k = 0; k < n_k_result_in; k++){
          col[k] = col_view[k];
        }
        double mean = gsl_stats_mean(col, 1, n_k_result_in);
        mu_result[i][j] = mean;
      }
      for(int k = 0; k < dim; k++){
        for(int l = 0; l < dim; l++){
          RcppGSL::vector_view<double> col_view_k = gsl_matrix_column(data_k, k);
          RcppGSL::vector_view<double> col_view_l = gsl_matrix_column(data_k, l);
          double col_k[n_k_result_in];
          double col_l[n_k_result_in];
          for(int m = 0; m < n_k_result_in; m++){
            col_k[m] = col_view_k[m];
            col_l[m] = col_view_l[m];
          }
          double cov = gsl_stats_covariance(col_k, 1,col_l, 1, n_k_result_in);
          sigma_result[i][k][l] = cov;
          sigma_result[i][l][k] = cov;
        }
      }
    }
  }

  count = 0;
  for(int i = 0; i <= k_count; i++){
    if(n_k_result[i] >= 1){
      result_arry(count, 0) = i;
      result_arry(count, 1) = n_k_result[i];
      for(int j = 0; j < dim; j++){
        result_arry(count, 2+j) = mu_result[i][j];
      }
      int count_in = 0;
      for(int j = 0; j < dim; j++){
        for(int k = 0; k < dim; k++){
          result_arry(count, 2 + dim + count_in) = sigma_result[i][j][k];
          count_in = count_in + 1;
        }
      }
      count_in = 0;
      count = count + 1;
    }
  }
  List result = List::create(Named("clusters") = result_arry, Named("result") = max, Named("z") = z);
 
  gsl_rng_free(r);
  data.free();
  data_result.free();
  mu_0.free();
  sigma_new.free();
  sample.free();
  work.free();
  max.free();
  max_value.free();
  n_k_result.free();
  z.free();
  z_count.free();
  mu_k.clear();
  sigma_k.clear();
  mu_result.clear();
  sigma_result.clear();
  n_k.clear();
  prob_k.clear();

  double entropy_all = 0;
  for(int i = 0; i < result_arry.nrow(); i++){
    entropy_all = entropy_all + (result_arry(i, 1) / data_length) * log2(result_arry(i ,1) / data_length);
  }
  Rcout << "Entropy All : "<< -1 * entropy_all << "\n";
  
  return result;
  exit(0);
}

