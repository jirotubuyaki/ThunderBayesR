// [[Rcpp::depends(RcppGSL)]]
//------------------------------------------------------------------------------Welcome
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

using namespace Rcpp;

double mean(RcppGSL::matrix<double> &data_in, int i_in){
  double data_sum = 0;
  int count = 0;
  for(int o = 0; o < data_in->size1; o++){
    data_sum = data_sum + data_in(o, i_in);
    count = count + 1; 
  }
  return data_sum / count;
}
double cov(RcppGSL::matrix<double> &data_in, int i_in, int j_in, double i_mean_in, double j_mean_in){
  double data_sum = 0;
  int count = 0;
  for(int o = 0; o < data_in->size1; o++){
    data_sum += (data_in(o, i_in) - i_mean_in) * (data_in(o, j_in) - j_mean_in);
    count = count + 1;
  }
  return data_sum /count;
}
double cal_mean(RcppGSL::matrix<double> &data_in, int dim_in, int j_in, RcppGSL::matrix<int> &z_in,int t_in, int k_in){
  double data_sum = 0;
  int count = 0;
  for(int o = 0; o < data_in->size1; o++){
    if(z_in(t_in, o) == j_in){
      data_sum = data_sum + data_in(o, k_in);
      count = count + 1; 
    }
  }

  return data_sum / count;
}

double cal_cov(RcppGSL::matrix<double> &data_in, int dim_in, int j_in, RcppGSL::matrix<int> &z_in,int t_in, int k_in, int l_in, double k_mean_in, double l_mean_in){
  double data_sum = 0;
  int count = 0;
  for(int o = 0; o < data_in->size1; o++){
    if(z_in(t_in, o) == j_in){
      data_sum += (data_in(o, k_in) - k_mean_in) * (data_in(o, l_in) - l_mean_in);
      count = count + 1;
    }
  }
  return data_sum /count;
}

double cal_mean_result(RcppGSL::matrix<double> &data_in, int dim_in, int j_in, RcppGSL::vector<double> &max_in, int k_in){
  double data_sum = 0;
  int count = 0;
  for(int o = 0; o < data_in->size1; o++){
    if(max_in[o] == j_in){
      data_sum = data_sum + data_in(o, k_in);
      count = count + 1; 
      }
   }
  return data_sum / count;
}

double cal_cov_result(RcppGSL::matrix<double> &data_in, int dim_in, int j_in, RcppGSL::vector<double> &max_in, int k_in, int l_in, double k_mean_in, double l_mean_in){
  double data_sum = 0;
  int count = 0;
  for(int o = 0; o< data_in->size1; o++){
    if(max_in[o] == j_in){
      data_sum += (data_in(o, k_in) - k_mean_in) * (data_in(o, l_in) - l_mean_in);
      count = count + 1;
    }
  }
  return data_sum / count;
}

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
  RcppGSL::matrix<double> sigma_new(dim, dim);
  for(int i = 0; i < dim; i++){
    double i_mean = mean(data, i);
    mu_0[i] = i_mean;
    for(int j = 0; j < dim; j++){
      double i_mean = mean(data, i);
      double j_mean = mean(data, j);
      sigma_new(i, j) = cov(data, i, j, i_mean, j_mean);
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
      Rcout << "Warning: the covariance-matrix is not positive definite." << "\n";
      for(int l = 0; l < dim; l++){
        for(int m = 0; m < dim; m++){
          sigma_new(l, m) = 1;
        }
      }
    }
  }
  eval.free();
  evec.free();
  if(flag == 0){
    gsl_linalg_cholesky_decomp(sigma_new);
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
        if(isnan(prob_k_tmp)){
          prob_k[j] = 0;
          Rcout << "Warning: Joined Probability is 0.\n";
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
      if(isnan(prob_k_tmp)){
        prob_k[k_count + 1] = 0;
        Rcout << "Warning: Joined Probability is 0.\n";
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
            sigma_k[k_count][k][l] = sigma_new(k, l) / k_count_used + (sigma_new(k, l) / k_count_used) * ((data_length / k_count_used) / data_length);
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
      k_count_used = k_count_used + 1;
      RcppGSL::matrix<double> sigma_v_in(dim, dim);
      for(int k = 0;k < dim; k++){
        double k_mean = cal_mean(data,dim, j, z, t, k);
        mu_k[j][k] = k_mean;
        for(int l = 0; l < dim; l++){
          double l_mean = cal_mean(data, dim, j, z, t, l);
          sigma_v_in(k, l) = cal_cov(data, dim, j, z, t, k, l, k_mean, l_mean);
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
          Rcout << "Warning: the covariance-matrix is not positive definite." << "\n";
          for(int l = 0; l < dim; l++){
            for(int m = 0; m < dim; m++){
              if(k_count_used_before == 0){
                k_count_used_before = 1;
              }
              sigma_k[j][l][m] = sigma_new(l, m) / k_count_used_before + (sigma_new(l, m) / k_count_used_before) * ((data_length / k_count_used_before) / data_length);
            }
          }
        }
      }
      eval.free();
      evec.free();
      if(flag == 0){
        gsl_linalg_cholesky_decomp(sigma_v_in);
        for(int k = 0; k < dim; k++){
          for(int l = 0; l < dim; l++){
            sigma_k[j][k][l] = sigma_v_in(k, l);
          }
        }
      }
      sigma_v_in.free();
    }
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
          if(isnan(prob_k_tmp)){
            prob_k[j] = 0;
            Rcout << "Warning: Joined Probability is 0.\n";
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
      if(isnan(prob_k_tmp)){
        prob_k[k_count + 1] = 0;
        Rcout << "Warning: Joined Probability is 0.\n";
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
            sigma_k[k_count][k][l] = sigma_new(k, l) / k_count_used + (sigma_new(k, l) / k_count_used) * ((data_length / k_count_used) / data_length);
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
        RcppGSL::matrix<double> sigma_v_in(dim, dim);
        for(int k = 0; k < dim; k++){
          double k_mean = cal_mean(data, dim, j, z, t, k);
          mu_k[j][k] = k_mean;
          for(int l = 0; l < dim; l++){
            double l_mean = cal_mean(data, dim, j, z, t, l);
            sigma_v_in(k, l) = cal_cov(data, dim, j, z, t, k, l, k_mean, l_mean);
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
            Rcout << "Warning: the covariance-matrix is not positive definite." << "\n";
            for(int l = 0; l < dim; l++){
              for(int m = 0; m < dim; m++){
                if(k_count_used_before == 0){
                  k_count_used_before = 1;
                }
                sigma_k[j][l][m] = sigma_new(l, m) / k_count_used_before + (sigma_new(l, m) / k_count_used_before) * ((data_length / k_count_used_before) / data_length);
              }
            }
          }
        }
        eval.free();
        evec.free();
        if(flag == 0){
          gsl_linalg_cholesky_decomp(sigma_v_in);
          for(int k = 0; k < dim; k++){
            for(int l =0; l < dim; l++){
              sigma_k[j][k][l] = sigma_v_in(k, l);
            }
          }
        }
        sigma_v_in.free();
      }
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
  RcppGSL::vector<double> n_k_result(k_count + 1);
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
    if(n_k_result[i] >= 1){
      for(int j = 0; j < dim; j++){
        for(int k = 0; k < dim; k++){
          if(j <= k){
            double j_mean = cal_mean_result(data_result, dim, i, max, j);
            double k_mean = cal_mean_result(data_result, dim, i, max, k);
            mu_result[i][j] = j_mean;
            mu_result[i][k] = k_mean;
            sigma_result[i][j][k] = cal_cov_result(data_result, dim, i, max, j, k, j_mean, k_mean);
            sigma_result[i][k][j] = cal_cov_result(data_result, dim, i, max, j, k, j_mean, k_mean);
          }
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

