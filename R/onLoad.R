#-----------------------------------------------------------------------Welcome
#               Chinese Restaurant Process Clustering
#                                                                    2018/07/28
#                                                                 Masashi Okada
#                                                      okadaalgorithm@gmail.com
#
#
#
# Description: Chinese Restaurant Process is used in order to compose Dirichlet
# Process. The Clustering which uses Chinese Restaurant Process does not need to
# decide a number of clusters in advance. This algorithm automatically adjusts
# it. And this package calculates ambiguity as entropy.

#---------------------------------------------------------------------------End

#-----------------------------------------------------------------------Library
library(rscala)
library(ggplot2)
library(GGally)
library(mvtnorm)
library(MASS)
#---------------------------------------------------------------------------End
#' Loading Scala Jar file
#' @import rscala
#' @param libname: libname
#' @param pkgname: pkagname
#' @export
.onLoad <- function(libname, pkgname) {
  scalaPackage(c(pkgname, "CRPClustering"))
}

#' Markov chain Monte Carlo methods for CRP Clustering
#' @import grid
#' @import png
#' @import rscala
#' @importFrom utils read.table
#' @param data : a matrix of data for clustering. Row is each data_i and column is dimensions of each data_i.
#' @param mu : a vector of center points of data. If data is 3 dimensions, a vector of 3 elements like c(2,4,2).
#' @param sigma_table : a numeric of CRP variance.
#' @param alpha : a numeric of a CRP concentrate rate.
#' @param ro_0 : a numeric of a CRP mu change rate.
#' @param burn_in :  an iteration integer of burn in.
#' @param iteration : an iteration integer.
#' @return result : an array expresses cluster number and joined data number and cluster's mean and variance matrix.
#' @examples
#' data <- matrix(c(1.8,1.9,2.1,2.5,5.6,5.2,6,6.1), 4, 2)
#' result <- crp_gibbs(
#'                       data,
#'                       mu=c(0,0),
#'                       sigma_table=1,
#'                       alpha=1,
#'                       ro_0=1,
#'                       burn_in=100,
#'                       iteration=1000
#'                      )
#' @export
crp_train<- function(data, mu=c(0,0), sigma_table=1, alpha=1, ro_0=1, burn_in=100, iteration=1000) {
  if(is.matrix(data) == FALSE){
    print("Error: Input data type is not matrix.")
    return(FALSE)
  }
  if(is.vector(mu) == FALSE){
    print("Error: Input mu type is not vector.")
    return(FALSE)
  }
  if(is.numeric(sigma_table) == FALSE){
    print("Error: Input sigma_table type is not numeric.")
    return(FALSE)
  }
  if(sigma_table <= 0){
    print("Error: Input sigma_table is less than or equal to zero.")
    return(FALSE)
  }
  if(is.numeric(alpha) == FALSE){
    print("Error: Input alpha type is not numeric.")
    return(FALSE)
  }
  if(alpha <= 0){
    print("Error: Input alpha is less than or equal to zero.")
    return(FALSE)
  }
  if(is.numeric(ro_0) == FALSE){
    print("Error: Input ro_0 type is not numeric.")
    return(FALSE)
  }
  if(ro_0 <= 0){
    print("Error: Input ro_0 is less than or equal to zero.")
    return(FALSE)
  }
  if(burn_in <= 0){
    print("Error: Input burn_in is less than or equal to zero.")
    return(FALSE)
  }
  if(iteration <= 0){
    print("Error: Input iteration is less than or equal to zero.")
    return(FALSE)
  }
  if(ncol(data) != length(mu)){
    print("Error: Input data dimension do not equal to input mu dimention.")
    return(FALSE)
  }

  data_length <- nrow(data)
  dim <- ncol(data)
  burn_in <- as.integer(burn_in)
  iteration <- as.integer(iteration)

  scala(packages="CRPDensity", assign.name='e')
  out <- e(data_length=data_length, dim=dim, data=data, mu=mu, sigma_table=sigma_table, alpha=alpha, ro_0=ro_0, burn_in=burn_in, iteration=iteration) * 'CRPDensity.CRPDensity.crp_gibbs(data_length, dim, data, mu, sigma_table, alpha, ro_0, burn_in, iteration)'
  close(e)

  result <- data.frame(ClusterNumber = out[ , 1], TotalNumber = out[ , 2] ,ResultMean = out[ , 3 : (3 + dim - 1)], ResultVariance = out[ , (3 + dim) : (3 + dim * 2 + 1)])

  entropy_all <- 0
  for(j in 1 : nrow(result)) {
    entropy_all <- entropy_all + (out[j, 2] / data_length) * log2(out[j ,2] / data_length)
  }
  print(paste("Entropy All : ", -1 * entropy_all))
  print(result)

  return(out)
}
#' Predict The Clusters' Probability that Data Joins in. 
#' @import mvtnorm
#' @import MASS
#' @param data: a matrix of data for clustering. Row is each data_i and column is dimensions of each data_i.
#' @param result: an array expresses the return of crp_train method. 
#' @return predict:an array expresses a joined cluster number and joined probability for each cluster. 
#' @export
crp_predict<- function(data=data,result = result){
  dim <- ncol(data)
  data_length <- nrow(data)
  cluster_length <- nrow(result)

  predict_result <- array(0,dim=c(data_length, cluster_length +1))
  cluster_max <- array(0,dim=c(data_length))
  cluster_max_value <- array(0,dim=c(data_length))
  
  count = 1
  for(i in 1 : data_length){
    for(k in 1 : cluster_length){
      mu <- result[k, 3 : (3 + dim -1)]
      sigma <- matrix(result[k, (3 + dim) : (2 + dim + dim * dim )], nrow = dim, ncol = dim)
      predict_result[i, count + 1] <-  dmvnorm(data[i,], mu, sigma)
      if(cluster_max_value[i] <= predict_result[i, count + 1]){
        predict_result[i, 1] <- count
        cluster_max_value[i] <- predict_result[i, count + 1]
      }
      count = count + 1  
    }
    count <- 1
  }
  return (predict_result)
}

#' CRP Plot Matrix Visualization
#' @import ggplot2
#' @import GGally
#' @param data: a matrix of data for clustering. Row is each data_i and column is dimensions of each data_i.
#' @param predict: an array of return from crp_predict method.
#' @export
crp_plot <- function(data_input = data_input, predict =  predict){
  dim <- ncol(data_input)
  data <- data.frame(data_input,cluster=predict[ , 1])
  p_ <- GGally::print_if_interactive
  map <- ggpairs(data, columns = c(1:dim),
    diag=list(mapping = aes(color = factor(cluster),alpha=0.125), continuous=wrap("barDiag", binwidth=0.3)),
    lower = list(mapping = aes(color = factor(cluster)), continuous = wrap(lower_func)),
    upper = list(mapping = aes(color = factor(cluster), alpha=0.6), continuous = wrap("density", size=.35))
  )
  p_(map)
}

lower_func <- function(data, mapping){
  ggplot(data = data, mapping = mapping)+
  geom_point(size=.3, alpha=0.5)
}
