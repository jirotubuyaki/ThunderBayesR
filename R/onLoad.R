#-----------------------------------------------------------------------Welcome
#               Chinese Restaurant Process Clustering
#                                                                    2020/04/12
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
library(Rcpp)
library(RcppGSL)
library(ggplot2)
library(GGally)
library(MASS)
library(matrixcalc)
#---------------------------------------------------------------------------End

#' Markov chain Monte Carlo methods for CRP Clustering
#' @import grid
#' @import png
#' @import Rcpp
#' @import RcppGSL
#' @import ggplot2
#' @import GGally
#' @importFrom utils read.table
#' @param data: a data.frame of data for clustering. Row is each data_i and column is dimensions of each data_i.
#' @param alpha : a numeric of a CRP concentrate rate.
#' @param burn_in :  an abandoned iteration integer.
#' @param iteration : an iteration integer.
#' @param plot : a logical type of whether plot a result or not.
#' @return result : a list has three elements. The "clusters" is cluster number and joined data number and cluster's mean and variance matrix. The "max" is the cluster number for data i join in. The "z" is the iteration history for each data i join in paticular cluster.
#' @examples
#' data <- array(0, dim=c(30, 2))
#' data[1, ] <- c(-2.43185475495409,2.03311203531621)
#' data[2, ] <- c(-0.408783769317068,1.75003135626229)
#' data[3, ] <- c(1.55675133967144,-1.56420523659905)
#' data[4, ] <- c(0.657368060264562,-0.97383866031164)
#' data[5, ] <- c(0.068212889245427,1.13418709295999)
#' data[6, ] <- c(-1.73561815639666,1.81918235905252)
#' data[7, ] <- c(1.03672457449158,-1.8569658734557)
#' data[8, ] <- c(1.76100672727452,0.0761038984873157)
#' data[9, ] <- c(1.65809552931208,-1.68283298969143)
#' data[10, ] < c(-1.07764013453211,-0.32226716632212)
#' data[11, ] <- c(-0.261606415434224,1.75394513146391)
#' data[12, ] <- c(-0.70752394485502,-0.259888201772335)
#' data[13, ] <- c(-1.38617236009739,-1.24305620615396)
#' data[14, ] <- c(1.51663782772836,0.0161396983547837)
#' data[15, ] <- c(-0.914042151747574,1.76495756094862)
#' data[16, ] <- c(0.282796296711756,-0.0492279088948679)
#' data[17, ] <- c(1.08831769386705,-0.954851525684512)
#' data[18, ] <- c(-0.932745904717591,0.762387679372797)
#' data[19, ] <- c(1.69617665862324,-1.11969182200371)
#' data[20, ] <- c(-0.767903781929682,1.19342570049071)
#' data[21, ] <- c(0.401780436172324,-0.457100405625718)
#' data[22, ] <- c(0.98090840279728,-0.597487493647301)
#' data[23, ] <- c(-1.29713416781756,0.765401326146141)
#' data[24, ] <- c(-1.7620625656782,2.02686171889867)
#' data[25, ] <- c(1.68478051448753,-0.806918914294913)
#' data[26, ] <- c(0.466622923887724,0.197650126048092)
#' data[27, ] <- c(1.4851646799543,-1.53289806708663)
#' data[28, ] <- c(-2.12370907063395,1.6471140089684)
#' data[29, ] <- c(-0.660332309091402,1.73989289688085)
#' data[30, ] <- c(-0.957127602683051,-0.0156523076019355)
#' data <- round(data, 3)
#' result <- crp_train(
#'                       data,
#'                       alpha=1,
#'                       burn_in=10,
#'                       iteration=100,
#'                       plot=TRUE
#'                      )
#' @export
crp_train<- function(data = data, alpha=1, burn_in=100, iteration=1000, plot=TRUE) {
  if(is.numeric(alpha) == FALSE){
    print("Error: Input alpha type is not numeric.")
    return(FALSE)
  }
  if(alpha <= 0){
    print("Error: Input alpha is less than or equal to zero.")
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
  if(iteration <= burn_in){
    print("Error: Input iteration is more than burn_in.")
    return(FALSE)
  }
  if(is.logical(plot) == FALSE){
    print("Error: Input plot type is not logical.")
    return(FALSE)
  }
  data_result <- as.matrix(data)
  data <- scale(data)
  data <-as.matrix(data)

  data_length <- nrow(data)
  dim <- ncol(data)
  burn_in <- as.integer(burn_in)
  iteration <- as.integer(iteration)
  out <- gibbs(data_length, dim, data,data_result, alpha, burn_in, iteration)
  clusters <- out$clusters
  joined <- out$result
  z <- as.matrix(out$z)
  z_output <- z[burn_in:iteration, 1:data_length]
  result_frame <- data.frame(ClusterNumber = clusters[ , 1], TotalNumber = clusters[ , 2] ,ResultMean = clusters[ , 3 : (3 + dim - 1)], ResultVariance = clusters[ , (3 + dim) : (3 + dim * 2 + 1)])
  output <- list(result = result_frame, z = z_output)

  print(result_frame)
  
  if(plot == TRUE){
    dim <- ncol(data_result)
    data_result_frame <- data.frame(data_result,cluster=joined)
    p_ <- GGally::print_if_interactive
    map <- ggpairs(data_result_frame, columns = c(1:dim),
                   diag=list(mapping = aes(color = factor(cluster),alpha=0.125), continuous=wrap("barDiag", binwidth=0.3)),
                   lower = list(mapping = aes(color = factor(cluster)), continuous = wrap("points", size=0.42,stroke = 0.3)),
                   upper = list(mapping = aes(color = factor(cluster), alpha=0.6), continuous = wrap("density", size=.35))
    )
    p_(map)
    remove(p_)
    remove(map)
    remove(data_result_frame)
  }
  remove(alpha)
  remove(burn_in)
  remove(iteration)
  remove(data_length)
  remove(dim)
  remove(plot)
  remove(out)
  remove(clusters)
  remove(joined)
  remove(z)
  remove(z_output)
  remove(result_frame)
  remove(data)
  remove(data_result)

  return(output)
}
#' Predict The Clusters' Probability that Data Joins in. 
#' @import mvtnorm
#' @import MASS
#' @import matrixcalc
#' @param data: a data.frame of data for clustering. Row is each data_i and column is dimensions of each data_i.
#' @param result: a list returned from crp_train method.
#' @return predict: an array expresses a joined cluster number and joined probability in each cluster. 
#' @export
crp_predict<- function(data = data, result = result){
  data <- as.matrix(data)
  result_result <- result$result
  result_result <- as.matrix(result_result)
  dim <- ncol(data)
  data_length <- nrow(data)
  cluster_length <- nrow(result_result)
  


  predict_result <- array(0,dim=c(data_length, cluster_length + 1))
  cluster_max <- array(0,dim=c(data_length))
  cluster_max_value <- array(0,dim=c(data_length))
  name <- seq(0, length= cluster_length + 1)
  name[1] = "ClusterNumber"
  for(i in 1 : cluster_length){
     name[i + 1] = paste("ClusterProbability-", result_result[i, 1]) 
  }
  colnames(predict_result) <- name 

  count = 1
  for(i in 1 : data_length){
    for(k in 1 : cluster_length){
      mu <- result_result[k, 3 : (3 + dim -1)]
      sigma <- matrix(result_result[k, (3 + dim) : (2 + dim + dim * dim )], nrow = dim, ncol = dim)
      predict_result[i, count + 1] <-  dmvnorm(data[i,], mu, sigma)
      if(cluster_max_value[i] <= predict_result[i, count + 1]){
        predict_result[i, 1] <- result_result[count, 1]
        cluster_max_value[i] <- predict_result[i, count + 1]
      }
      count = count + 1  
    }
    count <- 1
  }
  remove(cluster_max)
  remove(cluster_max_value)
  remove(mu)
  remove(sigma)
  remove(count)
  remove(data)
  remove(result)
  remove(result_result)
  remove(cluster_length)
  remove(data_length)
  remove(dim)
  remove(i)
  remove(k)

  return (predict_result)
}

#' CRP Plot Matrix Visualization
#' @import ggplot2
#' @import GGally
#' @param data: a data.frame of data for clustering. Row is each data_i and column is dimensions of each data_i.
#' @param predict : an array returned from crp_predict method.
#' @export
crp_plot <- function(data = data, predict =  predict){
  data <- as.matrix(data)
  dim <- ncol(data)
  data_frame <- data.frame(data,cluster=predict[ , 1])
  p_plot_ <- GGally::print_if_interactive
  map_plot <- ggpairs(data_frame, columns = c(1:dim),
                 diag=list(mapping = aes(color = factor(cluster), alpha=0.125), continuous=wrap("barDiag", binwidth=0.3)),
                 lower = list(mapping = aes(color = factor(cluster)), continuous = wrap("points", size=0.42,stroke = 0.3)),
                 upper = list(mapping = aes(color = factor(cluster), alpha=0.6), continuous = wrap("density", size=.35))
  )
  p_plot_(map_plot)
  remove(dim)
  remove(data)
  remove(data_frame)
  remove(predict)
  remove(p_plot_)
  remove(map_plot)
}
#' Visualization Probability of the Clusters that Data i Join in.
#' @import ggplot2
#' @param i: a number of each data i.
#' @param result: a list returned from crp_train method. 
#' @export
crp_plot_z <- function(i = 1, result =  result){
  z <- result[["z"]]
  z_i <- matrix(0, nrow=max(z[, i]) + 1 , ncol=1)
  for(j in 1 : nrow(z)){
    z_i[z[j, i] + 1, 1] = z_i[z[j, i] + 1, 1] + 1;
  }
  count <- 1
  for(j in 1 : nrow(z)){
    if(z_i[z[j, i] + 1, 1] != 0){
      count = count + 1
    }
  }
  z_label <- matrix(0, nrow=count, ncol=1)
  z_result <- matrix(0, nrow=count, ncol=1)
  count <- 1
  for(j in 1 : nrow(z_i)){
    if(z_i[j, 1] != 0){
      z_label[count, 1] <- j - 1 
      z_result[count, 1] = z_i[j, 1]
      count <- count + 1
    }
  }
  df <- data.frame(cluster = z_label[1:count, 1], count = z_result[1:count, 1])
  g <- ggplot(df, aes(x = factor(cluster), y = count)) 
  g <- g + geom_bar(stat = "identity", aes(colour=factor(cluster), fill=factor(cluster)))
  plot(g)
  remove(count)
  remove(z_label)
  remove(z_result)
  remove(df)
  remove(g)
  remove(result)
}










