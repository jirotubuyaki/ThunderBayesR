# CRPClustering: An R package for Bayesian Nonparametric Chinese Restaurant Process Clustering with Entropy 
[![Build Status](https://travis-ci.org/jirotubuyaki/Jdmbs.svg?branch=master)](https://travis-ci.org/jirotubuyaki/Jdmbs)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/Jdmbs)](http://cran.r-project.org/package=Jdmbs)
[![](http://cranlogs.r-pkg.org/badges/grand-total/Jdmbs)](http://cran.rstudio.com/web/packages/Jdmbs/index.html)
## Abstract
Clustering is a scientific method which finds the clusters of data. And lots of methods are traditionally researched for long terms. Bayesian Nonparametric is a statistics which can treat models having infinite parameters. Chinese Restaurant Process is used in order to compose Dirichlet Process. The clustering which use Chinese Restaurant Process does not need to decide a number of clusters in advance. This algorithm automatically adjusts it. And this package can calculate clusters in addition to entropy as ambiguity of clusters.

## Introduction
Clustering is a traditional method in order to find clusters of data. And lots of methods are invented for several decades. The most popular method is called as K-mean[@Hartigan1979]. Kmean is a algorithmic way in order to search clusters of data. But its method needs to decide a number of clusters in advance. So if the data is both high dimensions and complex, deciding accurate numbers of clusters is difficult. And normal bayesian methods are too. For that reason, Bayesian Nonparametric methods is gradually important as computers are faster. In this package, we implement a Chinese Restaurant Process Clustering (CRP)[@Pitman1995]. CRP can compose infinite dimensional parameters as Dirichlet Process. It acts like a customers who sit at tables in restaurant and has probability to sit at a new table. As result, Its model always automates clustering. And we add the method which calculates entropy[@Elliott1999] of clusters into this package. It can check ambiguity of result. Then we explain the clustering model and how to use it in detail. finally, an example is explained.

## Background
#### Chinese Restaurant Process
Chinese Restaurant Process is a metaphor looks like customers sit at a table in Chinese restaurant. All customers except for $x_i$ have allready sat at finite tables. New customer $x_i$ will sit at either a table which other customers have allready sat at or a new table. New customer tends to sit at a table which has a number of customers more than other tables. Probability equation is given by    

![equa](./readme_images/equation_1.png "eque")

$n^{\backslash i}_{k}$ express the number of customers at table $k$. And $\alpha$ is a concentration parameter.

#### Markov chain Monte Carlo methods for CRP Clustering
Markov chain Monte Carlo methods (MCMC) is algorithmic methods[@Liu1994] to sample from posterior distribution. If conditional posterior distribution is given by models, it is the best way in order to acquire parameters as posterior distribution. Algorithm for this package is given by  

lots of iterations continue on below.

i) sampling $z_i$ for each i (i = 1,2, ・・・,n)

![equa](./readme_images/equation_2.png "eque")

ii) sampling $u_k$ for each k (k = 1,2, ・・・,∞)

![equa](./readme_images/equation_3.png "eque")

First parts of iterations as burn in have error range. For that reason, burn in parts are abandoned.

#### Cluster Entropy
Entropy express ambiguity of clustering. As the result of simulation, data $x_i$ joins in particular table. From the total numbers $n_k$ of particular table $k$ at the last iteration, probability $p_k$ at each cluster $k$ is calculated. Entropy equation is given by

![equa](./readme_images/equation_4.png "eque")


## Installation
If download from GitHub, you can use devtools by the commands:

```
> library(devtools)
> install_github("jirotubuyaki/CRPClustering")
```

Once the packages are installed, it needs to be made accessible to the current R session by the commands:

```
> library(CRPClustering)
```

For online help facilities or the details of a particular command (such as the function crp_gibbs) you can type:

```
> help(package="CRPClustering")
```

## Method
#### Calculating Methods for CRP Clustering 
  
```
> z_result <- crp_gibbs(path="./data/data.csv", mu=c(0,0), sigma=0.5,
                        sigma_table=12, burn_in=10, iteration=100)

```
  
This method calculates CRP Clustering.  
・ path : a data path and file. data.csv is a comma separate file. And cols is a dimension of data.  
・ mu : a floating-point type of center points of data. If data is 3 dim, a vector of 3 elements like "c(2,4,7)"  
・ sigma : a floating-point type of data variance.  
・ sigma_table : a floating-point type of table position variance.  
・ burn_in : iteration numbers of burn in.  
・ iteration : iteration numbers.   
・ z_result : a vector and express cluster numbers for each data $i$ 

#### Visualization Methods
  
```
> crp_graph_2d(path="./data/data.csv", z_result)
```
  
This method exhibits a two dimensional graph for the method "crp_gibbs".  
・ path : a data path and file. data.csv is a comma separate file. And cols is a dimension of data.  
・ z_result : the output of the method "crp_gibbs". It contains a number of cluster for each data.  

## Example
Data is generated from three Normal distributions and mu_0 = (-1,1) , mu_1 = (-1.3,-1.3) , mu_2 = (1, -1) and sigma_0 = 0.3$ , sigma_1 = 0.02 , sigma_2 = 0.3. The result is plotted as graph and each data joined in any clusters. The graph is given by below. 

![equa](./readme_images/figure_1.png "eque")

Figure 1. CRP Clustering Result

## Conclusions
Chinese Restaurant Process Clustering was implemented and explained how to use it. Computers resources is limited. Computer processing power is a the most important problem. And several improvements are planed. Please send suggestions and report bugs to okadaalgorithm@gmail.com.

## Acknowledgments
This activity would not have been possible without the support of my family and friends. To my family, thank you for lots of encouragement for me and inspiring me to follow my dreams. I am especially grateful to my parents, who supported me all aspects.  

## References
Hartigan, M. A., J. A.; Wong. 1979. “Algorithm as 136: A K-Means Clustering Algorithm,” Journal of the Royal Statistical Society, SeriesC. 28 (1): 100–108. JSTOR 2346830.  
Liu, Jun S. 1994. “The Collapsed Gibbs Sampler in Bayesian Computations with Applications to a Gene Regulation Problem,” Journal of the American Statistical Association 89 (427): 958–966.  
Pitman, Jim. 1995. “Exchangeable and Partially Exchangeable Random Partitions,” Probability Theory and Related Fields 102 (2): 145–158.  
Yngvason, Elliott H. Lieb; Jakob. 1999. “The Physics and Mathematics of the Second Law of Thermodynamics,” Physics Reports Volume:310 Issue:1 1–96.  
