package CRPClustering
import breeze.linalg._
import breeze.numerics._
import breeze.stats._
import scala.collection.mutable._
import breeze.stats.distributions.MultivariateGaussian
import breeze.stats.distributions.Multinomial
import org.apache.commons.math3.stat.correlation.Covariance
import org.apache.commons.math3.stat.StatUtils

object CRPClustering extends App{
	def crp_gibbs(data_length: Int, dim: Int, data_raw: Array[Array[Double]], mu: Array[Double], sigma_table: Double, alpha: Double, ro_0: Double, burn_in: Int, iteration: Int) : Array[Array[Double]] = {
		var data = DenseMatrix.zeros[Double](data_length, dim)
		for(i <- 0 to data_length-1){
			for(j <- 0 to dim-1){
				data(i, j) = data_raw(i)(j)
			}
		}
		var mu_0 = DenseVector.zeros[Double](dim)
		for(i <- 0 to dim - 1){
			mu_0(i) = mu(i)
		}

		var sigma_new = DenseMatrix.zeros[Double](dim, dim)
		diag(sigma_new) := sigma_table

		var data_k = DenseMatrix.zeros[Double](2, dim)
		var mu_k = DenseMatrix.zeros[Double](2, dim)	

		var sigma_k: collection.mutable.Map[Int, DenseMatrix[Double]] = collection.mutable.Map(1 -> DenseMatrix.zeros[Double](dim, dim))

		var z = DenseMatrix.zeros[Int](iteration, data_length)
		var n_k: Array[Int] = Array.ofDim[Int](1)
		var entropy_all: Double = 0

		var k_count: Int = 0
		var add_i = Array(Array(Int, Array.ofDim[Double](dim,dim)))
		var add = DenseMatrix.zeros[Double](1, dim)
		var add_v = DenseVector.zeros[Double](1)
		var add_n_k = Array.ofDim[Int](1)

		var prob_k = DenseVector.zeros[Double](2)
		var prob_k_CRP  = DenseVector.zeros[Double](2)

		println("____ Simulation Start. ____")
		for(i <- 0 to data_length - 1){
			if(i == 0){
				z(0, 0) = 0
				var mgauss = new MultivariateGaussian(mu_0, sigma_new)
				var sample = mgauss.sample(1)
				for(k <- 0 to dim - 1){
					mu_k(0, k) = sample(0)(k)
				}
				n_k(0) = 1
				data_k(0, ::) :+= data(0, ::)
				k_count = 0
				sigma_k.update(0, sigma_new)
				for(k <- 0 to dim - 1) mu_k(0, k) = data(0, k)
			}
			else{
				for(j <- 0 to k_count){
					var mgauss = new MultivariateGaussian(mu_k(j, ::).t, sigma_new)
					prob_k(j) = mgauss.pdf(data(i, ::).t)
					prob_k_CRP(j) = n_k(j) / (data_length + alpha)
					prob_k(j) = prob_k(j) * prob_k_CRP(j)
				}
				var mgauss = new MultivariateGaussian(mu_0, sigma_new)
				var sample = mgauss.sample(1)
				for(k <- 0 to dim - 1){
					mu_k(k_count + 1, k) = sample(0)(k)
				}
				var mgauss_new = new MultivariateGaussian(mu_k(k_count + 1, ::).t, sigma_new)
				prob_k(k_count + 1) = mgauss_new.pdf(data(i, ::).t)
				prob_k(k_count + 1) = prob_k(k_count + 1) * (alpha / (data_length + alpha))
				var mult = new Multinomial(prob_k)
				var sample_2 = mult.sample(1)
				if(k_count + 1 == sample_2(0)){
					k_count += 1
					data_k = DenseMatrix.vertcat(data_k, add)
					mu_k = DenseMatrix.vertcat(mu_k, add)
 					sigma_k = sigma_k + (k_count -> DenseMatrix.zeros[Double](dim, dim))
					n_k = Array.concat(n_k, add_n_k)
					prob_k = DenseVector.vertcat(prob_k, add_v)
					prob_k_CRP = DenseVector.vertcat(prob_k_CRP, add_v)
					n_k(k_count) = 0
					sigma_k.update(k_count.toInt, sigma_new)
					for(k <- 0 to dim - 1) mu_k(k_count, k) = data(i, k)
				}
				z(0, i) = sample_2(0)
				for(k <- 0 to dim - 1) data_k(sample_2(0), k) :+= data(i, k)
				n_k(sample_2(0)) += 1
			}
		}
		for(j <- 0 to k_count){
			if(n_k(j) >= dim){
				var data_f: Map[Int, Array[Double]] = Map()
				for(k <- 0 to dim - 1){
					var data_f_row: Array[Double] = Array.ofDim[Double](1)
					var count:Int = 0
					for(m <- 0 to data_length - 1){
						if(z(0, m) == j){
							data_f_row = Array.concat(data_f_row, Array.ofDim[Double](1))
							data_f_row(count) = data(m, k)
							count = count + 1
						}			
					}
					count = 0
					data_f = data_f + (k -> data_f_row)
				}
				var sigma_k_m = DenseMatrix.zeros[Double](dim, dim)
				for(k <- 0 to dim - 1){
					for(l <- 0 to dim - 1){
						var cov = new Covariance()
						var check: Array[Double] = extract_f(data_f.get(k))
						sigma_k_m(k, l) = (n_k(j) / (n_k(j) + ro_0)) * cov.covariance(extract_f(data_f.get(k)), extract_f(data_f.get(l))) + (ro_0 / (n_k(j) + ro_0)) * (mean(extract_f(data_f.get(k))) - mu_0(k)) * (mean(extract_f(data_f.get(l))) - mu_0(l))
					}
				}
				sigma_k.update(j.toInt, sigma_k_m)

				for(k <- 0 to dim - 1){
					mu_k(j, k) = (n_k(j) / (n_k(j) + ro_0)) * (data_k(j, k) / n_k(j)) + (ro_0 / (n_k(j) + ro_0)) * mu_0(k)
				}
			}
		}
		println("____ Iteration Start. ____")
		for(t <- 1 to iteration - 1){
			for(i <- 0 to data_length - 1){
				for(k <- 0 to dim - 1){
					data_k(z(t - 1, i), k) = data_k(z(t - 1, i), k) - data(i, k)
				}
				n_k(z(t - 1, i)) = n_k(z(t - 1, i)) - 1
				for(j <- 0 to k_count){
					if(n_k(j) == 0){
						prob_k(j) = 0
					}
					else{
						var mgauss = new MultivariateGaussian(mu_k(j, ::).t, extract(sigma_k.get(j)))
						prob_k(j) = mgauss.pdf(data(i, ::).t)
						prob_k_CRP(j) = n_k(j) / (data_length - 1 + alpha)
						prob_k(j) = prob_k(j) * prob_k_CRP(j)							
					}
				}
				var mgauss = new MultivariateGaussian(mu_0, sigma_new)
				var sample = mgauss.sample(1)
				for(k <- 0 to dim - 1) mu_k(k_count + 1, k) = sample(0)(k)
				var mgauss_new = new MultivariateGaussian(mu_k(k_count + 1, ::).t, sigma_new)
				prob_k(k_count + 1) = mgauss_new.pdf(data(i, ::).t)
				prob_k(k_count + 1) = prob_k(k_count + 1) * (alpha / (data_length - 1 + alpha))
				var mult = new Multinomial(prob_k)
				var sample_2 = mult.sample(1)

				if(k_count + 1 == sample_2(0)){
					k_count += 1
					data_k = DenseMatrix.vertcat(data_k, add)
					mu_k = DenseMatrix.vertcat(mu_k, add)
 					sigma_k = sigma_k + (k_count -> DenseMatrix.zeros[Double](dim, dim))
					n_k = Array.concat(n_k, add_n_k)
					prob_k = DenseVector.vertcat(prob_k, add_v)
					prob_k_CRP = DenseVector.vertcat(prob_k_CRP, add_v)
					sigma_k.update(k_count.toInt, sigma_new)
					for(k <- 0 to dim - 1) mu_k(k_count, k) = data(i, k)
					n_k(k_count) = 0

				}
				z(t, i) = sample_2(0)
				for(k <- 0 to dim - 1) data_k(sample_2(0), k) :+= data(i, k)
				n_k(sample_2(0)) = n_k(sample_2(0)) + 1
			}
			for(j <- 0 to k_count){
				if(n_k(j) >= dim){
					var data_f: Map[Int, Array[Double]] = Map()
					for(k <- 0 to dim - 1){
						var data_f_row: Array[Double] = Array.ofDim[Double](0)
						var count:Int = 0
						for(l <- 0 to data_length - 1){
							if(z(t, l) == j){
								data_f_row = Array.concat(data_f_row, Array.ofDim[Double](1))
								data_f_row(count) = data(l, k)
								count = count + 1
							}
							
						}
						count = 0
						data_f = data_f + (k -> data_f_row)
					}
					var sigma_k_m = DenseMatrix.zeros[Double](dim, dim)
					for(k <- 0 to dim - 1){
						for(l <- 0 to dim - 1){
							var cov = new Covariance()
							sigma_k_m(k, l) = (n_k(j) / (n_k(j) + ro_0)) * cov.covariance(extract_f(data_f.get(k)), extract_f(data_f.get(l))) + (ro_0 / (n_k(j) + ro_0)) * (mean(extract_f(data_f.get(k))) - mu_0(k)) * (mean(extract_f(data_f.get(l))) - mu_0(l))
						}
					}
					sigma_k.update(j.toInt, sigma_k_m)

					for(k <- 0 to dim - 1){
						mu_k(j, k) = (n_k(j) / (n_k(j) + ro_0)) * (data_k(j, k) / n_k(j)) + (ro_0 / (n_k(j) + ro_0)) * mu_0(k)
					}
				}
			}
			println("Iteration: " + t + " Table_number: " + k_count)
		}
		var z_count = Array.ofDim[Int](data_length, k_count + 1)
		var max = Array.ofDim[Int](data_length)
		var max_value = Array.ofDim[Int](data_length)
		for(i <- 0 to data_length - 1){
			max(i) = 0
			max_value(i) = 0
			for(j <- 0 to k_count){
				z_count(i)(j) = 0
			}
			for(t <- burn_in to iteration - 1){
				z_count(i)(z(t, i)) = z_count(i)(z(t, i)) + 1
			}
			for(j <- 0 to k_count){
				if(n_k(j) != 0){
					if(max_value(i) <= z_count(i)(j)){
						max(i) = j
						max_value(i) = z_count(i)(j)
					}
				}
			}	
		}
		var n_k_result: Array[Int] = Array.ofDim[Int](k_count + 1)
		for(j <- 0 to k_count){
			n_k_result(j) = 0
		}
		for(i <- 0 to data_length - 1){
			n_k_result(max(i)) = n_k_result(max(i)) + 1
		}

		var mu_result: Array[Array[Double]] = Array.ofDim[Double](1, dim)
		var count: Int = 0
		for(j <- 0 to k_count){
			if(n_k_result(j) >= 1){
				if(count != 0) {
					mu_result = Array.concat(mu_result, Array.ofDim[Double](1, dim))
				}
				for(k <- 0 to dim - 1){
					mu_result(count)(k) = mu_k(j, k)
				}
				count = count + 1
			}
		}
		var result_arry: Array[Array[Double]] = Array.ofDim(count, 2 + dim + dim * dim)
		var sigma_result: Map[Int, Array[Array[Double]]] = Map(0 -> Array.ofDim[Double](dim, dim))
		var count_sigma: Int = 0
		for(j <- 0 to k_count){
			if(n_k_result(j) >= 1){
				var data_f: Map[Int, Array[Double]] = Map()
				for(k <- 0 to dim - 1){
					var data_f_row: Array[Double] = Array.ofDim[Double](0)
					var count_data_f:Int = 0
					for(l <- 0 to data_length - 1){
						data_f_row = Array.concat(data_f_row, Array.ofDim[Double](1))
						data_f_row(count_data_f) = data(l, k)
						count_data_f = count_data_f + 1	
					}
					count_data_f = 0
					data_f = data_f + (k -> data_f_row)
				}
				var sigma_k_m = Array.ofDim[Double](dim, dim)
				for(k <- 0 to dim - 1){
					for(l <- 0 to dim - 1){
						var cov = new Covariance()
						sigma_k_m(k)(l) = (n_k_result(j) / (n_k_result(j) + ro_0)) * cov.covariance(extract_f(data_f.get(k)), extract_f(data_f.get(l))) + (ro_0 / (n_k_result(j) + ro_0)) * (mean(extract_f(data_f.get(k))) - mu_0(k)) * (mean(extract_f(data_f.get(l))) - mu_0(l))
					}
				}
				sigma_result = sigma_result + (count_sigma -> sigma_k_m)
				count_sigma = count_sigma + 1
			}
		}
		count = 0
		for(j <- 0 to k_count){ 
			if(n_k_result(j) >= 1){
				result_arry(count)(0) = j
				result_arry(count)(1) = n_k(j)
				for(k <- 0 to dim - 1){
					result_arry(count)(2 + k)	= mu_k(j, k)
				}
				var count_in = 0
				var sigma_k_in = extract(sigma_k.get(j))
				for(k <- 0 to dim - 1){
					for(l <- 0 to dim - 1){	
						result_arry(count)( 2 + dim + count_in) = sigma_k_in(k, l)
						count_in += 1
					}
				}
				count_in = 0
				count += 1
			}
		}
		return result_arry
  	}
	def extract(x: Option[DenseMatrix[Double]]):DenseMatrix[Double] = x match {
		case Some(s) => s
	}
	def extract_f(x: Option[Array[Double]]):Array[Double] = x match {
		case Some(s) => s
	}
	def extract_result(x: Option[Array[Array[Double]]]):Array[Array[Double]] = x match {
		case Some(s) => s
	}
}
