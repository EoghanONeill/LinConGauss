

#' @title Find Shift function implemented in C++
#'
#' @description Find a shift such that a fraction rho of X fall into the resulting domain
#' @import Rcpp
#' @param rho Fraction of X that fall in the resulting domain
#' @param X A matrix. Each column is a sample from a (possibly constrained) multivariate Gaussian.
#' @param A A matrix that together with b defines a set of constraints such that x is in the constrained domain if all(A%*%x + b >0 ). Number of rows is dimension of space in which samples are obtained. number of columns is number of constraints.
#' @param b A vector that together with A defines a set of constraints such that x is in the constrained domain if all(A%*%x + b >0 ). Number of elements is number of constraints.
#' @return A list is returned containing the following elements:
#' \item{gamma_scalar}{Scalar. Shift such that rho*N samples lie in the domain.}
#' \item{rho_hat}{Scalar. True fraction in domain after shift by gamma_scalar. Can deviate from rho.}
#' \item{inside_inds}{Indices of the columns in the constraint after shift by gamma_scalar.}
#' @useDynLib LinConGauss, .registration = TRUE
#' @export
FindShiftFast <- function(rho,X,A,b){

  reslist <- FindShift_cpp(rho, X, A, b)

  ret_list <- list()

  ret_list$gamma_scalar <- reslist[[1]]
  ret_list$rho_hat_out <- reslist[[2]]
  ret_list$inside_inds_out <- reslist[[3]]
  ret_list$inside_bins_out <- reslist[[4]]


  return(ret_list)



}



##///////////////////////////////////////////////////////////////////////////


#' @title Elliptical slice sampling for a linearly constrained standard normal distribution (Algorithm 2 of Gessner et al.) implemented in C++
#'
#' @description Elliptical slice sampling for a linearly constrained standard normal distribution (Algorithm 2 of Gessner et al.).
#' @param A A matrix that together with b defines a set of constraints such that x is in the constrained domain if all(A%*%x + b >0 ). Number of rows is dimension of space in which samples are obtained. number of columns is number of constraints.
#' @param b A vector that together with A defines a set of constraints such that x is in the constrained domain if all(A%*%x + b >0 ). Number of elements is number of constraints.
#' @param N Number of samples to draw from the constrained domain.
#' @param x0 Initial value. Must be a vector that lies inside the constrained domain.
#' @export
#' @return A matrix with number of columns equal to the number of samples, N, and number of rows equal to the the dimension of the space in which samples are obtained.
LinESSFast <- function(A,b,N,x0, nskip = 0){


  retmatrix <- LinESS_cpp(A,
                          b,
                          N,
                          x0,
                          nskip)



  return(retmatrix)



}


##///////////////////////////////////////////////////////////////////////////



#' @title Subset simulation for linear constraints (Algorithm 3 of Gessner et al. 2020) implemented in C++
#'
#' @description Subset simulation for linear constraints (Algorithm 3 of Gessner et al. 2020).
#' @param A A matrix that together with b defines a set of constraints such that x is in the constrained domain if all(A%*%x + b >0 ). Number of rows is dimension of space in which samples are obtained. number of columns is number of constraints.
#' @param b A vector that together with A defines a set of constraints such that x is in the constrained domain if all(A%*%x + b >0 ). Number of elements is number of constraints.
#' @param N Number of samples to draw from the constrained domain.
#' @param rho Fraction of X that fall in the resulting domain
#' @param nskip Number of samples to skip in order to get more independent samples
#' @export
#' @return A list is returned containing the following elements:
#' \item{logZ}{Scalar. Approximation of log of probability (multi-dimensional integral) for the constrained domain. Biased estimate.}
#' \item{shift_seq}{Vector of shifts defining nested domains that contain the constrained domain defined by A and b.}
SubsetSimFast <- function(A,b,N,rho = 0.5,nskip = 0){


  reslist <- SubsetSim_cpp(A,
                           b,
                           N,
                           rho,
                           nskip)


  ret_list <- list()

  ret_list$logZ <- reslist[[1]][1,1]
  ret_list$shift_seq <- reslist[[2]][,1]

  return(ret_list)


}



##///////////////////////////////////////////////////////////////////////////


#' @title The Holmes-Diaconis-Ross algorithm applied to linearly constrained Gaussians (Algorithm 1 of Gessner et al.) implemented in C++.
#'
#' @description The Holmes-Diaconis-Ross algorithm applied to linearly constrained Gaussians (Algorithm 1 of Gessner et al.).
#' @param A A matrix that together with b defines a set of constraints such that x is in the constrained domain if all(A%*%x + b >0 ). Number of rows is dimension of space in which samples are obtained. number of columns is number of constraints.
#' @param b A vector that together with A defines a set of constraints such that x is in the constrained domain if all(A%*%x + b >0 ). Number of elements is number of constraints.
#' @param gamma_vec Vector of shifts defining nested domains that contain the constrained domain defined by A and b.
#' @param N Number of samples to draw from the constrained domain.
#' @export
#' @return Approximation of log of probability (multi-dimensional integral) for the constrained domain.
HDR_algoFast <- function(A, b, gamma_vec, N, nskip = 0){

  logZ <- HDR_algo_cpp(A,
                       b,
                       gamma_vec,
                       N,
                       nskip)

  return(logZ)


}



##///////////////////////////////////////////////////////////////////////////





