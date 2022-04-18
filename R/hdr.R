

#' @title Find Shift function
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
FindShift <- function(rho,X,A,b){

  #element m,n of resulting array is t(a_m)%*%x_n
  #add b_m to each element of m^th row
  #equivalent to adding vector b to each column


  # atxb_mat is a matrix of size ncol(A) by ncol(X) ( = N)
  # print(" X = ")
  # print( X)

  atxb_mat <- sweep(t(A) %*% X,1,b,"+")

  #obtain the minimum for each column
  #i.e. minimum over rows m

  # mins_atxb_mat is a vector of ncol(X) ( = N) when we use apply(,2,)
  # print(" ncol(atxb_mat) = ")
  # print( ncol(atxb_mat))
  #
  # print(" nrow(atxb_mat) = ")
  # print( nrow(atxb_mat))

  mins_atxb_mat <- apply(atxb_mat, 2, min)

  # print(" length(mins_atxb_mat) = ")
  # print( length(mins_atxb_mat))


  #Then take the negative and sort in ascending order
  #these are teh shifts sorted in ascending order

  # gamma_vec is a vector of ncol(X) ( = N)
  # print("mins_atxb_mat = ")
  # print(mins_atxb_mat)
  gamma_vec <- sort(-mins_atxb_mat,decreasing = FALSE)

  #there are faster functions
  #for example, from the Rfast package eachcol.apply
  # eachcol.apply(t(A) %*% X,b,indices = NULL,oper = "+")

  #obtain floor of rho*N
  N <- ncol(X)
  rhoN_ind <- floor(rho*N)


  # print(" ncol(X) = ")
  # print( ncol(X))
  # print(" ncol(atxb_mat) = ")
  # print( ncol(atxb_mat))

  # print("length(gamma_vec) = ")
  # print(length(gamma_vec))
  # print("gamma_vec = ")
  # print(gamma_vec)


  # print("rhoN_ind = ")
  # print(rhoN_ind)


  gamma_scalar <- (gamma_vec[rhoN_ind] + gamma_vec[rhoN_ind+1])/2


  # print("gamma_scalar = ")
  # print(gamma_scalar)

  #inside if shift values greater than zero?
  #or greater than shift? (i.e. gamma?)

  # or mins_atxb_mat + gamma_scalar > 0?

  # rho_hat <- sum(mins_atxb_mat >0)

  #presumably must be number in domain after shift by gamma_scalar
  rho_hat <- sum(mins_atxb_mat + gamma_scalar >0)/N # scalar
  inside_inds <- which(mins_atxb_mat + gamma_scalar >0) # vector of length 0 <= length() < N

  # print("rho_hat = ")
  # print(rho_hat)

  # print("inside_inds = ")
  # print(inside_inds)

  ret_list <- list()


  ret_list$gamma_scalar <- gamma_scalar
  ret_list$rho_hat <- rho_hat
  ret_list$inside_inds <- inside_inds

  return(ret_list)



}



##///////////////////////////////////////////////////////////////////////////

##///////////////////////////////////////////////////////////////////////////

##///////////////////////////////////////////////////////////////////////////

##///////////////////////////////////////////////////////////////////////////



xtheta <- function(theta, x0,nu){
  return(x0*as.numeric(cos(theta))+nu*as.numeric(sin(theta)))
}


#' @title Elliptical slice sampling for a linearly constrained standard normal distribution (Algorithm 2 of Gessner et al.)
#'
#' @description Elliptical slice sampling for a linearly constrained standard normal distribution (Algorithm 2 of Gessner et al.).
#' @param A A matrix that together with b defines a set of constraints such that x is in the constrained domain if all(A%*%x + b >0 ). Number of rows is dimension of space in which samples are obtained. number of columns is number of constraints.
#' @param b A vector that together with A defines a set of constraints such that x is in the constrained domain if all(A%*%x + b >0 ). Number of elements is number of constraints.
#' @param N Number of samples to draw from the constrained domain.
#' @param x0 Initial value. Must be a vector that lies inside the constrained domain.
#' @export
#' @return A matrix with number of columns equal to the number of samples, N, and number of rows equal to the the dimension of the space in which samples are obtained.
LinESS <- function(A,b,N,x0, nskip = 0){

  #initial vector needs to be in domain
  #ensure function not clearly defined

  #Resample if not in domain?
  if(all(t(A)%*% x0 + b > 0)){
    # in domain
    # print("in domain")
  }else{
    stop("Error. Initial vector not in domain")

    # perhaps take random samples until a point is found in the doain
    # ElilipticalSliceSampler python function might do something like this
    # at least if there is not any initial value

    while(!(all(t(A)%*% x0 + b > 0))){
      print("not in domain. Sample x0 again")

      x0 <- rnorm(length(x0))
      # x0 <- xnew
    }
  }


  # missing a while loop?
  # while not self.is_converged()

  #initalize sample array
  X <- matrix(NA,
              nrow = nrow(A),
              ncol = N)


  M <- ncol(A)

  # print("line 120")

  for(n in 1:(N*nskip)){

    # print("line 124")

    nu <- rnorm(length(x0))

    #it does not appear to be necessary to define a separate xtheta function

    # print("line 130")

    #
    theta_jhalf_vec_plus <- rep(NA,M)
    theta_jhalf_vec_minus <- rep(NA,M)

    #not clear if require negative and positive values
    # length = ncol(A) or double this

    #CHECK PYTHON CODE
    # print("line 140")

    #in original paper, this is indexed by m
    # m = 1,...,M

    # print("b =")
    # print(b)
    #
    # print("A =")
    # print(A)
    #
    # print("x0 =")
    # print(x0)
    #
    # print("nu =")
    # print(nu)


    for(j in 1:M){
      r_j <- sqrt(  (t(A[,j])%*%x0)^2  + (t(A[,j])%*%nu)^2 )

      # print("j =")
      # print(j)
      #
      # print("r_j =")
      # print(r_j)
      #
      # print("-  b[j]/ r_j =")
      # print(-  b[j]/ r_j)
      #
      # print("(t(A[,j])%*%nu =")
      # print(t(A[,j])%*%nu)
      #
      # print("r_j +  (t(A[,j])%*%x0)  =")
      # print(r_j +  (t(A[,j])%*%x0) )
      #
      # print("acos(-  b[j]/ r_j) =")
      # print(acos(-  b[j]/ r_j))
      #
      # print("atan( (t(A[,j])%*%nu)/ (r_j +  (t(A[,j])%*%x0) )) =")
      # print(atan( (t(A[,j])%*%nu)/ (r_j +  (t(A[,j])%*%x0) )))


      theta_jhalf_vec_plus[j] <- acos(-  b[j]/ r_j) +
        2*atan( (t(A[,j])%*%nu)/ (r_j +  (t(A[,j])%*%x0) ))

      theta_jhalf_vec_minus[j] <- -acos(-  b[j]/ r_j) +
        2*atan( (t(A[,j])%*%nu)/ (r_j +  (t(A[,j])%*%x0) ))

      #just keep finite values ?
      #if theta is less than zero, add 2 pi? before sort

      # print("xtheta(theta =  theta_jhalf_vec_plus[j], x0 = x0, nu = nu ) = ")
      # print(xtheta(theta =  theta_jhalf_vec_plus[j], x0 = x0, nu = nu ))
      #
      # print("t(A[,j]) = ")
      # print(t(A[,j]))


      #test that it actually is a solution
      # print("t(A[,j]) %*% xtheta(theta =  theta_jhalf_vec_plus[j], x0 = x0, nu = nu ) + b  = ")
      # print(t(A[,j]) %*% xtheta(theta =  theta_jhalf_vec_plus[j], x0 = x0, nu = nu ) + b[j]  )
      #
      # print("t(A[,j]) %*% xtheta(theta =  theta_jhalf_vec_minus[j], x0 = x0, nu = nu ) + b  = ")
      # print(t(A[,j]) %*% xtheta(theta =  theta_jhalf_vec_minus[j], x0 = x0, nu = nu ) + b[j]  )


      theta_temp_plus <- acos(-  b[j]/( sign(t(A[,j])%*%x0)*r_j) ) - atan( - (t(A[,j])%*%nu)/ ( (t(A[,j])%*%x0) ))

      theta_temp_minus <- -acos(-  b[j]/( sign(t(A[,j])%*%x0)*r_j) ) - atan( - (t(A[,j])%*%nu)/ ( (t(A[,j])%*%x0) ))

      # print("theta_temp_plus = ")
      #
      # print(theta_temp_plus)
      #
      #
      # print("theta_temp_minus = ")
      #
      # print(theta_temp_minus)


      # print("t(A[,j]) %*% xtheta(theta =  theta_temp_plus, x0 = x0, nu = nu ) + b[j]  = ")
      # print(t(A[,j]) %*% xtheta(theta =  theta_temp_plus, x0 = x0, nu = nu ) + b[j]  )
      #
      # print("t(A[,j]) %*% xtheta(theta =  theta_temp_minus, x0 = x0, nu = nu ) + b[j]  = ")
      # print(t(A[,j]) %*% xtheta(theta =  theta_temp_minus, x0 = x0, nu = nu ) + b[j]  )




      #python code uses sort,
      #maybe  M by 2, and sorts each row (pair)
      # no, actually creates one vector first with is.finite


    }

    # print("line 231")
    # print("theta_jhalf_vec_plus = ")
    # print(theta_jhalf_vec_plus)
    #
    # print("theta_jhalf_vec_minus = ")
    # print(theta_jhalf_vec_minus)

    theta_jhalf_vec <- c(theta_jhalf_vec_plus[is.finite(theta_jhalf_vec_plus)],
                         theta_jhalf_vec_minus[is.finite(theta_jhalf_vec_minus)])

    # print("theta_jhalf_vec = ")
    # print(theta_jhalf_vec)

    # print("line 244")

    # theta_jhalf_vec <- theta_jhalf_vec[is.finite(theta_jhalf_vec)]

    #not sure how precise pi needs to be
    theta_jhalf_vec <- theta_jhalf_vec + (theta_jhalf_vec <0)*2*pi



    theta_jhalf_vec <-  sort(theta_jhalf_vec, decreasing = FALSE)

    # print("theta_jhalf_vec = ")
    # print(theta_jhalf_vec)

    #ell(x) is prod( t(A) %*% x + b > 0 )

    #the above lines equivalent to self.intersection_angles()
    #line 47 in active_intersections.py



    #begin find active intersections
    delta_theta <- 10^{-6} *2 * pi

    #index angles

    # compute indices of angles on the ellipse that are on
    # the boundary  of the integration domain

    #this for loop implements line 49 of active_intersections.py

    # if((length(theta_jhalf_vec) ==0)){
    #   stop("length(theta_jhalf_vec) ==0) , all NaNs.")
    # }


    active_directions <- rep(0, length(theta_jhalf_vec))


    if(length(theta_jhalf_vec) != 0){
      #find active intersections


      for(l in 1:length(theta_jhalf_vec)){

        # print("xtheta( theta_jhalf_vec[l] + delta_theta, x0, nu) = ")
        # print(xtheta( theta_jhalf_vec[l] + delta_theta, x0, nu))
        #
        # print("t(A) = ")
        # print(t(A))
        #
        # print("theta_jhalf_vec = ")
        # print(theta_jhalf_vec)
        # print("b = ")
        # print(b)


        # print("theta_jhalf_vec[l] + delta_theta = ")
        #
        # print(theta_jhalf_vec[l] + delta_theta)
        #
        #
        # print("theta_jhalf_vec[l] - delta_theta = ")
        #
        # print(theta_jhalf_vec[l] - delta_theta)


        active_directions[l] <- prod( t(A) %*% xtheta( theta_jhalf_vec[l] + delta_theta, x0, nu) + b > 0 ) -
          prod( t(A) %*% xtheta( theta_jhalf_vec[l] - delta_theta, x0, nu) + b > 0 )


        # print("t(A) %*% xtheta( theta_jhalf_vec[l] + delta_theta, x0, nu) + b = ")
        # print(t(A) %*% xtheta( theta_jhalf_vec[l] + delta_theta, x0, nu) + b)
        #
        #
        # print("t(A) %*% xtheta( theta_jhalf_vec[l] - delta_theta, x0, nu) + b = ")
        # print(t(A) %*% xtheta( theta_jhalf_vec[l] - delta_theta, x0, nu) + b)
        #
        # print("active_directions = ")
        # print(active_directions)
        #
        # print("line 325")


      }

      #not sure about the following, translated from python
      # theta_active <- theta_jhalf_vec[(active_directions != 0) & is.finite(active_directions)]

      theta_active <- theta_jhalf_vec[active_directions != 0]

      # print("theta_active = ")
      # print(theta_active)

      while( (length(theta_active)  %% 2)  == 1){
        delta_theta <- 10^{-1} * delta_theta

        # print("line 341")


        active_directions <- rep(0, length(theta_jhalf_vec))
        for(l in 1:length(theta_jhalf_vec)){


          active_directions[l] <- prod( t(A) %*% xtheta( theta_jhalf_vec[l] + delta_theta, x0, nu) + b > 0 ) -
            prod( t(A) %*% xtheta( theta_jhalf_vec[l] - delta_theta, x0, nu) + b > 0 )

        }

        theta_active <- theta_jhalf_vec[active_directions != 0]

      }

    }#end if statement for not length zero (i.e. not all NaN)


    # if((length(theta_jhalf_vec) ==0)){
    #   stop("length(theta_jhalf_vec) ==0) , all NaNs.")
    # }

    if((length(theta_jhalf_vec) ==0)){
      theta_active <- theta_jhalf_vec
    }

    ellipse_in_domain <- TRUE

    #this part is not clear
    if((length(theta_active) ==0)| length(theta_jhalf_vec) == 0 ) {
      # if((length(theta_active) ==0) ) {
      theta_active <- c(0,2*pi) # perhaps t(c(0,2*pi))

      # print("line 372")
      #
      # print("x0 = ")
      # print(x0)
      #
      # print("nu = ")
      # print(nu)


      #python code uses runif
      #arbitarily using 0,5 here for debugging purposes
      # if(   prod( t(A) %*% xtheta(2 * pi * runif(1) , x0, nu) + b > 0 ) == 0    ){
      if(   prod( t(A) %*% xtheta(2 * pi * 0.5 , x0, nu) + b > 0 ) == 0    ){
        #entire ellipse is outside of the domain

        # print("prod( t(A) %*% xtheta(2 * pi * runif(1) , x0, nu) + b > 0 ) = ")
        # print(prod( t(A) %*% xtheta(2 * pi * 0.5 , x0, nu) + b > 0 ))
        #
        # print("t(A) %*% xtheta(2 * pi * runif(1) , x0, nu) + b = ")
        # print(t(A) %*% xtheta(2 * pi * 0.5 , x0, nu) + b)

        ellipse_in_domain <- FALSE

      }

    }else{ #check python code
      if( (active_directions[active_directions !=0])[1] ==-1  ){

        theta_active <- c(theta_active[2:length(theta_active)] ,  theta_active[1])

      }
    }#end of find active intersections


    # print("theta_active  = ")
    #
    # print(theta_active)

    slices <- theta_active

    # print("line 408")
    #
    # print("slices  = ")
    #
    # print(slices)


    rotation_angle <- slices[1]
    slices <- slices - rotation_angle

    rotated_slices <- slices + (slices < 0)*2*pi

    # print("rotated slices = ")
    #
    # print(rotated_slices)

    #not clear from python code if guaranteed to have only one row
    rotated_slices <- matrix(rotated_slices, ncol = 2, byrow = TRUE)


    # print("rotated slices = ")
    #
    # print(rotated_slices)


    #slice sample (angle sampler) consists of:
    #active_intersections
    #rotation_angle
    #rotates_slices

    if( ellipse_in_domain == FALSE  ){
      stop("At least one point should be in the domain!")
    }




    #slice_sampler.draw_angle
    lengths_temp <- rotated_slices[,2] - rotated_slices[,1]
    cum_len <- c(0,cumsum(lengths_temp))
    l_temp <- cum_len[length(cum_len)] # l_temp = sum(lengths_temp)


    # print("cum_len = ")
    # print(cum_len)
    #
    # print("cum_len[length(cum_len)] = ")
    # print(cum_len[length(cum_len)])


    #line 9 in pseudo code
    sample_temp <- l_temp*runif(1) # random angle #line 9 of algorithm 2, LinESS?

    #which slices are we in?
    #for each value of sample, return the index of the first element of cum_len  greater than or equal
    #then subtract one (to give index of last element less than sample)
    idx <- which(  cum_len >= sample_temp  )[1] - 1

    #line 10
    theta_u <- rotated_slices[idx,1] + sample_temp - cum_len[idx] + rotation_angle
    #end of draw angle
    #might be transformed, i.e. up to and including line 10

    # print("theta_u = ")
    # print("theta_u")

    xtemp <- xtheta(theta_u, x0, nu)


    if(nskip ==0){
      X[,n] <- xtemp
    }else{
      if( (n %% nskip) == 0){
        # print("Column update, n = ")
        # print(n)
        X[,n/nskip] <- xtemp
      }
    }


    x0 <- xtemp

    # print("attempt new x0 = ")
    # print("x0")
    if( prod( t(A)%*% x0 + b > 0 )   == 0 ){
      print("Line 479 New inital value not in domain")
    }

    # print("line 482 new x0 value = ")
    # print(x0)


    # print("line 486")


    while(  prod((t(A)%*% x0 + b > 0)) ==0    ){

      # print("Point outside domain, resample")
      # nu <- rnorm(1)
      nu <- rnorm(length(x0))


      #it does not appear to be necessary to define a separate xtheta function


      #
      theta_jhalf_vec_plus <- rep(NA,M)
      theta_jhalf_vec_minus <- rep(NA,M)

      #not clear if require negative and positive values
      # length = M or double this

      #CHECK PYTHON CODE

      #in original paper, this is indexed by m
      # m = 1,...,M
      for(j in 1:M){
        r_j <- sqrt(  (t(A[,j])%*%x0)^2  + (t(A[,j])%*%nu)^2 )

        theta_jhalf_vec_plus[j] <- acos(-  b[j]/ r_j) +
          2*atan( (t(A[,j])%*%nu)/ (r_j +  (t(A[,j])%*%x0) ))

        theta_jhalf_vec_minus[j] <- -acos(-  b[j]/ r_j) +
          2*atan( (t(A[,j])%*%nu)/ (r_j +  (t(A[,j])%*%x0) ))

        #just keep finite values ?
        #if theta is less than zero, add 2 pi? before sort

        #python code uses sort,
        #maybe  M by 2, and sorts each row (pair) # no, actually creates one vector first with is.finite


      }


      theta_jhalf_vec <- c(theta_jhalf_vec_plus[is.finite(theta_jhalf_vec_plus)],
                           theta_jhalf_vec_minus[is.finite(theta_jhalf_vec_minus)])

      #not sure how precise pi needs to be
      theta_jhalf_vec <- theta_jhalf_vec + (theta_jhalf_vec <0)*2*pi

      theta_jhalf_vec <-  sort(theta_jhalf_vec, decreasing = FALSE)

      #ell(x) is prod( t(A) %*% x + b > 0 )


      #begin find active intersections
      delta_theta <- 10^{-6} *2 * pi

      #index angles

      # compute indices of angles on the ellipse that are on
      # the boundary  of the integration domain

      active_directions <- rep(0, length(theta_jhalf_vec))
      for(l in 1:length(theta_jhalf_vec)){


        active_directions[l] <- prod( t(A) %*% xtheta( theta_jhalf_vec[l] + delta_theta, x0, nu) + b > 0 ) -
          prod( t(A) %*% xtheta( theta_jhalf_vec[l] - delta_theta, x0, nu) + b > 0 )

      }

      #not sure about the following, translated from python
      theta_active <- theta_jhalf_vec[active_directions != 0]

      while( (length(theta_active)  %% 2)  == 1){
        delta_theta <- 10^{-1} * delta_theta

        active_directions <- rep(0, length(theta_jhalf_vec))
        for(l in 1:length(theta_jhalf_vec)){


          active_directions[l] <- prod( t(A) %*% xtheta( theta_jhalf_vec[l] + delta_theta, x0, nu) + b > 0 ) -
            prod( t(A) %*% xtheta( theta_jhalf_vec[l] - delta_theta, x0, nu) + b > 0 )

        }

        theta_active <- theta_jhalf_vec[active_directions != 0]

      }

      if((length(theta_jhalf_vec) ==0)){
        theta_active <- theta_jhalf_vec
      }

      ellipse_in_domain <- TRUE

      #this part is not clear
      if(length(theta_active) ==0){
        theta_active <- c(0,2*pi) # perhaps t(c(0,2*pi))

        # print("line 536")

        # if(   prod( t(A) %*% xtheta(2 * pi * runif(1) , x0, nu) + b > 0 ) == 0    ){
        if(   prod( t(A) %*% xtheta(2 * pi * 0.5 , x0, nu) + b > 0 ) == 0    ){
          #entire ellipse is outside of the domain
          # print("line 541")
          # print("prod( t(A) %*% xtheta(2 * pi * runif(1) , x0, nu) + b > 0 ) = ")
          # print(prod( t(A) %*% xtheta(2 * pi * 0.5 , x0, nu) + b > 0 ))
          #
          # print("t(A) %*% xtheta(2 * pi * runif(1) , x0, nu) + b = ")
          # print(t(A) %*% xtheta(2 * pi * 0.5 , x0, nu) + b)
          ellipse_in_domain <- FALSE

        }

      }else{
        if( (active_directions[active_directions !=0])[1] ==-1  ){

          theta_active <- c(theta_active[2:length(theta_active)] ,  theta_active[1])

        }
      }#end of find active intersections


      slices <- theta_active
      rotation_angle <- slices[1]
      slices <- slices - rotation_angle

      rotated_slices <- slices + (slices < 0)*2*pi


      # print("theta_active = ")
      # print(theta_active)
      #
      # print("slices = ")
      # print(slices)
      #
      # print("rotated_slices = ")
      # print(rotated_slices)

      #not clear from python code if guranteed to have only one row
      rotated_slices <- matrix(rotated_slices, ncol = 2,byrow = TRUE)



      # print("rotated_slices = ")
      # print(rotated_slices)

      #slice sample (angle sampler) consists of:
      #active_intersections
      #rotation_angle
      #rotates_slices

      # print("line 633")

      if( ellipse_in_domain == FALSE  ){
        stop("At least one point should be in the domain!")
      }


      #slice_sampler.draw_angle
      lengths_temp <- rotated_slices[,2] - rotated_slices[,1]
      cum_len <- c(0,cumsum(lengths_temp))
      l_temp <- cum_len[length(cum_len)]

      # print("cum_len = ")
      # print(cum_len)
      #
      # print("cum_len[length(cum_len)] = ")
      # print(cum_len[length(cum_len)])

      sample_temp <- l_temp*runif(1) # random angle #line 9 of algorithm 2, LinESS?


      # print("sample_temp = ")
      # print(sample_temp)

      #which slices are we in?
      #for each value of sample, return the index of the first element of cum_len  greater than or equal
      #then subtract one (to give index of last element less than sample)
      idx <- which(  cum_len >= sample_temp  )[1] - 1


      theta_u <- rotated_slices[idx,1] + sample_temp - cum_len[idx] + rotation_angle
      #end of draw angle
      #might be transformed, i.e. up to and including line 10

      # print("theta_u = ")
      # print(theta_u)

      xtemp <- xtheta(theta_u, x0, nu)


      if(nskip ==0){
        X[,n] <- xtemp
      }else{
        if( (n %% nskip) == 0){
          # print("Column update, n = ")
          # print(n)
          X[,n/nskip] <- xtemp
        }
      }


      x0 <- xtemp # this new initial value should be in the domain

      # print("attempt new x0 = ")
      # print("x0")

      if( prod( t(A)%*% x0 + b > 0 )   == 0 ){
        stop("Line 628. New inital value not in domain")
      }

      # print("line 631 new x0 value = ")
      print(x0)


    }



  }

  return(X)
}

#' @title Subset simulation for linear constraints (Algorithm 3 of Gessner et al. 2020)
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
SubsetSim <- function(A,b,N,rho = 0.5,nskip = 0){
  #A written with dimensions of Gessner et al

  # Ntotal <- N*(nskip+1)

  #N samples??
  X <- matrix(rnorm(nrow(A)*N),
              nrow = nrow(A),
              ncol = N)

  res_shift <- FindShift(rho,X,A,b)

  gamma_scalar <- res_shift$gamma_scalar # scalar
  rho_hat <- res_shift$rho_hat  # scalar
  inside_inds <- res_shift$inside_inds # vector of length 0 <= length() <= N


  # in the python code, log nesting factor appears to be
  #  log(rho_hat) minus log(N)

  # print("rho =  ")
  # print(rho)
  #
  # print("rho_hat =  ")
  # print(rho_hat)


  logZ <- log(rho_hat) #scalar

  shift_seq <- gamma_scalar# scalar





  # x0 <- NA

  x0 <- X[, inside_inds[1]]

  while(gamma_scalar >0){

    # for x0, randomly choose one column of x for values inside constraint?
    # x0 <- X[ , sample(inside_inds,1), drop = FALSE ] # vector of dim nrow(X) by 1

    #create an initial value
    #not sure if this should be in teh other while loop
    found_sample <- FALSE

    if( prod( t(A)%*% x0 + b + gamma_scalar > 0 )   == 1 ){
      found_sample <- TRUE
    }

    while(found_sample == FALSE){
      print("searching for initial value x0")

      x0 <- rnorm(nrow(A))

      if( prod( t(A)%*% x0 + b + gamma_scalar > 0 )   == 1 ){
        found_sample <- TRUE
      }


    }

    # print("found initial value")
    #
    # print("gamma_scalar = ")
    # print(gamma_scalar)
    #
    # print("t(A)%*% x0 + b + gamma_scalar = ")
    # print(t(A)%*% x0 + b + gamma_scalar)
    #
    # print("x0 = ")
    # print(x0)
    # print("gamma_scalar = ")
    # print(gamma_scalar)

    X <- LinESS(A, b + gamma_scalar, N, x0, nskip)

    # print("gamma_scalar = ")
    # print(gamma_scalar)

    res_shift <- FindShift(rho,X,A,b)
    gamma_scalar <- res_shift$gamma_scalar
    rho_hat <- res_shift$rho_hat
    inside_inds <- res_shift$inside_inds

    # print("inside_inds = ")
    # print(inside_inds)
    # print("gamma_scalar = ")
    # print(gamma_scalar)

    x0 <- X[, inside_inds[1]]

    # print("rho =  ")
    # print(rho)
    #
    # print("rho_hat =  ")
    # print(rho_hat)

    logZ <- logZ + log(rho_hat)

    shift_seq <- c(shift_seq, gamma_scalar)


  }


  ret_list <- list()

  ret_list$logZ <- logZ
  ret_list$shift_seq <- shift_seq

  return(ret_list)



}


##///////////////////////////////////////////////////////////////////////////

##///////////////////////////////////////////////////////////////////////////

##///////////////////////////////////////////////////////////////////////////

##///////////////////////////////////////////////////////////////////////////


#' @title The Holmes-Diaconis-Ross algorithm applied to linearly constrained Gaussians (Algorithm 1 of Gessner et al.).
#'
#' @description The Holmes-Diaconis-Ross algorithm applied to linearly constrained Gaussians (Algorithm 1 of Gessner et al.).
#' @param A A matrix that together with b defines a set of constraints such that x is in the constrained domain if all(A%*%x + b >0 ). Number of rows is dimension of space in which samples are obtained. number of columns is number of constraints.
#' @param b A vector that together with A defines a set of constraints such that x is in the constrained domain if all(A%*%x + b >0 ). Number of elements is number of constraints.
#' @param gamma_vec Vector of shifts defining nested domains that contain the constrained domain defined by A and b.
#' @param N Number of samples to draw from the constrained domain.
#' @export
#' @return Approximation of log of probability (multi-dimensional integral) for the constrained domain.
HDR_algo <- function(A, b, gamma_vec, N, nskip = 0){

  X <- matrix(rnorm(nrow(A)*N),
              nrow = nrow(A),
              ncol = N)

  logZ <- 0

  for(t in 1:length(gamma_vec)){

    atxb_mat <- sweep(t(A) %*% X,1,b,"+")

    #obtain the minimum for each column
    #i.e. minimum over rows m
    mins_atxb_mat <- apply(atxb_mat, 2, min)

    L_t <-  X[, mins_atxb_mat + gamma_vec[t] > 0, drop = FALSE]

    #see compute_log_nesting_factor
    #might need to check log(ncol(L_t))
    logZ <- logZ + log(ncol(L_t))    - log(N)

    #not clear how x0 chosen

    #presumably somewhere in lines 40 to 55 of holmes_diaconis_ross.py

    #t^th column of initial values has to be in t^th nesting ?

    #???
    #it should be the t^th column of X_init
    #not clear how to set X_init
    x0 <- L_t[,1]

    # print("x0 = ")
    # print(x0)
    # print("before LinESS in loop iteration number t =")
    # print (t)

    X <- LinESS(A, b + gamma_vec[t], N, x0, nskip)

    # print("end loop iteration number t =")
    # print (t)
    #
    # print("logZ = ")
    # print(logZ)

  }

  return(logZ)
}






