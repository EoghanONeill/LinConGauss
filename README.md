
# LinConGauss

<!-- badges: start -->
<!-- badges: end -->

Samples and integrals for standard multivariate Gaussian distributions under linear domain constraints. This package is an implementaiton of the methods described by Gessner et al. (2019). 

Gessner, A., Kanjilal, O., & Hennig, P. (2020, June). Integrals over Gaussians under linear domain constraints. In International Conference on Artificial Intelligence and Statistics (pp. 2764-2774). PMLR.


## Installation

You can install the development version of LinConGauss like so:

``` r
library(devtools)
install_github("EoghanONeill/LinConGauss")
```

## Example

Example for bivariate Gaussian with 5 constraints.


``` r
library(LinConGauss)


A <- rbind(c(1,1,-1,0,0),
           c(-1,0,0,1,-1))

c2 <- c(0,1)
c1 <- c(1,3)

bvec <- c(0,
          -c1[1],
          c1[2],
          -c2[1],
          c2[2])


N <- 1000
rho <- 0.5

sub_ret <- LinConGauss::SubsetSim(A, bvec,N, rho )

hdr_ret <- LinConGauss::HDR_algo(A, bvec,c(sub_ret$shift_seq[1:(length(sub_ret$shift_seq)-1)],0), N )

#from subset sim
exp(sub_ret$logZ)

#from HDR
exp(hdr_ret)




sub_retFast <- LinConGauss::SubsetSimFast(A, bvec,N, rho )

hdr_retFast <- LinConGauss::HDR_algoFast(A, bvec,c(sub_ret$shift_seq[1:(length(sub_ret$shift_seq)-1)],0), N )

#from subset sim
exp(sub_retFast$logZ)

#from HDR
exp(hdr_retFast)


N_samps <- 1000000
samps_in <- 0

for(i in 1:N_samps){
  x0 <- rnorm(nrow(A))

  if( prod( t(A)%*% x0 + bvec  > 0 )   == 1 ){
    samps_in <- samps_in +1
  }


}

samps_in/N_samps


#exact answer

(pnorm(3)-pnorm(1))*(pnorm(1)-pnorm(0))



```

