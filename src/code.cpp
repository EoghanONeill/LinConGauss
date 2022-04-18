
# include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

// # include <math.h>

// using namespace Rcpp;




//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::field<arma::vec> FindShift_cpp(double rho,
                                     arma::mat X,
                                     arma::mat A,
                                     arma::vec b){

  // Rcpp::Rcout << "Line 22. \n" ;

  arma::mat atxb_mat = arma::trans(A) * X ;
  // arma::mat atxb_mat = A.t() * X ;
  atxb_mat.each_col() += b;
  // Rcpp::Rcout << "Line 27. \n" ;

  arma::vec mins_atxb_mat = arma::min(atxb_mat,0).t();   //Min of each col
  // Rcpp::Rcout << "Line 30. \n" ;

  // by default sorted in ascending order

  arma::vec gamma_vec = arma::sort( -1*mins_atxb_mat) ;

  // Rcpp::Rcout << "FindShift_cpp line 36. gamma_vec = " << gamma_vec  << ". \n" ;


  // double N = arma::conv_to<double>::from(X.n_cols);
  int N = X.n_cols;
  double rhoN = rho*N;
  // double rhoN_ind = std::floor(static_cast<double>(rhoN));

  int rhoN_ind = std::floor(rhoN) ;

  // double rhoN_ind = std::floor(0.5);
  // double  tempgam = arma::as_scalar(gamma_vec(rhoN_ind -1)) ;
  //   + arma::as_scalar(gamma_vec(rhoN_ind -1)))
  // Rcpp::Rcout << "Line 46. \n" ;

  double gamma_scalar = arma::as_scalar((gamma_vec(rhoN_ind -1) + gamma_vec(rhoN_ind ))/2);

  arma::uvec inside_bins = (mins_atxb_mat + gamma_scalar >0);
  double rho_hat = double(arma::sum(inside_bins))/double(N);

  //Note: indexed from zero, unlike in R
  arma::uvec inside_inds = arma::find(inside_bins == 1);
  // Rcpp::Rcout << "Line 55. \n" ;


  arma::vec gamscalar_out = {gamma_scalar};
  arma::vec rho_hat_out = {rho_hat};
  arma::vec inside_inds_out = arma::conv_to<arma::vec>::from(inside_inds);
  arma::vec inside_bins_out = arma::conv_to<arma::vec>::from(inside_bins);

  arma::field<arma::vec> ret_f(4);

  ret_f(0)  = gamscalar_out ;
  ret_f(1)  = rho_hat_out;
  ret_f(2)  = inside_inds_out;
  ret_f(3)  = inside_bins_out;


  return(ret_f);

}




//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec xtheta_cpp(double theta,
                     arma::vec x0,
                     arma::vec nu){

  arma::vec xtheta_out = x0*(std::cos(theta))+nu*(std::sin(theta)) ;

  return(xtheta_out);

}







//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat LinESS_cpp(arma::mat A,
                     arma::vec b,
                     int N,
                     arma::vec x0,
                     int nskip){


  // arma::arma_rng::set_seed(value) ;
  arma::arma_rng::set_seed_random();

  if(arma::all(arma::trans(A)* x0 + b > 0)){

    // Rcpp::Rcout << "In domain . \n" ;



  }else{
    Rcpp::Rcout <<  "Error. Initial vector not in domain. \n"  ;



    while(arma::prod(arma::trans(A)* x0 + b > 0) == 0){
      Rcpp::Rcout << "Error. Initial vector not in domain. \n"  ;

      x0 = arma::randn<arma::vec>(x0.n_elem);

    }
  }


  arma::mat X(A.n_rows, N, arma::fill::zeros);

  int M = A.n_cols;

  int Ntotal = N*nskip;

  for(int n=0; n<Ntotal;n++){

    // Rcpp::Rcout << "Line 139. n = " << n  << ". \n" ;

    arma::vec nu = arma::randn<arma::vec>(x0.n_elem);

    arma::vec theta_jhalf_vec_plus(M);
    theta_jhalf_vec_plus.fill(arma::datum::nan);

    arma::vec theta_jhalf_vec_minus(M);
    theta_jhalf_vec_minus.fill(arma::datum::nan);

    for(int j=0; j<M;j++){

      double r_j = arma::as_scalar(arma::sqrt( arma::pow(A.col(j).t() *x0,2) + arma::pow(A.col(j).t() *nu,2)      )    )  ;


      double theta_temp_plus = arma::as_scalar(std::acos( - b(j)/r_j) + 2*arma::atan( A.col(j).t() *nu /( r_j + A.col(j).t() *x0 ) ) ) ;

      double theta_temp_minus = arma::as_scalar(-1*std::acos( - b(j)/r_j) + 2*arma::atan( A.col(j).t() *nu /( r_j + A.col(j).t() *x0 ) ) ) ;

      // double theta_temp_plus2 = arma::as_scalar(arma::acos( - b(j)/( arma::sign(A.col(j).t() *x0)*r_j  )) -
      //   arma::atan( A.col(j).t() *nu /( A.col(j).t() *x0 ) ) ) ;
      //
      // double theta_temp_minus2 = arma::as_scalar(-1*arma::acos( - b(j)/( arma::sign(A.col(j).t() *x0)*r_j  )) -
      //   arma::atan( A.col(j).t() *nu /( A.col(j).t() *x0 ) ) ) ;

      // Rcpp::Rcout << "Line 164. theta_temp_plus = " << theta_temp_plus  << ". \n" ;
      //
      // Rcpp::Rcout << "Line 164. theta_temp_minus = " << theta_temp_minus  << ". \n" ;
      //
      // Rcpp::Rcout << "Line 164. theta_temp_plus2 = " << theta_temp_plus2  << ". \n" ;
      //
      // Rcpp::Rcout << "Line 164. theta_temp_minus2 = " << theta_temp_minus2  << ". \n" ;



      theta_jhalf_vec_plus(j) = theta_temp_plus;

      theta_jhalf_vec_minus(j) = theta_temp_minus;




    }


    arma::vec finite_tplus  = theta_jhalf_vec_plus.elem( arma::find_finite(theta_jhalf_vec_plus) ) ;

    arma::vec finite_tminus  = theta_jhalf_vec_minus.elem( arma::find_finite(theta_jhalf_vec_minus) );

    arma::vec theta_jhalf_vec = arma::join_cols( finite_tplus, finite_tminus) ;


    // Rcpp::Rcout << "Line 189. arma::datum::pi = " << arma::datum::pi  << ". \n" ;

    // Rcpp::Rcout << "Line 189. (theta_jhalf_vec < 0 )*2*arma::datum::pi = " << (theta_jhalf_vec < 0 )*2*arma::datum::pi  << ". \n" ;


    // Rcpp::Rcout << "Line 189. theta_jhalf_vec = " << theta_jhalf_vec  << ". \n" ;

    // theta_jhalf_vec = theta_jhalf_vec + (theta_jhalf_vec < 0 )*2*arma::datum::pi;
    theta_jhalf_vec.elem(arma::find(theta_jhalf_vec < 0)) = theta_jhalf_vec.elem(arma::find(theta_jhalf_vec < 0)) + 2*arma::datum::pi;

    // Rcpp::Rcout << "Line 194. theta_jhalf_vec = " << theta_jhalf_vec  << ". \n" ;


    theta_jhalf_vec = arma::sort(theta_jhalf_vec) ;

    // Rcpp::Rcout << "Line 198. theta_jhalf_vec = " << theta_jhalf_vec  << ". \n" ;


    // double delta_theta = arma::pow(10,-6) *2 * pi;
    double delta_theta = 0.0000001 *2 *arma::datum::pi;

    // Rcpp::Rcout << "Line 189. delta_theta = " << delta_theta  << ". \n" ;

    // Rcpp::Rcout << "Line 189. n = " << n  << ". \n" ;

    arma::vec active_directions(theta_jhalf_vec.n_elem , arma::fill::zeros) ;
    // Rcpp::Rcout << "Line 192. active_directions = " << active_directions << ". \n" ;

    arma::vec theta_active;
    // Rcpp::Rcout << "Line 192. theta_jhalf_vec = " << theta_jhalf_vec << ". \n" ;

    if(theta_jhalf_vec.n_elem   != 0){

      for(unsigned int l=0; l< theta_jhalf_vec.n_elem ;l++){
        // Rcpp::Rcout << "Line 199. (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b) = " << (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b) << ". \n" ;
        // Rcpp::Rcout << "Line 199. arma::prod(  (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b) >0   ) = " << arma::prod(  (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b) >0   ) << ". \n" ;
        // Rcpp::Rcout << "Line 199. arma::prod( (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b) >0   ) = " << arma::prod( (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b) >0   ) << ". \n" ;
        // Rcpp::Rcout << "Line 203. arma::prod(  (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b) >0   )  - arma::prod( (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b) >0   ) = " << arma::prod(  (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b) >0   )  - arma::prod( (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b) >0   ) << ". \n" ;

        // Rcpp::Rcout << "Line 192. x0 = " << x0 << ". \n" ;
        // Rcpp::Rcout << "Line 192. nu = " << nu << ". \n" ;
        // Rcpp::Rcout << "Line 192. theta_jhalf_vec(l) = " << theta_jhalf_vec(l) << ". \n" ;

        // Rcpp::Rcout << "Line 192.xtheta_cpp(theta_jhalf_vec(l)  , x0, nu) = " << xtheta_cpp(theta_jhalf_vec(l)  , x0, nu) << ". \n" ;

        // Rcpp::Rcout << "Line 217. arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l)  , x0, nu) +b = " << arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l)  , x0, nu) +b << ". \n" ;

        arma::vec abovevec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b)   ;
        arma::vec belowvec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b)   ;

        // double tempabove = arma::prod(  abovevec   );
        // double tempbelow = arma::prod( belowvec   );


        double tempabove = arma::all(abovevec >0 )  ;
        double tempbelow = arma::all(belowvec >0 )  ;

        // Rcpp::Rcout << "Line 210. abovevec = " << abovevec << ". \n" ;
        // Rcpp::Rcout << "Line 210. belowvec = " << belowvec << ". \n" ;
        //
        // Rcpp::Rcout << "Line 210. tempabove = " << tempabove << ". \n" ;
        // Rcpp::Rcout << "Line 210. tempbelow = " << tempbelow << ". \n" ;

        active_directions(l) = tempabove  - tempbelow;

        // active_directions(l) = arma::prod(  (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b) >0   )  - arma::prod( (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b) >0   );
        // Rcpp::Rcout << "Line 226. active_directions(l) = " << active_directions(l) << ". \n" ;
        // Rcpp::Rcout << "Line 227. l = " << l << ". \n" ;

      } // end loop over l

      theta_active  = theta_jhalf_vec.elem( arma::find(active_directions != 0 ));

      // Rcpp::Rcout << "Line 233. active_directions = " << active_directions << ". \n" ;
      // Rcpp::Rcout << "Line 233. theta_active = " << theta_active << ". \n" ;

      int active_length = theta_active.n_elem;

      int tempcount269 = 0;

      while( (active_length % 2) == 1  ){


        tempcount269 ++;
        if(tempcount269 > 100){
          Rcpp::Rcout << "Line 540. tempcount534 = " << tempcount269  << ". \n" ;
          Rcpp::Rcout << "Line 540. active_length = " << active_length  << ". \n" ;
          Rcpp::Rcout << "Line 540. (active_length % 2) = " << (active_length % 2)  << ". \n" ;

        }

        // delta_theta = arma::pow(10,-1) *2 * pi;
        delta_theta = 0.1 * delta_theta;

        arma::vec new_active_directions(theta_jhalf_vec.n_elem , arma::fill::zeros) ;

        active_directions = new_active_directions;

        // Rcpp::Rcout << "Line 247. active_directions = " << active_directions << ". \n" ;

        for(unsigned int l=0; l< theta_jhalf_vec.n_elem ;l++){
          // arma::uvec abovevec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b) >0  ;
          // arma::uvec belowvec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b) >0  ;
          //
          // double tempabove = arma::prod(  abovevec   );
          // double tempbelow = arma::prod( belowvec   );
          // // Rcpp::Rcout << "Line 210. tempabove = " << tempabove << ". \n" ;
          // // Rcpp::Rcout << "Line 210. tempbelow = " << tempbelow << ". \n" ;
          //
          // active_directions(l) = tempabove  - tempbelow;

          arma::vec abovevec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b)   ;
          arma::vec belowvec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b)   ;

          // double tempabove = arma::prod(  abovevec   );
          // double tempbelow = arma::prod( belowvec   );


          double tempabove = arma::all(abovevec >0 )  ;
          double tempbelow = arma::all(belowvec >0 )  ;

          // Rcpp::Rcout << "Line 210. abovevec = " << abovevec << ". \n" ;
          // Rcpp::Rcout << "Line 210. belowvec = " << belowvec << ". \n" ;
          //
          // Rcpp::Rcout << "Line 210. tempabove = " << tempabove << ". \n" ;
          // Rcpp::Rcout << "Line 210. tempbelow = " << tempbelow << ". \n" ;

          active_directions(l) = tempabove  - tempbelow;


          // active_directions(l) = arma::prod(arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b >0   )  -
          //   arma::prod(arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b >0   );

        } // end loop over l


        arma::vec new_theta_active  = theta_jhalf_vec.elem( arma::find(active_directions != 0 ));

        theta_active = new_theta_active ;

        active_length = theta_active.n_elem;

      } //end while loop


    } // end if statement theta_jhalf_vec.n_elem   != 0


    if(theta_jhalf_vec.n_elem   == 0){
      theta_active = theta_jhalf_vec;
    }

    // Rcpp::Rcout << "Line 282. n = " << n  << ". \n" ;

    bool ellipse_in_domain  = true ;


    if((theta_active.n_elem ==0  ) || (theta_jhalf_vec.n_elem ==0)){
      // Rcpp::Rcout << "Line 288. n = " << n  << ". \n" ;

      theta_active = {0, 2*arma::datum::pi} ;

      // arbitrary angle
      // 2*0.5**arma::datum::pi
      if(arma::prod(arma::trans(A)*xtheta_cpp( arma::datum::pi, x0, nu) +b >0   ) == 0 ){
        ellipse_in_domain = false;

      }


    }else{
      // Rcpp::Rcout << "Line 301. n = " << n  << ". \n" ;
      // Rcpp::Rcout << "Line 301. theta_active = " << theta_active << ". \n" ;
      //
      // Rcpp::Rcout << "Line 301. active_directions = " << active_directions << ". \n" ;
      // Rcpp::Rcout << "Line 301. arma::find(active_directions !=0) = " << arma::find(active_directions !=0) << ". \n" ;

      arma::vec nonzeroactive = active_directions.elem(arma::find(active_directions !=0) );

      double firstactivenum = nonzeroactive(0);

      if(firstactivenum == -1 ){
        // arma::uvec tempinds = arma::regspace(0,  theta_active.n_elem - 1);

        // arma::vec tempactive(1);
        // tempactive(0) = theta_active(0);
        // Rcpp::Rcout << "Line 272. n = " << n  << ". \n" ;

        theta_active = arma::join_cols( theta_active.tail(theta_active.n_elem - 1) ,  theta_active.subvec(0,0)  );


      }

    } //end find active intersections

    // Rcpp::Rcout << "Line 278. n = " << n  << ". \n" ;

    arma::vec slices = theta_active ;

    double rotation_angle = slices(0) ;

    slices = slices - rotation_angle ;

    // arma::vec rotated_slices = slices + (slices < 0)*2*arma::datum::pi ;
    arma::vec rotated_slices = slices ;
    rotated_slices.elem(arma::find(rotated_slices < 0)) = rotated_slices.elem(arma::find(rotated_slices < 0)) + 2*arma::datum::pi;



    arma::mat t_rotated_mat = arma::reshape(rotated_slices, 2, rotated_slices.n_elem/2) ;

    arma::mat rotated_slicemat = arma::trans(t_rotated_mat);


    if(ellipse_in_domain == false){
      throw std::range_error("Line 395 Error. At least one point should be in the domain! \n");

      Rcpp::Rcout <<  "Error. At least one point should be in the domain! \n"  ;

    }
    // Rcpp::Rcout << "Line 297. n = " << n  << ". \n" ;

    arma::vec lengths_temp = rotated_slicemat.col(1) - rotated_slicemat.col(0);

    arma::vec tempzero = arma::zeros<arma::vec>(1) ;
    arma::vec cum_len = arma::join_cols( tempzero , arma::cumsum(lengths_temp) );

    double l_temp = cum_len(cum_len.n_elem - 1);

    //this should probably be edited for:
    // 1. setting a random seed
    // 2. better pseudo random number generation

    double rand_temp = arma::randu();
    double sample_temp = l_temp * rand_temp  ;

    arma::uvec tempcuminds = arma::find(cum_len >= sample_temp);
    arma::uword idx = tempcuminds(0) - 1;


    double theta_u = (rotated_slicemat.col(0))(idx) + sample_temp - cum_len(idx) + rotation_angle;


    arma::vec xtemp = xtheta_cpp(theta_u, x0, nu) ;

    if(nskip ==0){
      X.col(n/nskip) = xtemp ;
    }else{
      if( (n % nskip) == 0){
        X.col(n/nskip) = xtemp ;
      }
    }

    x0 =xtemp ;
    // Rcpp::Rcout << "Line 320. n = " << n  << ". \n" ;

    int tempcount430 = 0;

    while(arma::prod(arma::trans(A)* x0 + b > 0) ==0){
      // Rcpp::Rcout << "Line 406. arma::trans(A)* x0 + b > 0 = " << (arma::trans(A)* x0 + b > 0)  << ". \n" ;
      tempcount430 ++;
      if(tempcount430 > 100){
        Rcpp::Rcout << "Line 540. tempcount430 = " << tempcount430  << ". \n" ;
        // Rcpp::Rcout << "Line 540. active_length = " << active_length  << ". \n" ;
        // Rcpp::Rcout << "Line 540. (active_length % 2) = " << (active_length % 2)  << ". \n" ;

      }


      arma::vec nu = arma::randn<arma::vec>(x0.n_elem);

      arma::vec theta_jhalf_vec_plus(M);
      theta_jhalf_vec_plus.fill(arma::datum::nan);

      arma::vec theta_jhalf_vec_minus(M);
      theta_jhalf_vec_minus.fill(arma::datum::nan);

      for(int j=0; j<M;j++){

        double r_j = arma::as_scalar(arma::sqrt( arma::pow(A.col(j).t() *x0,2) + arma::pow(A.col(j).t() *nu,2)      )    )  ;


        double theta_temp_plus = arma::as_scalar(std::acos( - b(j)/r_j) + 2*arma::atan( A.col(j).t() *nu /( r_j + A.col(j).t() *x0 ) ) ) ;

        double theta_temp_minus = arma::as_scalar(-1*std::acos( - b(j)/r_j) + 2*arma::atan( A.col(j).t() *nu /( r_j + A.col(j).t() *x0 ) ) ) ;

        // theta_temp_plus = arma::as_scalar(arma::acos( - b(j)/( arma::sign(A.col(j).t() *x0)*r_j  )) -
        //   2*arma::atan( A.col(j).t() *nu /( A.col(j).t() *x0 ) ) ) ;
        //
        // theta_temp_minus = arma::as_scalar(-1*arma::acos( - b(j)/( arma::sign(A.col(j).t() *x0)*r_j  )) -
        //   2*arma::atan( A.col(j).t() *nu /( A.col(j).t() *x0 ) ) ) ;

        theta_jhalf_vec_plus(j) = theta_temp_plus;

        theta_jhalf_vec_minus(j) = theta_temp_minus;




      }

      // Rcpp::Rcout << "Line 395. n = " << n  << ". \n" ;

      arma::vec finite_tplus  = theta_jhalf_vec_plus.elem( arma::find_finite(theta_jhalf_vec_plus) ) ;

      arma::vec finite_tminus  = theta_jhalf_vec_minus.elem( arma::find_finite(theta_jhalf_vec_minus) );

      arma::vec theta_jhalf_vec = arma::join_cols( finite_tplus, finite_tminus) ;

      // Rcpp::Rcout << "Line 403. n = " << n  << ". \n" ;

      // theta_jhalf_vec = theta_jhalf_vec + (theta_jhalf_vec < 0 )*2*arma::datum::pi;
      theta_jhalf_vec.elem(arma::find(theta_jhalf_vec < 0)) = theta_jhalf_vec.elem(arma::find(theta_jhalf_vec < 0)) + 2*arma::datum::pi;

      theta_jhalf_vec = arma::sort(theta_jhalf_vec) ;

      // double delta_theta = arma::pow(10,-6) *2 * pi;
      double delta_theta = 0.0000001 *2 *arma::datum::pi;


      // Rcpp::Rcout << "Line 412. n = " << n  << ". \n" ;

      arma::vec active_directions(theta_jhalf_vec.n_elem , arma::fill::zeros) ;

      arma::vec theta_active;

      if(theta_jhalf_vec.n_elem   != 0){

        for(unsigned int l=0; l< theta_jhalf_vec.n_elem ;l++){
          // arma::uvec abovevec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b) >0  ;
          // arma::uvec belowvec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b) >0  ;
          //
          // double tempabove = arma::prod(  abovevec   );
          // double tempbelow = arma::prod( belowvec   );
          // // Rcpp::Rcout << "Line 210. tempabove = " << tempabove << ". \n" ;
          // // Rcpp::Rcout << "Line 210. tempbelow = " << tempbelow << ". \n" ;
          //
          // active_directions(l) = tempabove  - tempbelow;

          arma::vec abovevec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b)   ;
          arma::vec belowvec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b)   ;

          // double tempabove = arma::prod(  abovevec   );
          // double tempbelow = arma::prod( belowvec   );


          double tempabove = arma::all(abovevec >0 )  ;
          double tempbelow = arma::all(belowvec >0 )  ;

          // Rcpp::Rcout << "Line 210. abovevec = " << abovevec << ". \n" ;
          // Rcpp::Rcout << "Line 210. belowvec = " << belowvec << ". \n" ;
          //
          // Rcpp::Rcout << "Line 210. tempabove = " << tempabove << ". \n" ;
          // Rcpp::Rcout << "Line 210. tempbelow = " << tempbelow << ". \n" ;

          active_directions(l) = tempabove  - tempbelow;

          // active_directions(l) = arma::prod(arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b >0   )  -
          //   arma::prod(arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b >0   );

        } // end loop over l

        theta_active  = theta_jhalf_vec.elem( arma::find(active_directions != 0 ));

        // Rcpp::Rcout << "Line 522. n = " << n  << ". \n" ;
        // Rcpp::Rcout << "Line 523. theta_active = " << theta_active  << ". \n" ;
        // Rcpp::Rcout << "Line 524. active_directions = " << active_directions  << ". \n" ;

        int active_length = theta_active.n_elem;

        int tempcount534 = 0;

        while( (active_length % 2) == 1  ){

          tempcount534 ++;
          if(tempcount534 > 100){
            Rcpp::Rcout << "Line 540. tempcount534 = " << tempcount534  << ". \n" ;
            Rcpp::Rcout << "Line 540. active_length = " << active_length  << ". \n" ;
            Rcpp::Rcout << "Line 540. (active_length % 2) = " << (active_length % 2)  << ". \n" ;

          }

          // delta_theta = arma::pow(10,-1) *2 * pi;
          delta_theta = 0.1 * delta_theta;

          arma::vec new_active_directions(theta_jhalf_vec.n_elem , arma::fill::zeros) ;

          active_directions = new_active_directions;

          for(unsigned int l=0; l< theta_jhalf_vec.n_elem ;l++){
            // arma::uvec abovevec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b) >0  ;
            // arma::uvec belowvec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b) >0  ;
            //
            // double tempabove = arma::prod(  abovevec   );
            // double tempbelow = arma::prod( belowvec   );
            // // Rcpp::Rcout << "Line 210. tempabove = " << tempabove << ". \n" ;
            // // Rcpp::Rcout << "Line 210. tempbelow = " << tempbelow << ". \n" ;
            //
            // active_directions(l) = tempabove  - tempbelow;

            arma::vec abovevec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b)   ;
            arma::vec belowvec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b)   ;

            // double tempabove = arma::prod(  abovevec   );
            // double tempbelow = arma::prod( belowvec   );


            double tempabove = arma::all(abovevec >0 )  ;
            double tempbelow = arma::all(belowvec >0 )  ;

            // Rcpp::Rcout << "Line 210. abovevec = " << abovevec << ". \n" ;
            // Rcpp::Rcout << "Line 210. belowvec = " << belowvec << ". \n" ;
            //
            // Rcpp::Rcout << "Line 210. tempabove = " << tempabove << ". \n" ;
            // Rcpp::Rcout << "Line 210. tempbelow = " << tempbelow << ". \n" ;

            active_directions(l) = tempabove  - tempbelow;

            // active_directions(l) = arma::prod(arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) + delta_theta , x0, nu) +b >0   )  -
            //   arma::prod(arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b >0   );

          } // end loop over l

          // Rcpp::Rcout << "Line 569. active_directions = " << active_directions  << ". \n" ;

          arma::vec new_theta_active  = theta_jhalf_vec.elem( arma::find(active_directions != 0 ));

          theta_active = new_theta_active ;

          active_length = theta_active.n_elem;

        } //end while loop


      } // end if statement theta_jhalf_vec.n_elem   != 0

      // Rcpp::Rcout << "Line 581. n = " << n  << ". \n" ;
      // Rcpp::Rcout << "Line 581. theta_active = " << theta_active  << ". \n" ;

      if(theta_jhalf_vec.n_elem   == 0){
        theta_active = theta_jhalf_vec;
      }
      // Rcpp::Rcout << "Line 590. theta_active = " << theta_active  << ". \n" ;

      bool ellipse_in_domain  = true ;


      if((theta_active.n_elem ==0  ) || (theta_jhalf_vec.n_elem ==0)){

        theta_active = {0, 2*arma::datum::pi} ;

        // arbitrary angle
        // 2*0.5**arma::datum::pi
        arma::vec abovevec2 = (arma::trans(A)*xtheta_cpp( arma::datum::pi, x0, nu) +b)   ;
        // arma::vec belowvec = (arma::trans(A)*xtheta_cpp(theta_jhalf_vec(l) - delta_theta , x0, nu) +b)   ;

        // double tempabove = arma::prod(  abovevec   );
        // double tempbelow = arma::prod( belowvec   );


        double tempabove2 = arma::all(abovevec2 >0 )  ;

        if(tempabove2 == 0 ){
          ellipse_in_domain = false;

        }


      }else{
        // Rcpp::Rcout << "Line 608. n = " << n  << ". \n" ;
        // Rcpp::Rcout << "Line 608. theta_active = " << theta_active << ". \n" ;
        //
        // Rcpp::Rcout << "Line 611. active_directions = " << active_directions << ". \n" ;
        // Rcpp::Rcout << "Line 612. arma::find(active_directions !=0) = " << arma::find(active_directions !=0) << ". \n" ;

        arma::vec nonzeroactive = active_directions.elem(arma::find(active_directions !=0) );

        // Rcpp::Rcout << "Line 616. arma::find(active_directions !=0) = " << arma::find(active_directions !=0) << ". \n" ;

        double firstactivenum = nonzeroactive(0);

        if(firstactivenum == -1 ){
          // arma::uvec tempinds = arma::regspace(0,  theta_active.n_elem - 1);

          // arma::vec tempactive(1);
          // tempactive(0) = theta_active(0);
          // Rcpp::Rcout << "Line 625. n = " << n  << ". \n" ;

          theta_active = arma::join_cols( theta_active.tail(theta_active.n_elem - 1) ,  theta_active.subvec(0,0)  );


        }

      } //end find active intersections

      // Rcpp::Rcout << "Line 634. n = " << n  << ". \n" ;

      arma::vec slices = theta_active ;

      double rotation_angle = slices(0) ;

      slices = slices - rotation_angle ;

      // arma::vec rotated_slices = slices + (slices < 0)*2*arma::datum::pi ;
      arma::vec rotated_slices = slices ;
      rotated_slices.elem(arma::find(rotated_slices < 0)) = rotated_slices.elem(arma::find(rotated_slices < 0)) + 2*arma::datum::pi;

      arma::mat t_rotated_mat = arma::reshape(rotated_slices, 2, rotated_slices.n_elem/2) ;

      arma::mat rotated_slicemat = arma::trans(t_rotated_mat);

      // Rcpp::Rcout << "Line 535. n = " << n  << ". \n" ;

      if(ellipse_in_domain == false){
        throw std::range_error("Line 659. Error. At least one point should be in the domain! \n");

        // Rcpp::Rcout <<  "Error. At least one point should be in the domain! \n"  ;

      }

      arma::vec lengths_temp = rotated_slicemat.col(1) - rotated_slicemat.col(0);

      arma::vec tempzero = arma::zeros<arma::vec>(1) ;
      arma::vec cum_len = arma::join_cols( tempzero , arma::cumsum(lengths_temp) );

      double l_temp = cum_len(cum_len.n_elem - 1);

      //this should probably be edited for:
      // 1. setting a random seed
      // 2. better pseudo random number generation


      // Rcpp::Rcout << "Line 686. n = " << n  << ". \n" ;

      double rand_temp = arma::randu();
      double sample_temp = l_temp * rand_temp  ;
      // Rcpp::Rcout << "Line 686. l_temp = " << l_temp  << ". \n" ;

      arma::uvec tempcuminds = arma::find(cum_len >= sample_temp);
      arma::uword idx = tempcuminds(0) - 1;


      double theta_u = (rotated_slicemat.col(0))(idx) + sample_temp - cum_len(idx) + rotation_angle;


      arma::vec xtemp = xtheta_cpp(theta_u, x0, nu) ;


      if(nskip ==0){
        X.col(n/nskip) = xtemp ;
      }else{
        if( (n % nskip) == 0){
          X.col(n/nskip) = xtemp ;
        }
      }

      x0 =xtemp ;


    } // end while loop

  } // end loop over n


  return(X);

}





//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::field<arma::vec>  SubsetSim_cpp(arma::mat A,
                                      arma::vec b,
                                      int N,
                                      double rho,
                                      int nskip){

  // Rcpp::Rcout << "Line 514. \n" ;


  // arma::arma_rng::set_seed(value) ;
  arma::arma_rng::set_seed_random();

  arma::mat X = arma::randn<arma::mat>(A.n_rows , N );

  // Rcpp::Rcout << "Line 521. \n" ;


  arma::field<arma::vec> res_shift = FindShift_cpp(rho,
                                                   X,
                                                   A,
                                                   b);


  // Rcpp::Rcout << "Line 536. \n" ;

  arma::vec gamma_vec  =  res_shift(0)  ;
  arma::vec rho_hat  =  res_shift(1)  ;
  arma::vec inside_inds  =  res_shift(2)  ;

  double logZ =  std::log(rho_hat(0));

  arma::vec shift_seq = gamma_vec;

  // Rcpp::Rcout << "Line 546 gamma_vec = " << gamma_vec  << ". \n" ;

  double gamma_scalar = gamma_vec(0);

  // Rcpp::Rcout << "Line 549. \n" ;

  // Rcpp::Rcout << "Line 551 logZ = " << logZ  << ". \n" ;

  // arma::vec x0 = arma::zeros(A.n_rows);
  arma::vec x0 = X.col(inside_inds(0));


  int tempcount = 0;
  while(gamma_scalar >0){


    tempcount += 1;

    // Rcpp::Rcout << "Line 549. \n" ;

    // Rcpp::Rcout << "Line 565. tempcount = " << tempcount  << ". \n" ;


    bool found_sample = false;

    if(arma::all(arma::trans(A)*x0 + b + gamma_scalar > 0   ) == 1 ){
      found_sample = true;

    }

    // Rcpp::Rcout << "Line 807. \n" ;


    while(found_sample == false){

      x0 = arma::randn(A.n_rows);

      if(arma::all(arma::trans(A)*x0 + b + gamma_scalar > 0   ) == 1 ){
        found_sample = true;

      }

    }

    // Rcpp::Rcout << "Line 821. \n" ;

    arma::mat X = LinESS_cpp(A,
                             b + gamma_scalar,
                             N,
                             x0,
                             nskip);

    // Rcpp::Rcout << "Line 825. \n" ;


    arma::field<arma::vec> res_shift = FindShift_cpp(rho,
                                                     X,
                                                     A,
                                                     b);

    // Rcpp::Rcout << "Line 576. \n" ;

    arma::vec gamma_vec  =  res_shift(0)  ;
    arma::vec rho_hat  =  res_shift(1)  ;
    arma::vec inside_inds  =  res_shift(2)  ;

    x0 = X.col(inside_inds(0));


    logZ = logZ + std::log(rho_hat(0));

    // Rcpp::Rcout << "Line 591. gamma_vec = " << gamma_vec  << ". \n" ;

    gamma_scalar = gamma_vec(0);

    shift_seq = arma::join_cols(shift_seq , gamma_vec.subvec(0,0) )  ;



  }

  // Rcpp::Rcout << "Line 592. \n" ;

  // Rcpp::Rcout << "Line 603.  logZ = " << logZ  << ". \n" ;
  // Rcpp::Rcout << "Line 603.  shift_seq = " << shift_seq  << ". \n" ;

  arma::vec logZ_out = {logZ};

  arma::field<arma::vec> ret_list (2);


  ret_list(0)  = logZ_out ;
  ret_list(1)  = shift_seq ;

  return(ret_list);


}




//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double HDR_algo_cpp(arma::mat A,
                    arma::vec b,
                    arma::vec gamma_vec,
                    int N,
                    int nskip){

  // arma::arma_rng::set_seed(value) ;
  arma::arma_rng::set_seed_random();

  arma::mat X = arma::randn<arma::mat>(A.n_rows , N );

  double logZ = 0;

  for(unsigned int t=0; t< gamma_vec.n_elem ;t++){

    arma::mat atxb_mat = arma::trans(A) * X ;
    // arma::mat atxb_mat = A.t() * X ;
    atxb_mat.each_col() += b;

    arma::vec mins_atxb_mat = arma::min(atxb_mat,0).t();   //Min of each col

    arma::mat L_t = X.cols(arma::find( mins_atxb_mat  + gamma_vec(t) > 0 )) ;

    logZ = logZ + std::log( L_t.n_cols ) - std::log(double(N));

    // aribtrary choice of column?
    arma::vec x0 = L_t.col(0);

    X = LinESS_cpp(A,
                   b + gamma_vec(t),
                   N,
                   x0,
                   nskip);


  } // end of loop over t


  return(logZ);

}



//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double mcint_lincon(int N,
                    arma::mat A,
                    arma::vec b){

  // arma::arma_rng::set_seed(value) ;
  arma::arma_rng::set_seed_random();

  int samps_in = 0;

  for(int n=0; n< N ; n++){

    arma::vec x0 = arma::randn<arma::vec>(A.n_rows);

    if(arma::prod(arma::trans(A)*x0 + b  > 0   ) == 1 ){
      samps_in = samps_in + 1;

    }

  }


  double prob = double(samps_in)/double(N);

  return(prob);

}








