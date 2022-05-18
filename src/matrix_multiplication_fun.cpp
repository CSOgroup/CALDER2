
	//http://arma.sourceforge.net/docs.html#stats_fns
	#include "RcppArmadillo.h"
	using namespace Rcpp;
	
	// [[Rcpp::export]]
	arma::mat matrix_multiplication_cpp(arma::mat A, arma::mat B)
	{
		return A*B;
	}

	// [[Rcpp::export]]
	arma::mat matrix_multiplication_sym_cpp(arma::mat A)
	{
		return A*A;
	}

	
	
