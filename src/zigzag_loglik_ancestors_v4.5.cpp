    // -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

    #include "RcppArmadillo.h"
    using namespace Rcpp;
    using namespace std;
    // #define M_PI  3.141592653589793238462643383280    /* pi */                    

    // [[Rcpp::export]]
    double loglik_lnorm_cpp( double sum_ln1, double sum_ln2, double p, double q ) // several times faster than fitdistr
    {
        if( sum_ln2 < 0 ) cout << "sum_ln2 not valid in loglik_lnorm_cpp\n";
        if( p < 0 ) cout << "p not valid in loglik_lnorm_cpp\n";
        if( q < 0 ) cout << "q not valid in loglik_lnorm_cpp\n";
            
        if( q <=1 ) return 0; // p = sum(x==0): the number of zero-values
        double alpha = 1.0*p/(p+q);
        double mu = sum_ln1/q;

        double sigma2 = sum_ln2/q + pow(mu, 2) - 2*mu*sum_ln1/q; // %: element-wise product 
        if( abs(sigma2) <=1E-10 ) 
        {
            return 0;
        }
        double loglik = -sum_ln1 - 0.5*q*(log(2*M_PI) + log(sigma2)) - (sum_ln2 + q*pow(mu, 2) - 2*mu*sum_ln1)/2/sigma2; 
        if( p==0 ) {return loglik;} else {return p*log(alpha) + q*log(1-alpha) + loglik;}
    }

    // [[Rcpp::export]]
    double loglik_lnorm_cpp_vec( arma::vec vec_values ) //log-likelihood of vec_values
    {
        int p, q;
        int n = vec_values.n_elem;
        if( n < 2 ) return 0; // added on 19-07-2018

        double sum_ln1, sum_ln2;
        arma::vec positive_vec = vec_values.elem(find(vec_values > 0));
        q = positive_vec.n_elem; // the number of positive values
        p = n - q; // the number of zeros
        if( q <= 1 ) return 0;
        
        sum_ln1 = sum(log(positive_vec));
        sum_ln2 = sum(pow(log(positive_vec), 2));
        
        return loglik_lnorm_cpp( sum_ln1, sum_ln2, p, q );
    }
    
    // [[Rcpp::export]]
    arma::mat get_A_len(arma::mat A) // get matrix A_len: A_len := A*(A>0). This is used for computing the number of positive elements in a rectangle region
    {
        int n_row = A.n_rows;
        arma::mat A_len=arma::zeros<arma::mat>(n_row, n_row); // for test
        arma::uvec ids = find(A > 0);
        arma::vec new_values = arma::ones<arma::vec>(ids.n_elem);
        A_len.elem(ids) = new_values;
        return A_len;
    }
    
    // [[Rcpp::export]]
    arma::mat get_A_ln1(arma::mat A) // log(A_ij)
    {
        int n_row = A.n_rows;
        arma::mat A_ln1=arma::zeros<arma::mat>(n_row, n_row); // for test
        arma::uvec ids = find(A > 0);
        arma::vec new_values = log(A.elem(ids));
        A_ln1.elem(ids) = new_values;
        return A_ln1;
    }
    
    // [[Rcpp::export]]
    arma::mat get_A_ln2(arma::mat A) // log(A_ij)^2
    {
        int n_row = A.n_rows;
        arma::mat A_ln2=arma::zeros<arma::mat>(n_row, n_row); // for test
        arma::uvec ids = find(A > 0);
        arma::vec new_values = pow(log(A.elem(ids)), 2);
        A_ln2.elem(ids) = new_values;
        return A_ln2;
    }
 
     // compute the loglik matrix   
    // [[Rcpp::export]]
    arma::mat loglik_lnorm_cpp_mat( arma::mat sum_ln1, arma::mat sum_ln2, arma::mat ps, arma::mat qs ) // several times faster than fitdistr
    {
        int n_row = sum_ln1.n_rows;
        int n_col = sum_ln1.n_cols;
        arma::mat loglik(n_row, n_col);
        for(int i=0; i<n_row; i++)
            for(int j=0; j<n_col; j++)
                loglik(i,j) = loglik_lnorm_cpp( sum_ln1(i,j), sum_ln2(i,j), ps(i,j), qs(i,j) );
        return loglik;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // compute the blockwise sum
    // [[Rcpp::export]] // use the same version name v4 instead of v4_5
    List zigzag_loglik_ancestors_v4_5( arma::mat A, int k, int min_n_bins=2 ){
        int n_row = A.n_rows;
        // int min_n_bins = 2;
        int p, q, max_mid_index, max_mid;
        arma::mat ps, loglik_tmp, loglik;
        arma::vec n_cells;
        arma::mat A_len = get_A_len(A); // A_len = A*(A > 0)
        arma::mat A_ln1 = get_A_ln1(A); // log(A_ij)
        arma::mat A_ln2 = get_A_ln2(A); // log(A_ij)^2
        StringMatrix ancestors(n_row, n_row);
        arma::mat L=arma::zeros<arma::mat>(n_row, n_row); // for test
        Rcpp::List res;

        for( int k= min_n_bins-1; k<=(2*min_n_bins - 2); k++ ) // other values of L are 0
        {
            for( int v=1; v<= (n_row - k); v++ )
            {
                arma::mat tmp_mat=A.submat(arma::span(v-1, v+k-1), arma::span(v-1, v+k-1)); // span(0,1) := 1:2 in R
                arma::vec upper_tri_vec = tmp_mat.elem(find(trimatu(tmp_mat)));
                L(v-1, v-1+k) = loglik_lnorm_cpp_vec( upper_tri_vec );
            }
        } // Checked to be the same in R
        // cout << "Finished initialize L\n";
        
        // initialize the rad_mat as the first off-diagonal values
        arma::mat rad_mat_current_ln1(n_row-1, 1);
        arma::mat rad_mat_current_ln2(n_row-1, 1);
        arma::mat rad_mat_current_len(n_row-1, 1);
        
        for(int i=0; i<(n_row-1); i++) // The first off-diagonal values
        {
            rad_mat_current_ln1(i, 0) = A_ln1(i, i+1);
            rad_mat_current_ln2(i, 0) = A_ln2(i, i+1);
            rad_mat_current_len(i, 0) = A_len(i, i+1);
        }
        
        // initialized to be two vertical cells (2 rows)
        arma::mat vertical_columns_next_ln1(2, n_row-2);
        arma::mat vertical_columns_next_ln2(2, n_row-2);
        arma::mat vertical_columns_next_len(2, n_row-2);

        for(int i=0; i<(n_row-2); i++) 
        {
            vertical_columns_next_ln1.col(i) = A_ln1( arma::span(i, i+1), i+2 );
            vertical_columns_next_ln2.col(i) = A_ln2( arma::span(i, i+1), i+2 );
            vertical_columns_next_len.col(i) = A_len( arma::span(i, i+1), i+2 );
        }
        
        for(int i=1; i<2; i++) // cumsum of the two vertical cells (i=1:1)
        {
            vertical_columns_next_ln1.row(i) = vertical_columns_next_ln1.row(i) + vertical_columns_next_ln1.row(i-1); //cumsum
            vertical_columns_next_ln2.row(i) = vertical_columns_next_ln2.row(i) + vertical_columns_next_ln2.row(i-1); //cumsum
            vertical_columns_next_len.row(i) = vertical_columns_next_len.row(i) + vertical_columns_next_len.row(i-1); //cumsum
        }
        
        arma::mat rad_mat_next_ln1 = rad_mat_current_ln1;
        arma::mat rad_mat_next_ln2 = rad_mat_current_ln2;
        arma::mat rad_mat_next_len = rad_mat_current_len;

        arma::mat vertical_columns_current_ln1 = vertical_columns_next_ln1; // this line just create the vertical_columns_current_ln
        arma::mat vertical_columns_current_ln2 = vertical_columns_next_ln2; // this line just create the vertical_columns_current_ln
        arma::mat vertical_columns_current_len = vertical_columns_next_len; // this line just create the vertical_columns_current_ln

        // Rcout << L << "\n";
        // cout << "Begin iteration:\n";
        
        // time complexity of this part: n^3
        // for(int shift=3; shift<=n_row; shift++)
        // each row of rad_mat represent one off-diagonal point
        for(int shift=3; shift<=k; shift++)
        {
            rad_mat_current_ln1 = rad_mat_next_ln1;
            rad_mat_current_ln2 = rad_mat_next_ln2;
            rad_mat_current_len = rad_mat_next_len;

            rad_mat_next_ln1 = arma::mat( n_row-shift+1, shift-1);
            rad_mat_next_ln2 = arma::mat( n_row-shift+1, shift-1);
            rad_mat_next_len = arma::mat( n_row-shift+1, shift-1);
            
            n_cells = arma::zeros<arma::vec>(shift-1); // size of each rectangle
            for( int i=0; i< (shift-1); i++ ) n_cells(i) = (i+1)*(shift-i-1); 
            arma::mat rad_mat_next_len_all = (arma::ones<arma::mat>(n_row-shift+1, shift-1))*(arma::diagmat(n_cells)); // In R: rep(vec, n_row times) // Schur product: element-wise multiplication of two objects
            
            for(int i=1; i<=(n_row-shift+1); i++)
            {
                // next = current + vertical_columns_next_ln values
                rad_mat_next_ln1.submat(i-1, 0, i-1, shift-2-1) = rad_mat_current_ln1( i-1, arma::span(0, shift-2-1) ) + vertical_columns_next_ln1( arma::span(0, shift-2-1), i-1 ).t();
                rad_mat_next_ln2.submat(i-1, 0, i-1, shift-2-1) = rad_mat_current_ln2( i-1, arma::span(0, shift-2-1) ) + vertical_columns_next_ln2( arma::span(0, shift-2-1), i-1 ).t();
                rad_mat_next_len.submat(i-1, 0, i-1, shift-2-1) = rad_mat_current_len( i-1, arma::span(0, shift-2-1) ) + vertical_columns_next_len( arma::span(0, shift-2-1), i-1 ).t();

                rad_mat_next_ln1(i-1, shift-1-1) = vertical_columns_next_ln1(shift-1-1, i-1); // the last new element
                rad_mat_next_ln2(i-1, shift-1-1) = vertical_columns_next_ln2(shift-1-1, i-1); // the last new element
                rad_mat_next_len(i-1, shift-1-1) = vertical_columns_next_len(shift-1-1, i-1); // the last new element
            }
            

            
            ////////////////////////////////////// compute the vertical_columns_next values
            if(shift < n_row) //stop when shift=n
            {
                vertical_columns_current_ln1 = vertical_columns_next_ln1;
                vertical_columns_current_ln2 = vertical_columns_next_ln2;
                vertical_columns_current_len = vertical_columns_next_len;

                arma::mat first_row_ln1(1, n_row-shift);
                arma::mat first_row_ln2(1, n_row-shift);
                arma::mat first_row_len(1, n_row-shift);

                for(int i=0; i<(n_row-shift); i++) 
                {
                    first_row_ln1(0, i) = A_ln1(i, i+shift); // off-diagonal values to be appended to vertical_columns_next_ln
                    first_row_ln2(0, i) = A_ln2(i, i+shift); // off-diagonal values to be appended to vertical_columns_next_ln
                    first_row_len(0, i) = A_len(i, i+shift); // off-diagonal values to be appended to vertical_columns_next_ln
                }
                
                vertical_columns_next_ln1 = vertical_columns_current_ln1.submat(0, 1, shift-2, n_row-shift); // drop the first column
                vertical_columns_next_ln2 = vertical_columns_current_ln2.submat(0, 1, shift-2, n_row-shift); // drop the first column
                vertical_columns_next_len = vertical_columns_current_len.submat(0, 1, shift-2, n_row-shift); // drop the first column

                vertical_columns_next_ln1 = arma::join_cols(first_row_ln1, vertical_columns_next_ln1);
                vertical_columns_next_ln2 = arma::join_cols(first_row_ln2, vertical_columns_next_ln2);
                vertical_columns_next_len = arma::join_cols(first_row_len, vertical_columns_next_len);

                for(int i=1; i<shift; i++) 
                {
                    vertical_columns_next_ln1.row(i) = vertical_columns_next_ln1.row(i) + first_row_ln1; // cumsum
                    vertical_columns_next_ln2.row(i) = vertical_columns_next_ln2.row(i) + first_row_ln2; // cumsum
                    vertical_columns_next_len.row(i) = vertical_columns_next_len.row(i) + first_row_len; // cumsum
                }   
            }
            
            ////////////////////////////////////// compute the L and ancestors
            // if(shift >= 4)
            if(shift >= 2*min_n_bins) //            
            {
                ps = rad_mat_next_len_all - rad_mat_next_len; // number of positive values
                loglik_tmp = loglik_lnorm_cpp_mat( rad_mat_next_ln1, rad_mat_next_ln2, ps, rad_mat_next_len );
                // loglik = loglik_tmp.submat(0, 1, n_row-shift, shift-2-1); // remove first and last col because of the min_n_bins. SHOULD BE MODIFIED. submat: X.submat( first_row, first_col, last_row, last_col ), http://arma.sourceforge.net/docs.html#submat
                loglik = loglik_tmp.submat(0, min_n_bins-1, n_row-shift, shift-2-(min_n_bins-1)); // 2018-11-14, remove first and last min_n_bins-1 cols because of the min_n_bins. SHOULD BE MODIFIED. submat: X.submat( first_row, first_col, last_row, last_col ), http://arma.sourceforge.net/docs.html#submat

                
                arma::mat cases(1, shift-2*min_n_bins+1); // shift=5: 1:2, i.e., two cases
                // arma::mat loglik(n_row-shift, shift-min_n_bins);
                for( int row=1; row<=(n_row-shift+1); row++ )
                {
                    p = row; //7
                    q = row + shift; // 7 + 4 = 11
                    cases = loglik( p-1, arma::span(0, shift-2*min_n_bins) ) + L(p-1, arma::span(p-1-1 + min_n_bins, q - min_n_bins-1-1) ) + (L(arma::span(p-1 + min_n_bins, q - min_n_bins-1), q-1-1)).t();
                    L(p-1, q-1-1) = cases.max();

                    max_mid_index = cases.index_max() + 1; //The c++ offset 1
                    max_mid = (p-1 + min_n_bins) + max_mid_index -1; // should minus one
                    // ancestor = paste(i, max_mid, max_mid+1, j, sep='-')
                    ancestors(p-1, q-1-1) = to_string(p) + "-" + to_string(max_mid) + "-" + to_string(max_mid+1) + "-" + to_string(q-1);
                    // cout << "cases:" << L(p-1, q-1-1) << "\n";
                }
            }
            
        }
        
        res["ancestors"] = ancestors;
        res["L"] = L;
        return res;
    }
    
    // compute the ancestors
    // [[Rcpp::export]]
    List compute_L( arma::mat A, arma::mat L, int k ) // A seems not needed here, and can be removed, Y.L, 2018-11-14 (Indeed this part is not used)
    {
        int n_row = A.n_rows;
        int min_n_bins = 2;
        StringMatrix ancestors(n_row, n_row);
    
        // should rewrite this part
        // for( int i= min_n_bins-1; i<=(2*min_n_bins - 2); i++ )
            // for( int j=i; j<n_row; j++ ) L(j-i, j) = 1;
        
        for( int shift=(2*min_n_bins - 1); shift<=(n_row-1); shift++ )  
        {
            int i, j, max_mid_index, max_mid;
            arma::mat cases(1, shift-2*min_n_bins+1);
            arma::mat prob_vec(n_row-shift, shift-min_n_bins);
            for( int row=1; row<=(n_row-shift); row++ )
            {
                    i = row;
                    j = row + shift;
                    cases = prob_vec( i-1, arma::span(0, j-2*min_n_bins-i+1) ) + L(i-1, arma::span(i-1-1 + min_n_bins, j - min_n_bins-1) ) + (L(arma::span(i + min_n_bins-1, j - min_n_bins), j-1)).t();
                    L(i-1, j-1) = cases.max();
                    max_mid_index = cases.index_max() + 1; //The c++ offset 1
                    max_mid = (i-1 + min_n_bins) + max_mid_index -1; // should minus one
                    // ancestor = paste(i, max_mid, max_mid+1, j, sep='-')
                    ancestors(i-1, j-1) = to_string(i) + "-" + to_string(max_mid) + "-" + to_string(max_mid+1) + "-" + to_string(j); // This 'to_string' may be optimized 
            }
        }
    return List::create(Named("L") = L, Named("ancestors") = ancestors);    
    }
    
