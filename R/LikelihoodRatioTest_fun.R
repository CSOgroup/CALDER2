	## The following code tries to evaluate the performance of different likelihood ratio-based tests
	## Yuanlong LIU
	## 01-03-2018
	
	############## model 1: keep the same structure ##############
	## In this case, a tree structure should be provided
	## Constraint: the top three levels have the same contact intensity parameter

	## 	01-03-2018
	p_likelihood_ratio <- function( A, head, mid, tail, num )
	{
		A_whole = A[head:tail, head:tail]
		A_a = A[head:mid, head:mid]
		A_b = A[(mid+1):tail, (mid+1):tail]
		A_ab = A[head:mid, (mid+1):tail]
		
		tri_Awhole = A_whole[upper.tri(A_whole, diag=FALSE)] 
		tri_Aa = A_a[upper.tri(A_a, diag=FALSE)] 
		tri_Ab = A_b[upper.tri(A_b, diag=FALSE)] 
		# gamma_fit(tri_Awhole, num)
		
		inner_a = get_prob( tri_Aa )
		inner_b = get_prob( tri_Ab )
		inter_ab = get_prob( A_ab )
		
		## negative binomial
		# inner_a = get_prob_ng( tri_Aa )
		# inner_b = get_prob_ng( tri_Ab )
		# inter_ab = get_prob_ng( A_ab )

		
		## log likelihood of H1
		LH1 = inner_a + inner_b + inter_ab
		LH0 = get_prob( tri_Awhole )
		
		Lambda = -2*( LH0 - LH1 )
		# cat(Lambda, '\n')
		df_h0 = 1*1
		df_h1 = 3*1
		df = df_h1 - df_h0
		
		p = pchisq(Lambda, df=df, lower.tail = FALSE, log.p = FALSE)
		info = list( Lambda=Lambda, p=p )
		return(info)
	}
	
	p_likelihood_ratio_nb <- function( A, head, mid, tail )
	{
		A_whole = A[head:tail, head:tail]
		A_a = A[head:mid, head:mid]
		A_b = A[(mid+1):tail, (mid+1):tail]
		A_ab = A[head:mid, (mid+1):tail]
		
		tri_Awhole = A_whole[upper.tri(A_whole, diag=TRUE)] 
		tri_Aa = A_a[upper.tri(A_a, diag=TRUE)] 
		tri_Ab = A_b[upper.tri(A_b, diag=TRUE)] 
		
		## negative binomial
		inner_a = get_prob_nb( tri_Aa )
		inner_b = get_prob_nb( tri_Ab )
		inter_ab = get_prob_nb( A_ab )

		
		## log likelihood of H1
		LH1 = inner_a + inner_b + inter_ab
		LH0 = get_prob_nb( tri_Awhole )
		
		Lambda = -2*( LH0 - LH1 )
		# cat(Lambda, '\n')
		n_parameters = 2
		df_h0 = 1*n_parameters
		df_h1 = 3*n_parameters
		df = df_h1 - df_h0
		
		p = pchisq(Lambda, df=df, lower.tail = FALSE, log.p = FALSE)
		info = list( Lambda=Lambda, p=p )
		return(info)
	}

	p_likelihood_ratio_norm <- function( A, head, mid, tail )
	{
		A_whole = A[head:tail, head:tail]
		A_a = A[head:mid, head:mid]
		A_b = A[(mid+1):tail, (mid+1):tail]
		A_ab = A[head:mid, (mid+1):tail]
		
		tri_Awhole = A_whole[upper.tri(A_whole, diag=TRUE)] 
		tri_Aa = A_a[upper.tri(A_a, diag=TRUE)] 
		tri_Ab = A_b[upper.tri(A_b, diag=TRUE)] 
		
		## norm
		inner_a = likelihood_norm( tri_Aa )
		inner_b = likelihood_norm( tri_Ab )
		inter_ab = likelihood_norm( A_ab )

		
		## log likelihood of H1
		LH1 = inner_a + inner_b + inter_ab
		LH0 = likelihood_norm( tri_Awhole )
		
		Lambda = -2*( LH0 - LH1 )
		# cat(Lambda, '\n')
		n_parameters = 2
		df_h0 = 1*n_parameters
		df_h1 = 3*n_parameters
		df = df_h1 - df_h0
		
		p = pchisq(Lambda, df=df, lower.tail = FALSE, log.p = FALSE)
		info = list( Lambda=Lambda, p=p )
		return(info)
	}

	
	
	p_likelihood_ratio_gamma <- function( A, head, mid, tail, n_parameters, imputation )
	{
		A_whole = A[head:tail, head:tail]
		A_a = A[head:mid, head:mid]
		A_b = A[(mid+1):tail, (mid+1):tail]
		A_ab = A[head:mid, (mid+1):tail]
		
		tri_Awhole = A_whole[upper.tri(A_whole, diag=TRUE)] 
		## added 25-03-2018. If no zero values, no imputation
		if( length(tri_Awhole==0) == 0 ) imputation = FALSE
		
		tri_Aa = A_a[upper.tri(A_a, diag=TRUE)] 
		tri_Ab = A_b[upper.tri(A_b, diag=TRUE)] 
		
		## norm
		inner_a = likelihood_gamma_mme( tri_Aa )
		inner_b = likelihood_gamma_mme( tri_Ab )
		inter_ab = likelihood_gamma_mme( A_ab )
		whole = likelihood_gamma_mme( tri_Awhole )
		
		if( imputation ) ## zero values are imputed by the estimated distribution based on non random values
		{
			inner_a = likelihood_gamma_mme( tri_Aa[tri_Aa!=0] )/length( tri_Aa[tri_Aa!=0] )*length( tri_Aa )
			inner_b = likelihood_gamma_mme( tri_Ab[tri_Ab!=0] )/length( tri_Ab[tri_Ab!=0] )*length( tri_Ab )
			inter_ab = likelihood_gamma_mme( A_ab )/length( A_ab[A_ab!=0] )*length( A_ab )
			whole = likelihood_gamma_mme( tri_Awhole[tri_Awhole!=0] )/length( tri_Awhole[tri_Awhole!=0] )*length( tri_Awhole )
			n_parameters = n_parameters - 1 ## the mixture parameter of 0 is not taken into account
		}

		
		## log likelihood of H1
		LH1 = inner_a + inner_b + inter_ab
		LH0 = whole
		
		Lambda = -2*( LH0 - LH1 )
		# cat(Lambda, '\n')
		df_h0 = 1*n_parameters
		df_h1 = 3*n_parameters
		df = df_h1 - df_h0
		
		p = pchisq(Lambda, df=df, lower.tail = FALSE, log.p = FALSE)
		info = list( Lambda=Lambda, p=p )
		return(info)
	}

	lognormal_mean_test <- function( cA, head, mid, tail )
	{
		A_whole = cA[head:tail, head:tail]
		A_a = cA[head:mid, head:mid]
		A_b = cA[(mid+1):tail, (mid+1):tail]
		A_ab = cA[head:mid, (mid+1):tail]
		
		tri_Aa = A_a[upper.tri(A_a, diag=TRUE)] 
		tri_Ab = A_b[upper.tri(A_b, diag=TRUE)] 
		
		tri_Aa_vec = as.vector( tri_Aa )
		tri_Ab_vec = as.vector( tri_Ab )
		A_ab_vec = as.vector( A_ab )
		
		tri_Aa_vec_p = tri_Aa_vec[tri_Aa_vec!=0]
		tri_Ab_vec_p = tri_Ab_vec[tri_Ab_vec!=0]
		A_ab_vec_p = A_ab_vec[A_ab_vec!=0]
		
		p_Aa = p_lognormal_mean( tri_Aa_vec_p, A_ab_vec_p )
		p_Ab = p_lognormal_mean( tri_Ab_vec_p, A_ab_vec_p )
		return( list(p_Aa=p_Aa, p_Ab=p_Ab) )
	}
	
	## https://www.jstor.org/stable/2533570
	## google: Methods for Comparing the Means of Two Independent Log-Normal Samples
	## 03-07-2018
	p_lognormal_mean <- function( vec_aa, vec_ab )
	{
		if( all(vec_aa==0) | all(vec_ab==0) ) return(0)
		
		n_aa = length(vec_aa)
		n_ab = length(vec_ab)
		fited_info_aa = MASS::fitdistr(vec_aa, 'lognormal') ## intra
		mu_aa = fited_info_aa$estimate[1]
		sd_aa = fited_info_aa$estimate[2]
		s2_aa = sd_aa^2

		fited_info_ab = MASS::fitdistr(vec_ab, 'lognormal') ## inter
		mu_ab = fited_info_ab$estimate[1]
		sd_ab = fited_info_ab$estimate[2]
		s2_ab = sd_ab^2
		
		z_score = ( (mu_aa - mu_ab) + (1/2)*(s2_aa - s2_ab) ) / sqrt( s2_aa/n_aa + s2_ab/n_ab + (1/2)*( s2_aa^2/(n_aa-1) + s2_ab^2/(n_ab-1) ) )
		p = pnorm( z_score, lower.tail=FALSE )
		return(p)
	}

	get_corner_xy <- function(A_whole)
	{
		n = nrow(A_whole)
		corner_size = floor(n/2)
		A_corner = A_whole[1:corner_size, (n - corner_size + 1 ):n]
		expected_high = A_corner[upper.tri(A_corner, diag=FALSE)]
		expected_low = A_corner[lower.tri(A_corner, diag=TRUE)]
		return(list(x=expected_low, y=expected_high))
	}

	get_half_mat_values <- function(mat)
	{
		n1 = nrow(mat)
		n2 = ncol(mat)
		delta = n1/n2
		rows = lapply( 1:n2, function(x) (1+ceiling(x*delta)):n1 )
		rows = rows[ sapply(rows, function(v) max(v) <= n1) ]
		flag = which(diff(sapply(rows, length)) > 0)
		if(length(flag)>0) rows = rows[ 1:min(flag) ]
		
		row_col_indices = cbind( unlist(rows), rep(1:length(rows), sapply(rows, length)))
		x = mat[row_col_indices]
		mat[row_col_indices] = NA
		y = na.omit(as.vector(mat))
		return(list(x=x, y=y))
	}


	get_half_mat_values_v2 <- function(mat)
	{
		## https://stackoverflow.com/questions/52990525/get-upper-triangular-matrix-from-nonsymmetric-matrix/52991508#52991508
		y = mat[nrow(mat) * (2 * col(mat) - 1) / (2 * ncol(mat)) - row(mat) > -1/2]
		x = mat[nrow(mat) * (2 * col(mat) - 1) / (2 * ncol(mat)) - row(mat) < -1/2]
		return(list(x=x, y=y))
	}

	p_wilcox_test_nested <- function( A, head, mid, tail, alternative, is_CD )
	{
		test_1 = p_wilcox_test( A, head, mid, tail, alternative, is_CD, only_corner=FALSE ) #: coner + inter
		if( (tail - head <= 4) | (test_1$p > 0.05) ) ## when > 0.05 or it is already small, do not cosider it as nested
		{
			info = list( Lambda=NULL, p=0.555555, mean_diff=0 )
			return(info)
		}

		## try_error happens when tad to test is too small. THEREFORE, ASSIGN P=0 TO THE TAD
		test_left_tad = try(p_wilcox_test( A, head, mid=ceiling((head+mid)/2), mid, alternative, is_CD=FALSE, only_corner=TRUE )) #: coner + inter
		if( class(test_left_tad)=="try-error" ) test_left_tad = list(p=0)
		test_right_tad = try(p_wilcox_test( A, mid+1, mid=ceiling((mid+1+tail)/2), tail, alternative, is_CD=FALSE, only_corner=TRUE )) #: coner + inter
		if( class(test_right_tad)=="try-error" ) test_right_tad = list(p=0)

		info = list( Lambda=NULL, p=max( test_1$p, min(test_left_tad$p, test_right_tad$p) ), mean_diff=0 )
		return(info)
	}



	p_wilcox_test = function( A, head, mid, tail, alternative, is_CD=FALSE, only_corner=FALSE ) ## only_corner tests if a domain is a TAD (no nesting)
        {
                A_whole = A[head:tail, head:tail]
                A_a = A[head:mid, head:mid]
                A_b = A[(mid+1):tail, (mid+1):tail]
                A_ab = A[head:mid, (mid+1):tail]

                tri_Awhole = A_whole[upper.tri(A_whole, diag=TRUE)] ## need to check whether diag should be false or true

                tri_Aa = A_a[upper.tri(A_a, diag=TRUE)] ## need to check whether diag should be false or true
                tri_Ab = A_b[upper.tri(A_b, diag=TRUE)] ## need to check whether diag should be false or true

                tri_Aa_vec = as.vector( tri_Aa )
                tri_Ab_vec = as.vector( tri_Ab )
                A_ab_vec = as.vector( A_ab )

                corner_mat_info = get_half_mat_values_v2(A_ab)
                A_ab_corner = as.vector(corner_mat_info$y)
                A_ab_ncorner = corner_mat_info$x
                p_corner = wilcox.test(x=A_ab_corner, y=A_ab_ncorner, alternative='greater', exact=F)
                p_inter = wilcox.test(x=c(tri_Ab_vec, tri_Aa_vec), y=A_ab_ncorner, alternative='greater', exact=F)
                # p_inter = wilcox.test(x=c(tri_Ab_vec, tri_Aa_vec), y=c(A_ab_ncorner, A_ab_corner), alternative='greater', exact=F)

                if(is_CD==FALSE) p = max(p_inter$p.value, p_corner$p.value) ## if the tested part is the CD
                if(is_CD==TRUE) p = p_inter$p.value ## if the tested part is the CD
                # if(only_corner==TRUE) p = p_corner$p.value
                mean_diff_inter = mean(A_ab_ncorner) - mean(c(tri_Ab_vec, tri_Aa_vec)) ## negative is good
                mean_diff_corner = mean(A_ab_ncorner) - mean(A_ab_corner) ## negative is good


                if(is_CD==FALSE) mean_diff = min(mean_diff_corner, mean_diff_inter) ## if the tested part is the CD
                if(is_CD==TRUE) mean_diff = mean_diff_inter

                if(min(length(tri_Ab_vec), length(tri_Ab_vec)) < 10) mean_diff = 100 ## when one of the two twins is too small. 10: dim(4,4)

                # mean_diff = mean(A_ab_corner) - mean(A_ab_ncorner)

                # p_test_Aa =  wilcox.test(x=A_ab_vec, y=tri_Aa_vec, alternative="less", exact=F)
                # p_test_Ab =  wilcox.test(x=A_ab_vec, y=tri_Ab_vec, alternative="less", exact=F)
                # p = wilcox.test(x=c(tri_Ab_vec, tri_Aa_vec), y=A_ab_vec, alternative=alternative, exact=F)

                # xy = get_corner_xy(A_whole)
                # p = wilcox.test(x=xy$x, y=xy$y, alternative=alternative, exact=F)

                # p = max(p_test_Aab$p.value, p_test_Ab$p.value)
                info = list( Lambda=NULL, p=p, mean_diff=mean_diff)
                return(info)
        }



	
	p_likelihood_ratio_lnorm <- function( A, head, mid, tail, n_parameters, imputation, imputation_num=1E2 )
	{
		likelihood_fun = likelihood_lnorm_mle
		A_whole = A[head:tail, head:tail]
		A_a = A[head:mid, head:mid]
		A_b = A[(mid+1):tail, (mid+1):tail]
		A_ab = A[head:mid, (mid+1):tail]
		
		tri_Awhole = A_whole[upper.tri(A_whole, diag=TRUE)] 
		## added 25-03-2018. If no zero values, no imputation
		no_zero_flag = 0
		if( sum(tri_Awhole==0) == 0 ) { no_zero_flag = 1; imputation = FALSE }
		
		tri_Aa = A_a[upper.tri(A_a, diag=TRUE)] 
		tri_Ab = A_b[upper.tri(A_b, diag=TRUE)] 
		
		tri_Aa_vec = as.vector( tri_Aa )
		tri_Ab_vec = as.vector( tri_Ab )
		A_ab_vec = as.vector( A_ab )

		mean_a = mean( tri_Aa_vec[tri_Aa_vec!=0] )
		mean_b = mean( tri_Ab_vec[tri_Ab_vec!=0] )
		mean_ab = mean( A_ab_vec[A_ab_vec!=0] )
		mean_diff = mean_ab - min(mean_a, mean_b)
		
		if( (all(tri_Aa_vec==0)) | (all(tri_Ab_vec==0)) | (all(A_ab_vec==0))  )
		{
			info = list( Lambda=NA, p=0, mean_diff=mean_diff )
			return(info)
		}

		## lnorm
		if(!imputation)
		{
			inner_a = likelihood_fun( tri_Aa )
			inner_b = likelihood_fun( tri_Ab )
			inter_ab = likelihood_fun( A_ab )
			whole = likelihood_fun( tri_Awhole )

			## log likelihood of H1
			LH1 = inner_a + inner_b + inter_ab
			LH0 = whole
			Lambda = -2*( LH0 - LH1 )
			
			if( no_zero_flag ) n_parameters = n_parameters - 1 ## no alpha parameter
			df_h0 = 1*n_parameters ##(theta_1 = theta_2 = theta_3)
			df_h1 = 3*n_parameters ##(theta_1, theta_2, theta_3)
			df = df_h1 - df_h0
		
			p = pchisq(Lambda, df=df, lower.tail = FALSE, log.p = FALSE)
			info = list( Lambda=Lambda, p=p, mean_diff=mean_diff )
			return(info)
		}
		
		if( imputation ) ## zero values are imputed by the estimated distribution based on non random values
		{
			vec_list = list( tri_Aa=tri_Aa, tri_Ab=tri_Ab, A_ab=A_ab )
			imputations = imputation_list( vec_list, imputation_num )
			
			inner_as = apply(imputations$tri_Aa, 1, likelihood_fun)
			inner_bs = apply(imputations$tri_Ab, 1, likelihood_fun)
			inter_abs = apply(imputations$A_ab, 1, likelihood_fun)
			
			wholes = apply(do.call( cbind, imputations ), 1, likelihood_fun)
			LH1s = inner_as + inner_bs + inter_abs
			LH0s = wholes
			Lambdas = -2*( LH0s - LH1s )
			Lambda = mean( Lambdas )
			n_parameters = n_parameters - 1 ## the mixture parameter is not taken into account
			
			# cat(Lambda, '\n')
			df_h0 = 1*n_parameters
			df_h1 = 3*n_parameters
			df = df_h1 - df_h0
			
			p = pchisq(Lambda, df=df, lower.tail = FALSE, log.p = FALSE)
			info = list( Lambda=Lambda, p=p, mean_diff=mean_diff )
			return(info)
			
		# if( imputation ) ## zero values are imputed by the estimated distribution based on non random values
		# {
			# inner_a = likelihood_lnorm( tri_Aa[tri_Aa!=0] )/length( tri_Aa[tri_Aa!=0] )*length( tri_Aa )
			# inner_b = likelihood_lnorm( tri_Ab[tri_Ab!=0] )/length( tri_Ab[tri_Ab!=0] )*length( tri_Ab )
			# inter_ab = likelihood_lnorm( A_ab )/length( A_ab[A_ab!=0] )*length( A_ab )
			# whole = likelihood_lnorm( tri_Awhole[tri_Awhole!=0] )/length( tri_Awhole[tri_Awhole!=0] )*length( tri_Awhole )
			# n_parameters = n_parameters - 1 ## the mixture parameter of 0 is not taken into account
		# }
		
		}
	}

	imputation_list <- function( vec_list, imputation_num )
	{
		imputations = lapply( vec_list, function(v) 
		{
			if(sum(v!=0)==1) {final_vec = matrix( v[v!=0], imputation_num, length(v) ); return(final_vec)}
			
			fit = fit_lnorm(v)
			set.seed(1)
			
			## THERE WILL BE ERROR IF IS.NA SDLOG
			if(!is.na(fit['sdlog'])) imputation_vec = matrix(rlnorm( sum(v==0)*imputation_num, fit['meanlog'], fit['sdlog'] ), nrow=imputation_num)
			if(is.na(fit['sdlog'])) stop("In function imputation_list, sdlog=NA encountered")

			ori_vec =  t(replicate(imputation_num, v[v!=0]))
			final_vec = cbind( ori_vec, imputation_vec )
		} )
		return( imputations )
	}


	fit_lnorm <- function(vec)
	{
		vec = vec[vec!=0]
		fit = MASS::fitdistr(vec, 'lognormal')$estimate
		return(fit)
	}

	