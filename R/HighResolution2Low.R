	## These functions try to aggregate bin-wise contact matrix into a domain-wise contact matrix

	HighResolution2Low_k_complete <-function( A, max_nbins )
	{
		nbins = nrow( A )
		
		max_nbins_bp = max_nbins
		max_nbinss = max_nbins + (-5:20)
		expand_fold_info = t(sapply(max_nbinss, fold_expand_rate, K=nbins))
		index = max(which(expand_fold_info[,'rescale_rate']==min(expand_fold_info[,'rescale_rate'])))
		max_nbins = max_nbinss[ index ]
		n2one_head = expand_fold_info[index, 'n2one_head']
		n2one_last = expand_fold_info[index, 'n2one_last']
		
		if(n2one_head==1) return(list())
		rescale_rate = n2one_last / n2one_head

		if(max_nbins!=max_nbins_bp) warning('Value of max_nbins has been adjusted from ', max_nbins_bp, ' to ', max_nbins, '. This results in rescale_rate as: ', rescale_rate, '\n')
		
		if(rescale_rate!=1) ## only do when rescale is not 1
		{
			aa_end = n2one_head*(max_nbins-1)
			bb_end = nrow(A)
			A_aa = A[1:aa_end, 1:aa_end]
			A_bb = A[(aa_end+1):bb_end, (aa_end+1):bb_end]
			A_ab = A[1:aa_end, (aa_end+1):bb_end, drop=FALSE]
			A_aa_low = HighResolution2Low_k(mat=A_aa, k=max_nbins-1)
			A_bb_low = sum( A_bb )
			
			## https://stackoverflow.com/questions/15265512/summing-every-n-points-in-r
			v = apply(A_ab, 1, sum)
			A_ab_low = unname(tapply(v, (seq_along(v)-1) %/% n2one_head, sum))
			
			A_ab_low_new = A_ab_low / n2one_last*n2one_head
			A_bb_low_new = A_bb_low / n2one_last^2*n2one_head^2
			A_low = rbind(cbind(A_aa_low, A_ab_low_new), c( A_ab_low_new, A_bb_low_new ))
			colnames(A_low) = NULL
		}
		
		if(rescale_rate==1) A_low = HighResolution2Low_k(mat=A, k=max_nbins) ## k is the dimension of A_final_low
		## A_low has no row/col name
		res = list(A_low=A_low, n2one_head=n2one_head, n2one_last=n2one_last)
		return(res)
	}	
	
	split_info <-function( K, max_nbins )
	{
		x = floor( K/max_nbins ) + 1
		n2 = max_nbins*x - K
		n1 = max_nbins - n2
		vec = c(rep(x, n1), rep(x-1, n2))
		split_vec = split(1:K, rep(1:max_nbins, vec))

		expand_fold_ratio = min(n1/n2, n2/n1)
		res = list(vec=vec, expand_fold_ratio=expand_fold_ratio, split_vec=split_vec)
		return(res)
	}
	
	
	## diagonal values are summed twice
	HighResolution2Low <-function( A, rescale, which=2 )
	{
		nbins = nrow(A)
		nbins_low = floor(nbins/rescale)
		A_low = matrix( , nbins_low, nbins_low)
		keep_rows = nbins - nbins%%nbins_low
		A_truncated = A[ 1:keep_rows, 1:keep_rows ]
		
		if(which==1) A_low = matsplitter_sum(A_truncated, keep_rows/nbins_low, keep_rows/nbins_low)
		if(which==2) A_low = mat_split_sum(A_truncated, keep_rows/nbins_low, keep_rows/nbins_low)

		return( A_low )
	}
	
	## https://stackoverflow.com/questions/24299171/function-to-split-a-matrix-into-sub-matrices-in-r
	## Function to split a matrix into sub-matrices in R
	matsplitter_sum <- function(M, r, c) 
	{
		rg <- (row(M)-1)%/%r + 1
		cg <- (col(M)-1)%/%c + 1
		rci <- (rg-1)*max(cg) + cg
		N <- prod(dim(M))/r/c
		cv <- unlist(lapply(1:N, function(x) M[rci==x]))
		dim(cv)<-c(r, c, N)
		res = matrix(apply(cv, 3, sum), nrow(M)/r, byrow=TRUE)
		return(res)
	}


	mat_split_sum <- function(M, r, c)
	{
	  nr <- ceiling(nrow(M)/r)
	  nc <- ceiling(ncol(M)/c)
	  newM <- matrix(NA, nr*r, nc*c)
	  newM[1:nrow(M), 1:ncol(M)] <- M

	  div_k <- kronecker(matrix(seq_len(nr*nc), nr, byrow = TRUE), matrix(1, r, c))
	  matlist <- split(newM, div_k)
	  res = matrix(sapply(matlist, sum), nrow(M)/r, byrow=TRUE)
	  return(res)
	}
	
	## this is much faster
	HighResolution2Low_k <- function(mat, k)
	{
		chunk2 <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE)) 
		n = nrow(mat)
		A_low = matrix( , k, k )
		indices = chunk2( 1:n, k )
		for(row_index in 1:k)
		{
			A_sub = mat[indices[[row_index]], ]
			for(col_index in 1:k) A_low[row_index, col_index] = sum( A_sub[ , indices[[col_index]]] )
		}
		return( A_low )
	}
	
	HighResolution2Low_k_rectangle <- function(mat, row_split, col_split, sum_or_mean=c('sum', 'mean', 'median', 'p_above_one'))
	{
		# if(remove_zero==TRUE) mat[mat==0] = NA
		k_row = length(row_split)
		k_col = length(col_split)
		
		A_low = matrix( , k_row, k_col )
		
		if(sum_or_mean=='sum')
		{
			for(row_index in 1:k_row)
			{
				A_sub = mat[row_split[[row_index]], , drop=FALSE ]
				for(col_index in 1:k_col) A_low[row_index, col_index] = sum( A_sub[ , col_split[[col_index]]] )
			}
		}

		if(sum_or_mean=='mean')
		{
			for(row_index in 1:k_row)
			{
				A_sub = mat[row_split[[row_index]], , drop=FALSE ]
				for(col_index in 1:k_col) A_low[row_index, col_index] = mean( A_sub[ , col_split[[col_index]]], na.rm=TRUE )
			}
		}

		if(sum_or_mean=='median')
		{
			for(row_index in 1:k_row)
			{
				A_sub = mat[row_split[[row_index]], , drop=FALSE ]
				for(col_index in 1:k_col) A_low[row_index, col_index] = median( A_sub[ , col_split[[col_index]]], na.rm=TRUE )
			}
		}

		if(sum_or_mean=='p_above_one')
		{
			for(row_index in 1:k_row)
			{
				A_sub = mat[row_split[[row_index]], , drop=FALSE ]
				for(col_index in 1:k_col) A_low[row_index, col_index] = mean( A_sub[ , col_split[[col_index]]] > 1 )
			}
		}
		# if(remove_zero==TRUE) A_low[is.na(A_low)] = 0
		return( A_low )
	}	
	
	## for a n x n matrix, this function generate a nxm matrix, with m*compress_size = n
	## this function is used for generating the nxm matrix, where the i,jth value is the contact value 
	## 08-08-2018  
	compress_mat_fast = function(input_mat, compress_size)
	{
		mat_block <- function(n, r) suppressWarnings( matrix(c(rep(1, r), rep(0, n)), n, n/r) )
		n = ncol(input_mat)
		mat2prod = mat_block(n, compress_size)
		# return(t(input_mat%*%mat2prod / compress_size))
		return(t(matrix_multiplication_cpp(input_mat, mat2prod)) / compress_size)
	}
	
	