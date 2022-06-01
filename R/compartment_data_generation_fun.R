	
	compress_mat_fast_tmp = function(input_mat, compress_size)
	{
	    mat_block <- function(n, r) suppressWarnings( matrix(c(rep(1, r), rep(0, n)), n, n/r) )
	    n = ncol(input_mat)
	    mat2prod = mat_block(n, compress_size)
	    # return(t(input_mat%*%mat2prod / compress_size))
	    return(t(matrix_multiplication_cpp(input_mat, mat2prod)))
	}

	mat_10to40kb = function(mat2compress, bin_size2look, bin_size_input) ## small bin size to bigger bin_size
	{
		compress_size = bin_size2look / bin_size_input
		len = nrow(mat2compress) - nrow(mat2compress)%%compress_size
		mat2compress = mat2compress[1:len, 1:len]
		c_mat_compressed = compress_mat_fast_tmp( as.matrix(mat2compress), compress_size=compress_size )
		mat_compressed = compress_mat_fast_tmp( c_mat_compressed, compress_size=compress_size )
		rownames(mat_compressed) = colnames(mat_compressed) = as.character( 1:nrow(mat_compressed) )
		return(mat_compressed)
	}

	contact_mat_processing_v2 = function(contact_tab_dump=NULL, contact_file_dump=NULL, contact_file_hic=NULL, chr_num, bin_size_input, bin_size2look, black_list_bins=NULL)
	{	
		compress_size = ifelse(bin_size2look < 40E3, 1, 1)
		zero_ratio = 0.01


		if(!is.null(contact_tab_dump)) contact_mat_raw = contact_tab_dump ## if contact matrix in data.frame or data.table of strawr format is provided
		if(!is.null(contact_file_dump)) contact_mat_raw = data.table::fread(contact_file_dump) ## if contact matrix in strawr format is provided

		if(!is.null(contact_file_hic)) ## if contact matrix in hic format is provided
		{
			bin_sizes_in_hic = strawr::readHicBpResolutions(contact_file_hic)
			chrs_in_hic = as.vector(strawr::readHicChroms(contact_file_hic)[[1]])
			chr2query = chrs_in_hic[match(toupper(chr_num), gsub('chr', '', toupper(chrs_in_hic), ignore.case=TRUE))]

			if(!(bin_size_input %in% bin_sizes_in_hic)) stop(sprintf('Your provided hic file only contains resolutions of: %s', paste0(bin_sizes_in_hic, collapse=' ')))
			if(!(chr2query %in% chrs_in_hic)) stop(sprintf('Your provided hic file only contains chrs of: %s', paste0(chrs_in_hic, collapse=' ')))

			## try different normalization to get available dataset

			contact_mat_raw = try(strawr::straw("KR", contact_file_hic, as.character(chr2query), as.character(chr2query), "BP", bin_size_input))
			if(class(contact_mat_raw)=='try-error' | (class(contact_mat_raw)!='try-error' & nrow(na.omit(contact_mat_raw)) < 100)) contact_mat_raw = try(strawr::dump("VC_SQRT", contact_file_hic, as.character(chr2query), as.character(chr2query), "BP", bin_size_input))
			if(class(contact_mat_raw)=='try-error' | (class(contact_mat_raw)!='try-error' & nrow(na.omit(contact_mat_raw)) < 100)) contact_mat_raw = try(strawr::dump("VC", contact_file_hic, as.character(chr2query), as.character(chr2query), "BP", bin_size_input))
			if(class(contact_mat_raw)=='try-error' | (class(contact_mat_raw)!='try-error' & nrow(na.omit(contact_mat_raw)) < 100)) stop(sprintf('Your provided hic file does not contain information given the bin_size=%s and any of the normalization method KR/VC/VC_SQRT', bin_size_input))	
			contact_mat_raw = data.table::as.data.table(contact_mat_raw)
		}

		colnames(contact_mat_raw) = c('pos_1', 'pos_2', 'val')	

		contact_mat = subset(contact_mat_raw, !is.na(val))
		contact_mat[,1] = contact_mat[,1]/bin_size_input
		contact_mat[,2] = contact_mat[,2]/bin_size_input

		if(!all(contact_mat[[2]] >= contact_mat[[1]])) stop('\nYour provided matrix does not represent an upper triangular matrix!\n\n')

		n_bins = max(max(contact_mat[[1]]), max(contact_mat[[2]])) + 1 ## should +1 because contact_mat index starts from 0 (bin 0 represents: 0-10E3, checked by looking at the juicebox map, 2018-11-19)
		mat_sparse = Matrix::Matrix(0, nrow=n_bins, ncol=n_bins)
		mat_sparse[cbind(contact_mat[[1]]+1, contact_mat[[2]]+1)] <- contact_mat[[3]]

		rownames(mat_sparse) = colnames(mat_sparse) = as.character( 1:nrow(mat_sparse) )
		
		########################################################## remove black listed regions

		if(length(black_list_bins) > 0)
		{
			black_list_bins = intersect(as.character(black_list_bins), rownames(mat_sparse))
			if(length(black_list_bins) > 0) mat_sparse[black_list_bins, ] = mat_sparse[, black_list_bins] = 0
		}

		########################################################## remove bins having too sparse contacts

		mat_sparse = Matrix::forceSymmetric(mat_sparse, uplo='U')
		if(bin_size2look!=bin_size_input) mat_sparse = mat_10to40kb( mat_sparse, bin_size2look, bin_size_input )
		mat_dense = remove_blank_cols(mat_sparse, sparse=TRUE, ratio=zero_ratio) ## has the same rows/cols as A
		if(nrow(mat_dense) < 100) mat_dense = remove_blank_cols(mat_sparse, sparse=TRUE, ratio=0) ## when all are dense
		while(min(apply(mat_dense, 1, sd))==0) ## sometimes after removing the cols / rows, the remained part will all be 0
		{
			mat_dense = remove_blank_cols(mat_dense, sparse=TRUE, ratio=1E-7) ## has the same rows/cols as A
			if(nrow(mat_dense) < 1) stop('Error in generating mat_dense at the data generating step')
		}

		##########################################################

		nrow2keep = nrow(mat_dense) - nrow(mat_dense)%%compress_size
		mat_dense_2_compress = mat_dense[, 1:nrow2keep]

		bin_names = rownames(mat_dense)

		mat_dense_compressed = compress_mat_fast( as.matrix(mat_dense_2_compress), compress_size=compress_size )
		colnames(mat_dense_compressed) = bin_names
		rm(mat_dense_2_compress); gc()
		
		# range(mat_dense_compressed)
		# # sum(mat_dense_compressed > 1000)
		# # mat_dense_compressed[mat_dense_compressed > 1000] = 1000
		# mat_dense_compressed_sparse = mat_dense_compressed
		mat_dense_compressed = as.matrix(mat_dense_compressed)
		mat_dense_compressed_log = log2(mat_dense_compressed + 1)
		
		# #########################################################
		# cat('compute correlation matrix ... ')

		cmat_dense_compressed_log = fast_cor(mat_dense_compressed_log)
		ccmat_dense_compressed_log = fast_cor(cmat_dense_compressed_log)

		# cat('compute correlation matrix done ... ')

		# #########################################################
		# # ccmat_dense_compressed_atanh = atanh(ccmat_dense_compressed - 1E-7)
		ccmat_dense_compressed_log_atanh = atanh(ccmat_dense_compressed_log / (1+1E-7))

		# rm(mat_dense_compressed, mat_dense_compressed_sparse, cmat_dense_compressed_log)
		gc()
		# #########################################################
		# cat('ready to compute compartment domains\n')

		out = list(mat_dense=mat_dense, atanh_score=ccmat_dense_compressed_log_atanh)

		return(out)
	}
	


	# if(!file.exists(contact_mat_file)) contact_mat_file = paste0('/mnt/ndata/Yuanlong/2.Results/1.Juicer/', CELL_LINE, '/contact_mat/mat_chr', chr, '_', bin_size_initial_kb, 'kb_ob.txt.gz')

	## bin_size_initial is the binsize of your input matrix, can be different from the bin_size of your planned analysis
	# contact_mat_processing = function(contact_mat_file, bin_size, bin_size_initial=bin_size)
	# {	
		
	# 	compress_size = ifelse(bin_size < 40E3, 1, 1)
	# 	zero_ratio = 0.01

	# 	combined_xk_oe_raw = data.table::fread(contact_mat_file)

	# 	## this code generates the compartment domains

	# 	combined_xk_oe_raw = subset(combined_xk_oe_raw, !is.na(V3))
	# 	combined_xk_oe_raw[,1] = combined_xk_oe_raw[,1]/bin_size_initial
	# 	combined_xk_oe_raw[,2] = combined_xk_oe_raw[,2]/bin_size_initial
	# 	combined_xk_oe = combined_xk_oe_raw

	# 	colnames(combined_xk_oe) = c('pos_1', 'pos_2', 'val')	
	# 	if(!all(combined_xk_oe[[2]] >= combined_xk_oe[[1]])) stop('\nYou provided matrix does not represent an upper triangular matrix!\n\n')

	# 	oe_size = max(max(combined_xk_oe[[1]]), max(combined_xk_oe[[2]])) + 1 ## should +1 because combined_xk_oe index starts from 0 (bin 0 represents: 0-10E3, checked by looking at the juicebox map, 2018-11-19)
	# 	mat_oe_sparse = Matrix::Matrix(0, nrow=oe_size, ncol=oe_size)
	# 	mat_oe_sparse[cbind(combined_xk_oe[[1]]+1, combined_xk_oe[[2]]+1)] <- combined_xk_oe[[3]]

	# 	rownames(mat_oe_sparse) = colnames(mat_oe_sparse) = as.character( 1:nrow(mat_oe_sparse) )
		
	# 	mat_oe_sparse = Matrix::forceSymmetric(mat_oe_sparse, uplo='U')
	# 	if(bin_size!=bin_size_initial) mat_oe_sparse = mat_10to40kb( mat_oe_sparse, bin_size, bin_size_initial )
	# 	A_oe = remove_blank_cols(mat_oe_sparse, sparse=TRUE, ratio=zero_ratio) ## has the same rows/cols as A
	# 	if(nrow(A_oe) < 100) A_oe = remove_blank_cols(mat_oe_sparse, sparse=TRUE, ratio=0) ## when all are dense
	# 	while(min(apply(A_oe, 1, sd))==0) ## sometimes after removing the cols / rows, the remained part will all be 0
	# 	{
	# 		A_oe = remove_blank_cols(A_oe, sparse=TRUE, ratio=1E-7) ## has the same rows/cols as A
	# 		if(nrow(A_oe) < 1) stop('ERROR IN GENERATING MEANINGFUL A_oe at the data generating step')
	# 	}

	# 	##########################################################

	# 	len = nrow(A_oe) - nrow(A_oe)%%compress_size
	# 	A_oe_2_compress = A_oe[, 1:len]

	# 	bin_names = rownames(A_oe)

	# 	A_oe_compressed = compress_mat_fast( as.matrix(A_oe_2_compress), compress_size=compress_size )
	# 	colnames(A_oe_compressed) = bin_names
	# 	rm(A_oe_2_compress); gc()
		
	# 	range(A_oe_compressed)
	# 	# # sum(A_oe_compressed > 1000)
	# 	# # A_oe_compressed[A_oe_compressed > 1000] = 1000
	# 	A_oe_compressed_sparse = A_oe_compressed
	# 	A_oe_compressed = as.matrix(A_oe_compressed)
	# 	A_oe_compressed_log = log2(A_oe_compressed + 1)
		
	# 	# #########################################################
	# 	# cat('compute correlation matrix ... ')

	# 	cA_oe_compressed_log = fast_cor(A_oe_compressed_log)
	# 	ccA_oe_compressed_log = fast_cor(cA_oe_compressed_log)

	# 	# cat('compute correlation matrix done ... ')

	# 	# #########################################################
	# 	# # ccA_oe_compressed_atanh = atanh(ccA_oe_compressed - 1E-7)
	# 	ccA_oe_compressed_log_atanh = atanh(ccA_oe_compressed_log / (1+1E-7))

	# 	# rm(A_oe_compressed, A_oe_compressed_sparse, cA_oe_compressed_log)
	# 	gc()
	# 	# #########################################################
	# 	# cat('ready to compute compartment domains\n')

	# 	out = list(A_oe=A_oe, atanh_score=ccA_oe_compressed_log_atanh)

	# 	return(out)
	# }

