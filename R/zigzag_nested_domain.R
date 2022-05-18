	## Yuanlong LIU
	## 12/04/2018
	## 06/06/2018
	## 16/06/2018
	## This code runs zigzag search on each compartment domain as a whole without dividing into 400

	## 09-10-2018
	
	## A should be already named
	## "zigzag" resembles computational path of finding the optimum  
	# HRG_zigzag_compartment_domain_main_fun <- function(A, res_dir, compartment_segs, allowed_max_nbins_seq, max_nbins_fine, chr, min_n_bins=2)
	HRG_zigzag_compartment_domain_main_fun <- function(A, res_dir, compartment_segs, allowed_max_nbins_seq, max_nbins_fine, min_n_bins=2)
	{
		arg_list = as.list(environment())
		
		res_folder = file.path(res_dir)
		dir.create(res_folder, recursive=TRUE, showWarnings = FALSE)
		
		total_execution_time_file = file.path(res_dir, 'total_execution.time')
		time_begin = Sys.time()
		cat('Execution begins:', as.character(time_begin), '\n', file=total_execution_time_file, append=TRUE)
	
		## check whether your input matrix is symmetric
		# A_sym = as.matrix( Matrix::forceSymmetric(data.matrix(A), uplo='U') )
		# A_sym = Matrix::forceSymmetric(data.matrix(A), uplo='U') ## keep sparse, 2018-11-11
		A_sym = Matrix::forceSymmetric(A, uplo='U') ## keep sparse, 2018-11-11

		tol = 100 * .Machine$double.eps
		max_diff = max(abs(A_sym - A))
		notSym_flag = max_diff > tol
		if( notSym_flag ) warning('Your input contact matrix is not symmetric. The maximum difference between symmetric values is: ', max_diff, '\nBy default, the contact profile used for downstream analysis is taken from the upper triangular part. To use the lower triangular part you can transpose the input contact matrix first')
		
		if( is.null(rownames(A)) | is.null(colnames(A))) stop('A should be named by the bin indices')
		
		# pA_sym = rm_zeros(A_sym) ## pA_sym: positive A
		pA_sym = remove_blank_cols(A_sym, sparse=TRUE, ratio=0)		
		n_zero_rows = nrow(A_sym) - nrow(pA_sym)
		zero_rows_flag = n_zero_rows > 0
		if( zero_rows_flag )
		{
			warning('There are ', n_zero_rows, ' rows/columns in your input contact matrix that have all their values being 0. These rows/columns are removed for downstream analysis')
			original_row_names = rownames(A_sym)
			kept_row_names = rownames(pA_sym)
		}

		# res_inner = rep(list(), nrow(compartment_segs)) ## for each compartment domain
		# for( i in 1:nrow(compartment_segs) )
		# {
		# 	seg = compartment_segs[i,1]:compartment_segs[i,2]
		# 	cat('Compute seg:', i, 'of length:', length(seg), '\n')

		# 	A_seg = as.matrix(pA_sym[seg, seg])
		# 	res_zigzag = zigzag_loglik_ancestors_v4(A_seg, nrow(A_seg))
		# 	res_outer = list(A=A_seg, L=res_zigzag$L, ancestors=res_zigzag$ancestors)
		# 	res_inner[[i]] = res_outer
		# 	cat('finished', '\n')
		# }

		## changed to paralell, 2018-11-11
		cat('\n')

		res_inner = foreach::foreach(i=1:nrow(compartment_segs)) %do%
		{
			seg = compartment_segs[i,1]:compartment_segs[i,2]
			cat('\r', sprintf('Find sub-domains in %d of %d CDs | length of current CD: %d bins', i, nrow(compartment_segs), length(seg)))

			A_seg = as.matrix(pA_sym[seg, seg])
			res_zigzag = zigzag_loglik_ancestors_v4_5(A_seg, nrow(A_seg), min_n_bins=min_n_bins)
			res_outer = list(A=A_seg, L=res_zigzag$L, ancestors=res_zigzag$ancestors, min_n_bins=min_n_bins)
			res_outer
			# res_inner[[i]] = res_outer
		}

		cat('\n')
		
		segmentss = compartment_segs
		res_info = list( arg_list=arg_list, pA_sym=pA_sym, A_final=pA_sym, segmentss=segmentss, res_inner=res_inner )
		# res_folder_final = file.path(res_dir, 'final')
		# dir.create(res_folder_final, recursive=TRUE, showWarnings = TRUE)
		# save(res_info, file=file.path(res_folder_final, 'res_info.Rdata'))
	
		time_finish = Sys.time()
		cat('Execution finishes:', as.character(time_finish), '\n\n', file=total_execution_time_file, append=TRUE)
		cat('Total execution time:', capture.output( time_finish - time_begin ), '\n\n', file=total_execution_time_file, append=TRUE)
	
		return( res_info )
	}
	
	
	
	
	
	
	
	