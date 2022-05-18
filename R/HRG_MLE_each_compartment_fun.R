	## Yuanlong LIU
	## 12/04/2018
	## 23/05/2018
	## 29/05/2018
	
	## A should be already named
	## A should be symmetric
	## Does not remove 0-rows/columns of A
	HRG_MLE_each_compartment_fun <- function(pA_sym, seg, allowed_max_nbins_seq=500:1000, max_nbins_fine=max(allowed_max_nbins_seq), ncores, res_dir, fast, adjust_segmentss_topdom=TRUE)
	{
		min_n_bins_outer=2
		min_n_bins_inner=2
		distr = 'lnorm'
		A = pA_sym[seg, seg]
		if( is.null(rownames(A)) | is.null(colnames(A))) stop('A should be named by the bin indices')
		arg_list = as.list(environment())
		res_folder = file.path(res_dir)
		dir.create(res_folder, recursive=TRUE, showWarnings = TRUE)
		
		total_execution_time_file = file.path(res_dir, 'total_execution.time')
		cat('n_core used:', ncores, '\n\n', file=total_execution_time_file, append=FALSE)
		time_begin = Sys.time()
		cat('Execution begins:', as.character(time_begin), '\n', file=total_execution_time_file, append=TRUE)
	
		## max_nbins is defined as the one from allowed_max_nbins_seq has the smallest residue, from which, the smallest is taken
		n_A = nrow(A)
		
		## IF THE SEGMENT IS ALREADY SMALL ENOUGH
		if( n_A <= max(allowed_max_nbins_seq) ) 
		{
			## This will be modified
			{
				dists = (2*min_n_bins_outer - 1):(nrow(A)-1)
				counts = sapply(dists, function(v) n_cells2compute( A, v, min_n_bins_outer ))
				split_info = split2BanlancedChunks_min_max(vec=counts, K=min(length(counts)-1, ncores))
				chunks = c(dists[1]-1, dists[split_info$break_points])
				res_outer = HRG_core( A=A, chunks=chunks, res_folder=NULL, min_n_bins=min_n_bins_outer, distr=distr, fast=fast )
			}
			
			## this part gets the hi_tree
			{
				hi_tree = get_tree_v0( res_outer )
				igraph::V(hi_tree)$left_rescaled = igraph::V(hi_tree)$left
				igraph::V(hi_tree)$right_rescaled = igraph::V(hi_tree)$right
				igraph::V(hi_tree)$width_rescaled = igraph::V(hi_tree)$right - igraph::V(hi_tree)$left + 1
				igraph::V(hi_tree)$name = paste('(',igraph::V(hi_tree)$left_rescaled, ',', igraph::V(hi_tree)$right_rescaled, ')', sep='')
			}
			
			full_tree = hi_tree
			res_info = list( arg_list=arg_list, trunk=NULL, segmentss_nadj=NULL, segmentss=NULL, res_outer=NULL, res_inner=list(res_outer) )
			res_info$full_tree = full_tree
			res_folder_final = file.path(res_dir, 'final')
			dir.create(res_folder_final, recursive=TRUE, showWarnings = TRUE)
			save(res_info, file=file.path(res_folder_final, 'res_info.Rdata'))
	
			time_finish = Sys.time()
			cat('Execution finishes:', as.character(time_finish), '\n\n', file=total_execution_time_file, append=TRUE)
			cat('Total execution time:', capture.output( time_finish - time_begin ), '\n\n', file=total_execution_time_file, append=TRUE)
			return(res_info)
		}
		
		## which one leads to the least residue
		max_nbins = allowed_max_nbins_seq[ min(which(n_A%%allowed_max_nbins_seq==min(n_A%%allowed_max_nbins_seq))) ]
		residue = n_A %% max_nbins
	
		residue_mat_info = get_least_residue_matrix(A, max_nbins=max_nbins, allowed_shifts=0)
		A_final = residue_mat_info$A_final
		n2one = residue_mat_info$n2one

		## A_low has no row/col name
		A_low = HighResolution2Low_k( A_final, k=max_nbins ) ## k is the dimension of A_final_low
		
		## Starting from here, the original row/col names are no longer used
		
		## this part run HRG_core on the outer part
		{
			dists = (2*min_n_bins_outer - 1):(nrow(A_low)-1)
			counts = sapply(dists, function(v) n_cells2compute( A_low, v, min_n_bins_outer ))
			split_info = split2BanlancedChunks_min_max(vec=counts, K=min(length(counts)-1, ncores))
			chunks = c(dists[1]-1, dists[split_info$break_points])
			res_outer = HRG_core( A=A_low, chunks=chunks, res_folder=NULL, min_n_bins=min_n_bins_outer, distr=distr, fast=fast )
		}

		## this part gets the hi_tree
		{
			hi_tree = get_tree_v0( res_outer )
			igraph::V(hi_tree)$left_rescaled = n2one*(igraph::V(hi_tree)$left - 1) + 1
			igraph::V(hi_tree)$right_rescaled = n2one*igraph::V(hi_tree)$right
			node2adjust = which(igraph::V(hi_tree)$right_rescaled==max(igraph::V(hi_tree)$right_rescaled)) ## append the residue bins here
			igraph::V(hi_tree)[node2adjust]$right_rescaled = n_A
			
			igraph::V(hi_tree)$width_rescaled = igraph::V(hi_tree)$right_rescaled - igraph::V(hi_tree)$left_rescaled + 1
			igraph::V(hi_tree)$name = paste('(',igraph::V(hi_tree)$left_rescaled, ',', igraph::V(hi_tree)$right_rescaled, ')', sep='')
		}

		# segmentss = get_segments(hi_tree, binsize_thresh=max_nbins_bp)	

		segmentss_nadj = get_segments(hi_tree, binsize_thresh=max_nbins_fine)
		
		## THIS PARAT ADJUST THE segmentss BASED ON TOPDOM BIN_SIGNAL
		## 14-05-2018
		## Also adjust the tree
		if(adjust_segmentss_topdom==TRUE)
		{
			shift = seg[1]-1
			segmentss = segmentss_adjust_topdom(pA_sym, segmentss_nadj+shift, n2one, ws=5:20) - shift
			for( i in 1:nrow(segmentss) )
			{
				index_left2adj = which(igraph::V(hi_tree)$left_rescaled == segmentss_nadj[i,1])
				index_right2adj = which(igraph::V(hi_tree)$right_rescaled == segmentss_nadj[i,2])
				igraph::V(hi_tree)[index_left2adj]$left_rescaled = segmentss[i,1]
				igraph::V(hi_tree)[index_right2adj]$right_rescaled = segmentss[i,2]
				igraph::V(hi_tree)$width_rescaled = igraph::V(hi_tree)$right_rescaled - igraph::V(hi_tree)$left_rescaled + 1
				igraph::V(hi_tree)$name = paste('(',igraph::V(hi_tree)$left_rescaled, ',', igraph::V(hi_tree)$right_rescaled, ')', sep='')
			}
		}
		
		trunk = hi_tree
		# save(segmentss, file=file.path(res_dir, 'outer', 'segmentss.Rdata'))
		
		res_inner = rep(list(list()), nrow(segmentss))
		for( i in 1:nrow(segmentss) )
		{
			cat(i,'\n')
			# if(i==2) while(as.numeric(substr(as.character(Sys.time()),12,13))!=20) {}
			index = segmentss[i,1]:segmentss[i,2]
			A_part = A[index, index]
			dists = (2*min_n_bins_inner - 1):(nrow(A_part)-1)
			counts = sapply(dists, function(v) n_cells2compute( A_part, v, min_n_bins_inner ))
			split_info = split2BanlancedChunks_min_max(vec=counts, K=min(length(counts)-1, ncores))
			chunks = c(dists[1]-1, dists[split_info$break_points])
			
			## A_part is named
			res_info_inner = HRG_core( A=A_part, chunks=chunks, res_folder=NULL, min_n_bins=min_n_bins_inner, distr=distr, fast=fast )
			res_inner[[i]] = res_info_inner      #minimum size filtering
			cat('finished', '\n')
		}
		
		res_info = list( arg_list=arg_list, trunk=trunk, segmentss_nadj=segmentss_nadj, segmentss=segmentss, res_outer=res_outer, res_inner=res_inner )
		res_folder_final = file.path(res_dir, 'final')
		dir.create(res_folder_final, recursive=TRUE, showWarnings = TRUE)
		save(res_info, file=file.path(res_folder_final, 'res_info.Rdata'))
	
		time_finish = Sys.time()
		cat('Execution finishes:', as.character(time_finish), '\n\n', file=total_execution_time_file, append=TRUE)
		cat('Total execution time:', capture.output( time_finish - time_begin ), '\n\n', file=total_execution_time_file, append=TRUE)


		## xenocfraf:
		# res_inner = res_inner[sapply(res_inner, length) > 0]
		branches = lapply( res_inner, get_tree_v0 )
		for( i in 1:length(branches) ) branches[[i]] = update_branch_name(branches[[i]], root_start=segmentss[i,1])
		
		full_tree = xenocraft( trunk, branches )
		names_tmp = do.call(rbind, strsplit(igraph::V(full_tree)$name, ','))
		igraph::V(full_tree)$left = substring(names_tmp[,1], 2)
		igraph::V(full_tree)$right = substring(names_tmp[,2], 1, nchar(names_tmp[,2])-1)
		if(!is_binary_tree(full_tree)) stop("Trunk + branches do not produce a binary tree")
		
		res_info$full_tree = full_tree
		save(res_info, file=file.path(res_folder_final, 'res_info.Rdata'))
	
		return( res_info )
	}
	
