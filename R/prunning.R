trim_tree_adaptive_top_down_v2 = function( tree, wilcox_p_thresh, mean_diff_thresh )
        {
                # leaves = get_leaves(tree)
                if(igraph::vcount(tree)==1) return(tree)
                # cat('I am in trim_tree_adaptive_top_down_v2\n')
                nsig_nodes = union( igraph::V(tree)[which(igraph::V(tree)$wilcox_p > wilcox_p_thresh)]$name, igraph::V(tree)[which(igraph::V(tree)$mean_diff > mean_diff_thresh)]$name )


                children_of_nsig = names(unlist(igraph::ego(tree, order=1, node=nsig_nodes, mode='out', mindist=1)))
                if(length(children_of_nsig)!=0) trimmed_tree = tree - children_of_nsig
                if(length(children_of_nsig)==0) trimmed_tree = tree

                comps = igraph::decompose(trimmed_tree)
                root_index = which(sapply( comps, function(comp) igraph::V(tree)[1]$name %in% igraph::V(comp)$name )==1)
                trimmed_tree = comps[[root_index]]
                if(!is_binary_tree( trimmed_tree )) stop("trim_tree_adaptive_top_down, the resulted tree is not a binary tree")
                return( trimmed_tree )
        }


prunning = function(branches, p0, to_correct=FALSE, width_thresh=-Inf, width_thresh_CD=5, boundary_signal_thresh=-Inf, return_which='TADs', top_down=FALSE, all_levels=FALSE, CD_border_adj=FALSE, peak_thresh=NULL, mean_diff_thresh)
	{

		size2correct = sum(sapply(branches, igraph::vcount)) - sum(sapply(branches, function(v) length(get_leaves(v))))
		p_thresh = p0/size2correct
		if(to_correct==FALSE) p_thresh=p0
		# if(!is.null(p0)) p_thresh = p0

		
		if(top_down==FALSE)
		{
			trimmed_branches = lapply( branches, trim_tree_adaptive, max_imp_p=p_thresh, max_nimp_p=Inf, width_thresh=width_thresh, boundary_signal_thresh=boundary_signal_thresh, peak_thresh=peak_thresh )
			# size2correct = sum(sapply(trimmed_branches, igraph::vcount)) - sum(sapply(trimmed_branches, function(v) length(get_leaves(v))))
			
			# for( i in 1:length( trimmed_branches ) )
			# {
				# trimmed_branch = trimmed_branches[[i]]
				# if(igraph::vcount(trimmed_branch) > 1) trimmed_branches[[i]] = lapply( trimmed_branches[i], trim_tree_adaptive, max_imp_p=p_thresh, max_nimp_p=Inf, width_thresh=width_thresh, boundary_signal_thresh=-1 )[[1]]
			# }
			

		}
		
		# if(top_down==TRUE) trimmed_branches = lapply( branches, trim_tree_adaptive_top_down, max_imp_p=p_thresh, max_nimp_p=Inf, width_thresh=width_thresh, boundary_signal_thresh=boundary_signal_thresh )
		if(top_down==TRUE) trimmed_branches = lapply( branches, function(branch) trim_tree_adaptive_top_down_v2(wilcox_p_thresh=p_thresh, mean_diff_thresh=mean_diff_thresh, tree=branch ))

		if(CD_border_adj==TRUE)
		{
			all_tads = get_adjusted_nested_TADs( trimmed_branches, width_thresh_CD, all_levels )
			return( all_tads )
		}
		
		## get all nested TADs in trimmed_branches
		if(all_levels==TRUE)
		{
			all_tads = data.frame(start_pos=numeric(), end_pos=numeric())
			widths = c(0, sapply(trimmed_branches, function(v) igraph::V(v)[1]$width))
			for(i in 1:length(trimmed_branches)) 
			{
				all_tads_i = get_all_tads_in_a_trimmed_branch(trimmed_branches[[i]], pos_shift=sum(widths[1:i]))
				all_tads = rbind(all_tads, all_tads_i)
			}
			return( all_tads )
		}
					
		if( return_which=='trimmed_branches' ) return( trimmed_branches )
		
		tad_sizes_ind = lapply( trimmed_branches, function(v) get_leaves(v, 'igraph')$width )
		tad_sizes = unlist(tad_sizes_ind)
		# tads = split(1:sum(tad_sizes), rep(seq_along(tad_sizes), tad_sizes))
		end_pos = cumsum(tad_sizes)
		start_pos = c(1, 1 + end_pos[-length(end_pos)])
		tads = data.frame(start_pos=start_pos, end_pos=end_pos)
		return( tads )
	}




	## This function combines prunning with branches of only one node
	prunning_hybrid <- function(branches, ...)
	{
		names(branches) = as.character(1:length(branches))
		normal_branches = branches[sapply( branches, function(v) class(v)=='igraph' )]
		unnormal_branches = branches[sapply( branches, function(v) class(v)!='igraph' )] ## that is reprsented as bin_start:bin_end

		trimmed_branches = prunning(normal_branches, return_which='trimmed_branches', ...)

		normal_tad_sizes_ind = lapply( trimmed_branches, function(v) get_leaves(v, 'igraph')$width )
		unormal_tad_sizes_ind = unnormal_branches
		tad_sizes_ind = c(normal_tad_sizes_ind, unormal_tad_sizes_ind)
		tad_sizes_ind = tad_sizes_ind[names(branches)]
		tad_sizes = unlist(tad_sizes_ind)

		# tads = split(1:sum(tad_sizes), rep(seq_along(tad_sizes), tad_sizes))
		end_pos = cumsum(tad_sizes)
		start_pos = c(1, 1 + end_pos[-length(end_pos)])
		tads = data.frame(start_pos=start_pos, end_pos=end_pos)
		return( tads )
	}


	get_all_tads_in_a_trimmed_branch <- function(trimmed_branch, pos_shift)
	{
		res = data.frame( start_pos=igraph::V(trimmed_branch)$left + pos_shift, end_pos=igraph::V(trimmed_branch)$right + pos_shift )
		res = res[order(res[,1], res[,2]), ]
		return(res)
	}
	
	prunning_bottom_up <- function(branches, p0=NULL, width_thresh)
	{
		size2correct = sum(sapply(branches, igraph::vcount)) - sum(sapply(branches, function(v) length(get_leaves(v))))
		
		p_thresh = 0.05/size2correct
		if(!is.null(p0)) p_thresh = p0
		
		trimmed_branches = lapply( branches, trim_tree_adaptive, max_imp_p=p_thresh, max_nimp_p=Inf, width_thresh=width_thresh )

		tad_sizes_ind = lapply( trimmed_branches, function(v) get_leaves(v, 'igraph')$width )
		tad_sizes = unlist(tad_sizes_ind)
		# tads = split(1:sum(tad_sizes), rep(seq_along(tad_sizes), tad_sizes))
		end_pos = cumsum(tad_sizes)
		start_pos = c(1, 1 + end_pos[-length(end_pos)])
		tads = data.frame(start_pos=start_pos, end_pos=end_pos)
		return( tads )
	}
	
	
	trim_tree_adaptive_bottom_up <- function( tree, which_p='imp_p' )
	{
		if(which_p=='imp_p') ps = sort(unique(igraph::V(tree)$imp_p), decreasing=TRUE)
		# if(which_p=='nimp_p') ps = sort(unique(igraph::V(tree)$nimp_p), decreasing=TRUE)
		# if(which_p=='both') ps = sort(unique(pmin(igraph::V(tree)$nimp_p, igraph::V(tree)$imp_p)), decreasing=TRUE)
		
		trimed_tree_current = tree
		trimmed_branch_bottom_up = vector('list', length(ps))
		for(i in 1:length(ps))
		{
			trimed_tree_current = trim_tree_adaptive( tree, L_diff_thresh=-Inf, max_imp_p=ps[i], max_nimp_p=Inf, width_thresh=-Inf )
			trimmed_branch_bottom_up[[i]] = trimed_tree_current
		}
		igraph::vcounts = sapply(trimmed_branch_bottom_up, igraph::vcount)
		ps = ps[!duplicated(igraph::vcounts)]
		trimmed_branch_bottom_up = trimmed_branch_bottom_up[!duplicated(igraph::vcounts)]
		res = list(ps=ps, trimmed_branch_bottom_up=trimmed_branch_bottom_up)
		return( res )
	}
	
	
	## get adjusted nested TADs
	get_adjusted_nested_TADs <- function( trimmed_branches, width_thresh_CD, all_levels )
	{
			widths = c(0, sapply(trimmed_branches, function(v) igraph::V(v)[1]$width))
			all_tads_i_list = lapply( 1:length(trimmed_branches), function(i) get_all_tads_in_a_trimmed_branch(trimmed_branches[[i]], pos_shift=sum(widths[1:i])))

			for(i in 1:length(trimmed_branches)) 
			{
				all_tads_i = all_tads_i_list[[i]]
				if( nrow(all_tads_i) <= 1 ) next
				## move the left-most border a little bit right if needed
				left_borders = unique(all_tads_i[,1])
				min_diff_left = left_borders[2] - left_borders[1]
				if( min_diff_left <= width_thresh_CD  )
				{
					all_tads_i[ all_tads_i==left_borders[1] ] = left_borders[2]
					all_tads_i = all_tads_i[ all_tads_i[,2] > all_tads_i[,1],  ] ## remove "negative" TADs
					all_tads_i = unique(all_tads_i[order(all_tads_i[,1], all_tads_i[,2]), ]) ## reorder the TADs

					all_tads_i_list[[i]] = all_tads_i

					## need to modify the right border of nested TADs in previous CD if the left border of this CD is modified
					if(i > 1)
					{
						## replace the max value of [i-1], i.e., the right most border, as the min of [i]-1, i.e., the left most border of [i]	
						all_tads_i_list[[i-1]][ all_tads_i_list[[i-1]]==max(all_tads_i_list[[i-1]]) ] = min(all_tads_i_list[[i]]) - 1
					}
				}

				if( nrow(all_tads_i) <= 1 ) next
				## move the right-most border a little bit left if needed
				right_borders = unique(rev(all_tads_i[,2]))
				min_diff_right = right_borders[1] - right_borders[2]
				if( min_diff_right <= width_thresh_CD  )
				{
					all_tads_i[ all_tads_i==right_borders[1] ] = right_borders[2]
					all_tads_i = all_tads_i[ all_tads_i[,2] > all_tads_i[,1],  ]
					all_tads_i = unique(all_tads_i[order(all_tads_i[,1], all_tads_i[,2]), ]) ## reorder the TADs
					
					all_tads_i_list[[i]] = all_tads_i

					if(i < length(trimmed_branches))
					{
						## replace the max value of [i-1], i.e., the right most border, as the min of [i]-1, i.e., the left most border of [i]	
						all_tads_i_list[[i+1]][ all_tads_i_list[[i+1]]==min(all_tads_i_list[[i+1]]) ] = max(all_tads_i_list[[i]]) + 1
					}
				}
			}

			if(!all_levels) all_tads_i_list = lapply( all_tads_i_list, function(v) data.frame(start_pos=head(v[,1],1), end_pos=tail(v[,2],1)) )

			all_tads = do.call(rbind, all_tads_i_list)
			colnames(all_tads) = c('start_pos', 'end_pos')
			return( all_tads )
	}
	
	
	
	