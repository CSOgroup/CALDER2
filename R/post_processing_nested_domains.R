
	get_cluster_index <- function(pos, initial_clusters, bin_names, bin_size)
	{
		compartment_segs = generate_compartment_segs( initial_clusters )
		CDs = get_original_tad_indices( bin_names, compartment_segs, bin_size=bin_size )
		cluster_index = which(apply( CDs, 1, function(v) (v[1] < pos)&(v[2]>pos)  )==1)
		return(cluster_index)
	}

	get_original_tad_indices_extra <- function(names_A_final, TADs, bin_size) ## to get the TADs for the extra edges
    {
        # start_pos = as.numeric(names_A_final[TADs$start_pos])
        end_pos = as.numeric(names_A_final[TADs$end_pos])
        # start_pos_ori = (start_pos - 1)*bin_size + 1
        end_pos_ori = end_pos*bin_size
        TADs = data.frame( start_pos_ori=end_pos_ori+1, end_pos_ori=end_pos_ori )
        return( TADs )
    }

	
	LikelihoodRatioTest <- function(sub_domains_raw, ncores, remove_zero=FALSE, distr, n_parameters=3, imputation_num=1E2, A_already_corrected=FALSE)
	{
    	`%dopar%` <- foreach::`%dopar%`
		`%do%` <- foreach::`%do%`

		pA_sym = sub_domains_raw$pA_sym
		if(A_already_corrected==TRUE) cpA_sym = pA_sym
		## CAN SPEEDUP correct_A_fast_divide_by_mean BECAUSE ONLY SOME OFF-DIAGNAL LINES NEED TO BE CORRECTED
		if(A_already_corrected==FALSE) cpA_sym = correct_A_fast_divide_by_mean(pA_sym, remove_zero=remove_zero) ## corrected pA_sym
		trees = foreach::foreach(j=1:length( sub_domains_raw$res_inner )) %do% 
		{
			name_index = rownames(sub_domains_raw$pA_sym)[sub_domains_raw$segmentss[j,1]:sub_domains_raw$segmentss[j,2]]
			sub_domains_raw$res_inner[[j]]$cA = as.matrix(cpA_sym[name_index, name_index]) ## the corrected A. Use this matrix is essential for reducing the distance-based false positves
			
			## for example, when the number of invovled bins in sub_domains_raw$res_inner[[j]] is too small
			# tmp = try(get_tree_decoration( sub_domains_raw$res_inner[[j]], distr=distr, n_parameters=n_parameters, imputation_num=imputation_num ))
			tmp = get_tree_decoration( sub_domains_raw$res_inner[[j]], distr=distr, n_parameters=n_parameters, imputation_num=imputation_num )
			# if( class(tmp)=="try-error" )
			if(class(tmp)!='igraph') ## if(tmp=='bad_tree')
			{
				root_tree = igraph::graph.empty() + 'root'
				igraph::V(root_tree)$width = nrow(sub_domains_raw$res_inner[[j]]$A)
				igraph::V(root_tree)$left = 1
				igraph::V(root_tree)$right = igraph::V(root_tree)$width
				tmp = root_tree ## represented by the length of the segment
			}
			tmp
		}
		return(trees)	
	}

	
	# pipline <- function(chr, TADs_raw, TADs_type, just_save_no_nest=FALSE, bin_size) ## TADs_tmp: start end
	create_TADs <- function(sub_domains_raw, chr, TADs_raw, TADs_type, just_save_no_nest=FALSE, bin_size) ## TADs_tmp: start end
	{
		chr_name = paste0('chr', chr)
		bin_size_kb = bin_size / 1E3

		if(just_save_no_nest)
		{
			TADs = get_original_tad_indices( rownames(sub_domains_raw$pA_sym), TADs_raw, bin_size=bin_size_kb*1E3 )
			TADs = cbind(chr_name, TADs)
			save( TADs, file=Tads_R_File )
		}
		TADs = get_original_tad_indices( rownames(sub_domains_raw$pA_sym), TADs_raw, bin_size=bin_size_kb*1E3 )
		# if(TADs_type=='extra') TADs = get_original_tad_indices_extra( rownames(sub_domains_raw$pA_sym), TADs_raw, bin_size=bin_size_kb*1E3 )

		# print(TADs)

		TADs = cbind(chr_name, TADs)
		return(TADs)
	}

	#####################################################################

	## INDEED, WHEN RESOLUTION == 40KB, DO NOT REMOVE
	clean_TADs_all = function(TADs_all_raw, CDs, bin_size)
	{
		TADs_all_end = TADs_all_raw[,2]
		CD_end = CDs[,2]
		outer_mat = outer(TADs_all_end, CD_end, "-")
		min_dist = apply(outer_mat, 1, function(v) min(abs(v)))
		# end2rm = TADs_all_end[(min_dist <= 3) & (min_dist > 0)] ## remove nested boundaris too close to CD boundary, at least 40kb		
		end2rm = TADs_all_end[(min_dist < 40E3/bin_size) & (min_dist > 0)] ## remove nested boundaris too close to CD boundary, at least 40kb		
		## INDEED, WHEN RESOLUTION == 40KB, DO NOT REMOVE


		TADs_all_head = TADs_all_raw[,1]
		CD_head = CDs[,1]
		outer_mat = outer(TADs_all_head, CD_head, "-")
		min_dist = apply(outer_mat, 1, function(v) min(abs(v)))
		head2rm = TADs_all_head[(min_dist < 40E3/bin_size) & (min_dist > 0)] ## remove nested boundaris too close to CD boundary	
		TADs_all = subset(TADs_all_raw, (!(start_pos %in% head2rm)) & (!(end_pos %in% end2rm)) )	
		return(TADs_all)
	}

	post_process_sub_domains = function(chr, sub_domains_raw, ncores, out_dir, bin_size)
	{
		
		distr=c('lnorm', 'wilcox')[2]
		remove_zero = FALSE
		n_parameters = 3
		imputation_num = 1E2
		A_already_corrected = FALSE

		decorated_branches = LikelihoodRatioTest(sub_domains_raw=sub_domains_raw, ncores=ncores, remove_zero=FALSE, distr=distr, n_parameters=n_parameters, imputation_num=imputation_num, A_already_corrected=A_already_corrected)
		chr_name = paste0('chr', chr)


		mean_diff_thresh = -0.1
		i = 1
		# for(i in 1:length(p0s))
		{
			# cat(i, '\n')
			# p0 = p0s[i]

			normal_decorated_branches = decorated_branches[sapply(decorated_branches, igraph::vcount) > 1]
			# ps = sort(unlist(lapply(normal_decorated_branches, function(v) igraph::V(v)$wilcox_p)), decreasing=FALSE) ## fdr correction
			# p0_adj = ps[min(which(p.adjust(ps, method = 'fdr') > p0))]
			# p0_adj = ps[min(which(p.adjust(ps, method = 'bonferroni') > p0))]
			p0_adj = 0.05

			TADs_all_raw = prunning(decorated_branches, to_correct=FALSE, p0=p0_adj, width_thresh=2, width_thresh_CD=width_thresh_CD, all_levels=TRUE, CD_border_adj=FALSE, top_down=TRUE, mean_diff_thresh=mean_diff_thresh)

			# CDs = prunning(decorated_branches, to_correct=FALSE, p0=p0_adj, width_thresh=2, width_thresh_CD=width_thresh_CD, all_levels=FALSE, CD_border_adj=FALSE, top_down=TRUE); cat(nrow(CDs), '\n')
			CDs = prunning(decorated_branches, to_correct=FALSE, p0=-1, width_thresh=2, width_thresh_CD=width_thresh_CD, all_levels=FALSE, CD_border_adj=FALSE, top_down=TRUE, mean_diff_thresh=mean_diff_thresh)
			# if(nrow(CDs)!=length(res$initial_clusters)) stop('nrow(CDs)!=length(res$initial_clusters)')

			TADs_all = clean_TADs_all(TADs_all_raw, CDs, bin_size=bin_size)

			# Tad_edges <- sort(unique(c(TADs_tmp[,1], (TADs_tmp[,2]+1))))
			TADs_all_edges <- sort(unique(TADs_all[,2]))
			CD_edges <- sort(unique(CDs[,2]))
			
			TADs_extra_edges = unique(setdiff(TADs_all_edges, CD_edges))
			TADs_extra = data.frame(start_pos=TADs_extra_edges, end_pos= TADs_extra_edges)


			TADs_pos_all = create_TADs(sub_domains_raw=sub_domains_raw, chr, TADs_all, TADs_type='ALL', bin_size=bin_size)
			TADs_pos_extra = create_TADs(sub_domains_raw=sub_domains_raw, chr, TADs_extra, TADs_type='extra', bin_size=bin_size)		
			TADs_pos_CD = create_TADs(sub_domains_raw=sub_domains_raw, chr, CDs, TADs_type='CD', bin_size=bin_size)

			# TADs_info = list(decorated_branches=decorated_branches, TADs_pos_all=TADs_pos_all, TADs_pos_CD=TADs_pos_CD)
			TADs_info = list(decorated_branches=decorated_branches, TADs_all=TADs_all, CDs=CDs, TADs_pos_all=TADs_pos_all, TADs_pos_CD=TADs_pos_CD, TADs_pos_extra=TADs_pos_extra)
		}


		options(stringsAsFactors=FALSE)
	    level_5 = TADs_info$TADs_pos_extra[, c('chr_name', 'end_pos_ori')]
	    colnames(level_5)[1:2] = c('chr', 'nested_boundary')
	    sub_domain_boundary_bed_file = paste0(out_dir, '/chr', chr, '_nested_boundaries.bed')

	    level_5_bed = data.frame(level_5, level_5[,2], '', '.', level_5, level_5[,2], '#000000')
	    op <- options(scipen=999)
	    write.table( level_5_bed, file=sub_domain_boundary_bed_file, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE )

	    on.exit(options(op))

		return(NULL)
	}


