		## Compute sub-compartments, V2. This version uses already computed compartments as reference to derive the A/B compartment direction
		## Yuanlong LIU, 16-07-2021

		## This code converts the compartment into table format

		make_comp_tab = function(comp_raw, bin_size)
	    {
	        # comp_raw[V1=="chrX"]$V1 = "chr23"
	        # comp_raw[[1]] = gsub('chr', '', comp_raw[[1]])
	        bin_indices_li = apply(comp_raw[, 2:3], 1, function(v) {bins = 1 + seq(v[1]-1, v[2], by=bin_size)/bin_size; bins[1:(length(bins) - 1)]})
	        expand_len = sapply(bin_indices_li, length)

	        n_domains = length(unique(comp_raw[[4]]))

	        continous_rank = (n_domains + 1 - rank(unique(comp_raw[[4]]))) / n_domains
	        names(continous_rank) = unique(comp_raw[[4]])

	        comp_tab = data.table::data.table(chr=rep(comp_raw[[1]], expand_len), bin_index=unlist(bin_indices_li), comp_name=rep(comp_raw[[4]], expand_len), comp_rank=rep(comp_raw[[5]], expand_len))
	        comp_tab$continous_rank = unname(continous_rank[comp_tab$comp_name])

	        comp_tab$pos_start = (comp_tab$bin_index - 1)*bin_size + 1
	        comp_tab$pos_end = (comp_tab$bin_index)*bin_size

	        return(comp_tab)
	    }

		make_comp_tab_calder = function(comp_raw, bin_size)
		{
			`%dopar%` <- foreach::`%dopar%`
			`%do%` <- foreach::`%do%`

			colnames(comp_raw)[1] = 'chr'
		    # comp_raw[V1=="chrX"]$V1 = "chr23"
		    # comp_raw[[1]] = as.numeric(gsub('chr', '', comp_raw[[1]]))

		    make_comp_tab_calder_chr = function(comp_raw_chr, bin_size)
		    {
			    bin_indices_li = apply(comp_raw_chr[, 2:3], 1, function(v) {bins = 1 + seq(v[1]-1, v[2], by=bin_size)/bin_size; bins[1:(length(bins) - 1)]})
			    expand_len = sapply(bin_indices_li, length)

			    n_domains = length(unique(comp_raw_chr[[4]]))

			    continous_rank = (n_domains + 1 - rank(unique(comp_raw_chr[[4]]))) / n_domains
			    names(continous_rank) = unique(comp_raw_chr[[4]]) 

			    comp_tab = data.table::data.table(chr=rep(comp_raw_chr[[1]], expand_len), bin_index=unlist(bin_indices_li), comp_name=rep(comp_raw_chr[[4]], expand_len), comp_rank=rep(comp_raw_chr[[5]], expand_len))
			    comp_tab$continous_rank = unname(continous_rank[comp_tab$comp_name])    

			    comp_tab$pos_start = (comp_tab$bin_index - 1)*bin_size + 1
			    comp_tab$pos_end = (comp_tab$bin_index)*bin_size 
			    return(comp_tab)  
		    }

		    comp_tab_ALL = foreach::foreach(comp_raw_chr=split(comp_raw, by="chr"), .combine=rbind) %do%
		    {
		    	make_comp_tab_calder_chr(comp_raw_chr, bin_size)
		    }

		    return(comp_tab_ALL)
		}

	    generate_domain_2D_bedpe = function(sub_compartment_file, bin_size, bedpe_file, n_sub=8)
		{
			generate_domain_2D_bedpe_chr = function(comp_raw)
			{
				comp_tab = c(comp_raw, bin_size=10E3)

				comp_tab$comp_name = substr(comp_tab$comp_name, 1, log2(n_sub)*2 - 1)
				n_row = nrow(comp_tab)
				comp_tab$is_boudary = c(1*(comp_tab$comp_name[2:n_row] != comp_tab$comp_name[1:(n_row-1)]), 0)

				pos_end = union( comp_tab[is_boudary==1]$pos_end, comp_tab[n_row, ]$pos_end )
				pos_start = union( comp_tab[1, ]$pos_start, pos_end[1:(length(pos_end) -1)] + 1 )
				pos_start = as.character(pos_start)
				pos_end = as.character(pos_end)

				chr_name = comp_tab$chr[1]
				bedpe = data.table::data.table(chr_name, pos_start, pos_end, chr_name, pos_start, pos_end)
				return(bedpe)
			}
			my_skip = as.numeric(strsplit(system(intern=TRUE, sprintf('zgrep -n "chr" %s', sub_compartment_file)), ':')[[1]][1])
			comp_raw_li = split(data.table::fread(sub_compartment_file, skip=my_skip), by="chr")
			bedpe = do.call(rbind, lapply(comp_raw_li, generate_domain_2D_bedpe_chr))

			data.table::fwrite(bedpe, file=bedpe_file, col.names=FALSE, sep="\t")
			return(bedpe)
		}


	


		## cool format is not accepted for the current version

		# cool2mat <- function(cool_file, chr_num, bin_size) ## convert cool format to long format / pos_1, pos_2, val / modified from https://rdrr.io/bioc/HiCcompare/src/R/hicpro2bedpe.R
		# {
		# 	dump <- rhdf5::h5dump(cool_file)
		# 	if(names(dump)=='resolutions')
		# 	{
		# 		cat('\nYour input contact matrix is in mcool format and contains resolutions of:', names(dump$resolutions), '\n')
				
		# 		which2keep = which(bin_size==as.numeric(names(dump$resolutions)))
		# 		dump = dump$resolutions[[which2keep]]
		# 	}

		# 	ids <- data.table::data.table(chr = dump$bins$chrom, start = dump$bins$start, id = seq(1, length(dump$bins$chrom), by = 1))
		# 	ids = ids[chr==chr_num][, c('id', 'start')]
			
		# 	# make sparse matrix

		# 	mat <- data.table::data.table(bin1 = dump$pixels$bin1_id, bin2 = dump$pixels$bin2_id, val = dump$pixels$count)
		# 	mat_chr = mat[(bin1 %in% ids$id) & (bin2 %in% ids$id)] ## keep only cis contacts
		# 	colnames(ids)[2] = 'pos_1'
		# 	contact_mat = left_join(mat_chr, ids, by = c('bin1' = 'id'))
		# 	colnames(ids)[2] = 'pos_2'
		# 	contact_mat <- left_join(contact_mat, ids, by = c('bin2' = 'id'))
		# 	contact_mat = contact_mat[, c('pos_1', 'pos_2', 'val')]
		# 	return(contact_mat)
		# }


		CALDER_CD_hierarchy_v2 = function(contact_tab_straw=NULL, contact_file_straw=NULL, contact_file_hic=NULL, chr, bin_size_input, bin_size2look, save_dir, save_intermediate_data=FALSE, swap_AB, ref_compartment_file, black_list_bins=NULL, annotation_track=NULL)
		{
			chr_num = gsub('chr', '', chr, ignore.case=TRUE)
			chr_name = paste0('chr', chr_num)

		    get_cor_with_ref = function(chr_bin_pc, bin_size, ref_compartment_file=NULL, annotation_track=NULL) ## correlation of PC1 with comp. domain rank of ref_genome
		    {
		        ext_chr_bin_pc = function(chr_bin_pc) ## spand chr_bin_pc using 5kb bin
		        {
		            bin_indices_li = unlist(apply(chr_bin_pc[, 4:5], 1, function(v) {bins = 1 + seq(v[1]-1, v[2], by=5E3)/5E3; list(bins[1:(length(bins) - 1)])}), recursive=FALSE)
		            expand_len = sapply(bin_indices_li, length)
		            chr_bin_pc_ext = data.table::data.table(chr=rep(chr_bin_pc[[1]], expand_len), bin_index=unlist(bin_indices_li),PC1_val=rep(chr_bin_pc$PC1_val, expand_len))

		            chr_bin_pc_ext$pos_start = (chr_bin_pc_ext$bin_index - 1)*5E3 + 1
		            chr_bin_pc_ext$pos_end = (chr_bin_pc_ext$bin_index)*5E3        
		            return(chr_bin_pc_ext)
		        }

		        ## function to generate binning scores // https://divingintogeneticsandgenomics.rbind.io/post/compute-averages-sums-on-granges-or-equal-length-bins/
			 	# annotation_track = data.table::data.table(chr=as.vector(GenomicRanges::seqnames(bw_val)), start=start(bw_val), end=end(bw_val), score=bw_val$score)


			 	get_binned_vals = function(annotation_track_chr, bin_size=5E3)
			 	{
			 		## helper to compute binned average // https://divingintogeneticsandgenomics.rbind.io/post/compute-averages-sums-on-granges-or-equal-length-bins/
			  		binnedMean <- function(bins, numvar, mcolname)
					{
					  stopifnot(is(bins, "GRanges"))
					  stopifnot(is(numvar, "RleList"))
					  stopifnot(identical(seqlevels(bins), names(numvar)))
					  bins_per_chrom <- split(ranges(bins), GenomicRanges::seqnames(bins))
					  sums_list <- lapply(names(numvar),
					      function(seqname) {
					          views <- Views(numvar[[seqname]],
					                         bins_per_chrom[[seqname]])
					          viewMeans(views)
					      })
					  new_mcol <- unsplit(sums_list, as.factor(GenomicRanges::seqnames(bins)))
					  mcols(bins)[[mcolname]] <- new_mcol
					  bins
					}


			 		GR = GenomicRanges::makeGRangesFromDataFrame(annotation_track_chr, keep.extra.columns=TRUE)
			 		GR_chrs = split(GR, GenomicRanges::seqnames(GR))
			 		seq_lens = sapply(GR_chrs, function(v) max(end(v)))

					GR_RleList = GenomicRanges::coverage(GR, weight="score")
					seq_info = GenomicRanges::seqinfo(GR_RleList)
					GenomeInfoDb::seqlengths(seq_info) = seq_lens

					bins = GenomicRanges::tileGenome(seq_info, tilewidth=bin_size, cut.last.tile.in.chrom=TRUE)
					bins = bins[width(bins)==bin_size]
					bin_val_tmp = binnedMean(bins, GR_RleList, "bin_val")
					bin_val_tmp = data.table::as.data.table(bin_val_tmp)
					bin_val = data.table::data.table(chr=bin_val_tmp$seqnames, bin_index=bin_val_tmp$end / bin_size, continous_rank=log2(1 + bin_val_tmp$bin_val - min(bin_val_tmp$bin_val)))
					return(bin_val)
			 	}


		        chr2query = chr_bin_pc$chr[1]

		        if(!is.null(ref_compartment_file)) 
		        {
		        	domain_rank = data.table::fread(ref_compartment_file)
			        colnames(domain_rank)[1] = 'chr'
			        domain_rank_chr = domain_rank[domain_rank[[1]]==chr2query, ]

			        ref_tab = make_comp_tab_calder(domain_rank_chr, bin_size=5E3) ## 5kb is convenient for most of the bin sizes
			        ref_tab = ref_tab[, c(1,2,5)]
		        }

		        if(!is.null(annotation_track)) 
		        {
		        	colnames(annotation_track) = c('chr', 'start', 'end', 'score')
		        	annotation_track$chr = gsub('chr', '', annotation_track$chr)
		        	annotation_track_chr = annotation_track[annotation_track$chr==chr_num, ]
		        	ref_tab = get_binned_vals(annotation_track_chr)
		        }


		        chr_name = gsub("chrchr", "chr", paste0("chr", ref_tab$chr))
		        # chr_name = ifelse(chr_name=="chr23", "chrX", chr_name)
		        ref_tab$chr = chr_name
		        colnames(ref_tab)[2] = "bin_index"

		        ################################# convert chr_bin_pc into 5kb

		        chr_bin_pc_ext = ext_chr_bin_pc(chr_bin_pc)
		        ref_and_pc = merge(ref_tab, chr_bin_pc_ext, by=c("chr", "bin_index"))
		        cor_with_ref = cor(method='spearman', ref_and_pc[, c("continous_rank", "PC1_val")])[1,2]

		        return(cor_with_ref)
		    }
	
			### The main function starts here
	    

	        time0 = Sys.time()
	        log_file = paste0(save_dir, '/chr', chr_num, '_log.txt')
	        warning_file = paste0(save_dir, '/WARNING_chr', chr_num, '.txt')
	        cor_log_file = paste0(save_dir, '/cor_with_ref.txt')
	        cat('\n')

	        cat('>>>> Begin process contact matrix and compute correlation score at:', as.character(Sys.time()), '\n', file=log_file, append=FALSE)
	        cat('>>>> Begin process contact matrix and compute correlation score at:', as.character(Sys.time()), '\n')
	        processed_data = contact_mat_processing_v2(contact_tab_straw, contact_file_straw=contact_file_straw, contact_file_hic=contact_file_hic, chr=chr_num, bin_size_input=bin_size_input, bin_size2look=bin_size2look, black_list_bins=black_list_bins)
	     

	        mat_dense = processed_data$mat_dense
	        ccmat_dense_compressed_log_atanh = processed_data$atanh_score

	        cat('\r', '>>>> Finish process contact matrix and compute correlation score at:', as.character(Sys.time()))
	        cat('>>>> Finish process contact matrix and compute correlation score at:', as.character(Sys.time()), '\n', file=log_file, append=TRUE)

	        p_thresh = ifelse(bin_size2look < 40000, 0.05, 1)
	        window.sizes = 3
	        compartments = vector("list", 2)
	        chr_name = paste0("chr", chr)

	        cat('>>>> Begin compute compartment domains and their hierachy at:', as.character(Sys.time()), '\n', file=log_file, append=TRUE)
	        cat('\r', '>>>> Begin compute compartment domains and their hierachy at:', as.character(Sys.time()))

	        compartments[[2]] = generate_compartments_bed(chr = chr, bin_size = bin_size2look, window.sizes = window.sizes, ccmat_dense_compressed_log_atanh, p_thresh, out_file_name = NULL, stat_window_size = NULL)
	        topDom_output = compartments[[2]]
	        bin_names = rownames(mat_dense)
	        mat_dense = as.matrix(mat_dense)
	        initial_clusters = apply(topDom_output$domain[, c("from.id", "to.id")], 1, function(v) v[1]:v[2])

	        if (sum(sapply(initial_clusters, length)) != max(unlist(initial_clusters))) {
	            stop(CELL_LINE, " initial_clusters error in topDom")
	        }

	        n_clusters = length(initial_clusters)
       		mat_dense_cluster_mean = HighResolution2Low_k_rectangle(mat_dense, initial_clusters, initial_clusters, sum_or_mean = "mean")
        
			trend_mean_list = lapply( 1:4, function(v) 1*(mat_dense_cluster_mean[, -(1:v)] > mat_dense_cluster_mean[, - n_clusters - 1 + (v:1)]) )
			trend_mean = do.call(cbind, trend_mean_list)
			c_trend_mean = cor(t(trend_mean))
			atanh_c_trend_mean= atanh(c_trend_mean / (1+1E-7))


			# if(to_scale)
			{
				trend_mean = scale(trend_mean)
				c_trend_mean = scale(c_trend_mean)
				atanh_c_trend_mean= scale(atanh_c_trend_mean)
			}


			PC_12_atanh = get_PCs(atanh_c_trend_mean, which=1:10)
			PC_12_atanh[, 2:10] = PC_12_atanh[, 2:10]/5 ## xx-xx-xxxx: compress PC2
			rownames(PC_12_atanh) = 1:nrow(PC_12_atanh)

			############################################################

			PC_direction = 1

			## switch PC direction based on correlation with "ground_truth"	
			{
				initial_clusters_ori_bins = lapply(initial_clusters, function(v) as.numeric(bin_names[v]))
				chr_bin_pc = data.table::data.table(chr=chr_name, bin=unlist(initial_clusters_ori_bins), PC1_val=rep(PC_12_atanh[,1], sapply(initial_clusters_ori_bins, length)))
				chr_bin_pc$start = (chr_bin_pc$bin - 1)*bin_size2look + 1
				chr_bin_pc$end = chr_bin_pc$bin*bin_size2look

				# chr_bin_pc_range = makeGRangesFromDataFrame(chr_bin_pc, keep.extra.columns=TRUE)
				# gene_info_chr = subset(gene_info, seqnames==chr_name)

				# refGR = chr_bin_pc_range
				# testGR = gene_info_chr
				# hits <- findOverlaps(refGR, testGR)
			 	# overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
			 	# overlaps_bins = unique(data.table::data.table(overlap_ratio=width(overlaps)/bin_size, bin=overlaps$bin))
			 	# bin_pc_gene_coverage = merge(chr_bin_pc, overlaps_bins, all.x=TRUE)
				# bin_pc_gene_coverage$overlap_ratio[is.na(bin_pc_gene_coverage$overlap_ratio)] = 0
				
				# gene_density_cor = cor(method='spearman', subset(bin_pc_gene_coverage, (PC1_val < quantile(PC1_val, 0.25)) | (PC1_val > quantile(PC1_val, 0.75)) , c('PC1_val', 'overlap_ratio')))[1,2]
				# if(abs(gene_density_cor) < 0.2) warning('correlation between gene density and PC1 is too weak')
				

				if(!is.null(ref_compartment_file)) cor_with_ref = try(get_cor_with_ref(chr_bin_pc, bin_size2look, ref_compartment_file=ref_compartment_file)) ## get correlation with supplied "ground truth or reference"
				if(is.null(ref_compartment_file)) cor_with_ref = try(get_cor_with_ref(chr_bin_pc, bin_size2look, annotation_track=annotation_track))



				if(class(cor_with_ref)=='try-error') cor_with_ref = 1 ## psudo cor
	        	if(!is.null(ref_compartment_file)) cat('\r', ">>>> Correlation between PC1 and reference compartment is :", format(abs(cor_with_ref), digits=5), '\n')
	        	if(is.null(ref_compartment_file)) cat('\r', ">>>> Correlation between PC1 and annotation_track is :", format(abs(cor_with_ref), digits=5), '\n')

				PC_direction = PC_direction*sign(cor_with_ref)
			    if(swap_AB==1) PC_direction = -PC_direction ## force swap PC direction if in some case the A/B direction is reverted
			    PC_12_atanh = PC_12_atanh*PC_direction
			}


			project_info = project_to_major_axis(PC_12_atanh)
			x_pro = project_info$x_pro
			
			############################################################

			hc_disect_kmeans_PC12 = bisecting_kmeans(PC_12_atanh[, 1:10, drop=FALSE])

			hc_hybrid_PC12 = hc_disect_kmeans_PC12

			{
				reordered_names = reorder_dendro(hc_hybrid_PC12, x_pro, aggregateFun=mean)
				hc_hybrid_PC12_reordered = dendextend::rotate(hc_hybrid_PC12, reordered_names)
				hc_hybrid_x_pro = hc_disect_kmeans_PC12
				reordered_names_x_pro = get_best_reorder(hc_hybrid_x_pro, x_pro)
				CALDER_hc = dendextend::rotate(hc_hybrid_x_pro, reordered_names_x_pro)	
			}

			############################################################
        	parameters = list(bin_size = bin_size2look, p_thresh = p_thresh)
			res = list( CALDER_hc=CALDER_hc, initial_clusters=initial_clusters, bin_names=bin_names, x_pro=x_pro, parameters=parameters)
			intermediate_data_file = paste0(save_dir, '/chr', chr, '_intermediate_data.Rds')
			
			hc = res$CALDER_hc
			hc_k_labels_full = try(get_cluser_levels(hc, k_clusters=Inf, balanced_4_clusters=FALSE)$cluster_labels)
			bin_comp = data.table::data.table(chr=chr, bin_index=res$bin_names, comp=rep(hc_k_labels_full, sapply(res$initial_clusters, length)))

			rownames(bin_comp) = NULL
			res$comp = bin_comp
			res$CDs = lapply(res$initial_clusters, function(v) res$bin_names[v])
			res$mat = mat_dense
			res$chr = chr
			generate_hierachy_bed(chr=chr, res=res, save_dir=save_dir, bin_size=bin_size2look)


	        cat('>>>> Finish compute compartment domains and their hierachy at: ', as.character(Sys.time()), '\n', file=log_file, append=TRUE)
	        cat('\r', '>>>> Finish compute compartment domains and their hierachy at: ', as.character(Sys.time()))

	       	if(!is.null(ref_compartment_file)) cat('Correlation between PC1 and reference compartment domain rank on this chr is: ', format(abs(cor_with_ref), 1, 5), '\n', file=log_file, append=TRUE)
	       	if(is.null(ref_compartment_file)) cat('Correlation between PC1 and annotation_track on this chr is: ', format(abs(cor_with_ref), 1, 5), '\n', file=log_file, append=TRUE)	       	

	       	# cat(sprintf("chr%s: %s\n", chr, format(abs(cor_with_ref), 1, 5)), file=cor_log_file, append=TRUE)
	       	# if(abs(cor_with_ref) < 0.3) cat('WARNING: correlation between PC1 and reference compartment domain rank on this chr is: ', format(cor_with_ref, digits=5), ', which is a bit low. Possible reason could be that this chromosome has some big structural variations (translocation, inversion for example). We suggest to overlay the compartment track with the hic map together with histone modification or gene expression track to double check the reliability of compartment calling on this chr.', '\n', file=warning_file)

	        time1 = Sys.time()
	        # delta_time  = gsub('Time difference of', 'Total time used for computing compartment domains and their hierachy:', print(time1 - time0))

	        delta_time <- time1 - time0
			timediff <- format(round(delta_time, 2), nsmall = 2)

	        cat('\n\n', 'Total time used for computing compartment domains and their hierachy:', timediff, '\n', file=log_file, append=TRUE)
	       	# if(abs(gene_density_cor) > 0.2) cat('The gene density and PC1 correlation on this chr is: ', substr(gene_density_cor, 1, 4), '\n', file=log_file, append=TRUE)

			############################################################
			intermediate_data = res
			if(save_intermediate_data==TRUE) saveRDS(intermediate_data, file=intermediate_data_file)
			# cat(intermediate_data_file)
			return(intermediate_data_file)
		}
