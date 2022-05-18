############################################################


get_gene_info <- function(genome) {

	if(genome == 'hg19') {
		gene_info_file = system.file("extdata", "TxDb.Hsapiens.UCSC.hg19.knownGene.rds", package = 'CALDER')
	} else if(genome == 'mm9') {
		gene_info_file = system.file("extdata", "TxDb.Mmusculus.UCSC.mm9.knownGene.rds", package = 'CALDER')
	} else {
		stop(paste0("Unknown genome (", genome, ")"))
	}

	
	gene_info = readRDS(gene_info_file)
}


CALDER_CD_hierarchy = function(contact_mat_file, 
							   chr, 
							   bin_size, 
							   save_dir, 
							   save_intermediate_data=FALSE,
							   genome = 'hg19')
{
    time0 = Sys.time()
    log_file = paste0(save_dir, '/chr', chr, '_log.txt')

    cat('\n')

    cat('>>>> Begin process contact matrix and compute correlation score at:', as.character(Sys.time()), '\n', file=log_file, append=FALSE)
    cat('>>>> Begin process contact matrix and compute correlation score at:', as.character(Sys.time()), '\n')
    processed_data = contact_mat_processing(contact_mat_file, bin_size=bin_size)
   
    A_oe = processed_data$A_oe
    ccA_oe_compressed_log_atanh = processed_data$atanh_score

    cat('\r', '>>>> Finish process contact matrix and compute correlation score at:', as.character(Sys.time()))
    cat('>>>> Finish process contact matrix and compute correlation score at:', as.character(Sys.time()), '\n', file=log_file, append=TRUE)

    p_thresh = ifelse(bin_size < 40000, 0.05, 1)
    window.sizes = 3
    compartments = vector("list", 2)
    chr_name = paste0("chr", chr)
    
    cat('>>>> Begin compute compartment domains and their hierachy at:', as.character(Sys.time()), '\n', file=log_file, append=TRUE)
    cat('\r', '>>>> Begin compute compartment domains and their hierachy at:', as.character(Sys.time()))

    compartments[[2]] = generate_compartments_bed(chr = chr, bin_size = bin_size, window.sizes = window.sizes, ccA_oe_compressed_log_atanh, p_thresh, out_file_name = NULL, stat_window_size = NULL)
    topDom_output = compartments[[2]]
    bin_names = rownames(A_oe)
    A_oe = as.matrix(A_oe)
    initial_clusters = apply(topDom_output$domain[, c("from.id", "to.id")], 1, function(v) v[1]:v[2])

    if (sum(sapply(initial_clusters, length)) != max(unlist(initial_clusters))) {
        stop(CELL_LINE, " initial_clusters error in topDom")
    }

    n_clusters = length(initial_clusters)
		A_oe_cluster_mean = HighResolution2Low_k_rectangle(A_oe, initial_clusters, initial_clusters, sum_or_mean = "mean")

	trend_mean_list = lapply( 1:4, function(v) 1*(A_oe_cluster_mean[, -(1:v)] > A_oe_cluster_mean[, - n_clusters - 1 + (v:1)]) )
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

	gene_info <- get_gene_info(genome)

	## switch PC direction based on gene density
	{
		initial_clusters_ori_bins = lapply(initial_clusters, function(v) as.numeric(bin_names[v]))
		chr_bin_pc = data.table::data.table(chr=chr_name, bin=unlist(initial_clusters_ori_bins), PC1_val=rep(PC_12_atanh[,1], sapply(initial_clusters_ori_bins, length)))
		chr_bin_pc$start = (chr_bin_pc$bin - 1)*bin_size + 1
		chr_bin_pc$end = chr_bin_pc$bin*bin_size
		chr_bin_pc_range = makeGRangesFromDataFrame(chr_bin_pc, keep.extra.columns=TRUE)
		gene_info_chr = subset(gene_info, seqnames==chr_name)

		refGR = chr_bin_pc_range
		testGR = gene_info_chr
		hits <- findOverlaps(refGR, testGR)
	    overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
	    overlaps_bins = unique(data.table::data.table(overlap_ratio=width(overlaps)/bin_size, bin=overlaps$bin))
	    bin_pc_gene_coverage = merge(chr_bin_pc, overlaps_bins, all.x=TRUE)
	    bin_pc_gene_coverage$overlap_ratio[is.na(bin_pc_gene_coverage$overlap_ratio)] = 0
		
	    gene_density_cor = cor(method='spearman', subset(bin_pc_gene_coverage, (PC1_val < quantile(PC1_val, 0.25)) | (PC1_val > quantile(PC1_val, 0.75)) , c('PC1_val', 'overlap_ratio')))[1,2]
	    if(abs(gene_density_cor) < 0.2) warning('correlation between gene density and PC1 is too weak')
	    PC_direction = PC_direction*sign(gene_density_cor)

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
	parameters = list(bin_size = bin_size, p_thresh = p_thresh)
	res = list( CALDER_hc=CALDER_hc, initial_clusters=initial_clusters, bin_names=bin_names, x_pro=x_pro, parameters=parameters)
	intermediate_data_file = paste0(save_dir, '/chr', chr, '_intermediate_data.Rds')
	
	hc = res$CALDER_hc
	hc_k_labels_full = try(get_cluser_levels(hc, k_clusters=Inf, balanced_4_clusters=FALSE)$cluster_labels)
	bin_comp = data.table::data.table(chr=chr, bin_index=res$bin_names, comp=rep(hc_k_labels_full, sapply(res$initial_clusters, length)))

	rownames(bin_comp) = NULL
	res$comp = bin_comp
	res$CDs = lapply(res$initial_clusters, function(v) res$bin_names[v])
	res$mat = A_oe
	res$chr = chr
	generate_hierachy_bed(chr=chr, res=res, save_dir=save_dir, bin_size=bin_size)


    cat('>>>> Finish compute compartment domains and their hierachy at: ', as.character(Sys.time()), '\n', file=log_file, append=TRUE)
    cat('\r', '>>>> Finish compute compartment domains and their hierachy at: ', as.character(Sys.time()))

   	if(abs(gene_density_cor) < 0.2) cat('WARNING: correlation between gene density and PC1 on this chr is: ', substr(gene_density_cor, 1, 4), ', which is a bit low', '\n', file=log_file, append=TRUE)

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
	return(intermediate_data)
}

CALDER_sub_domains = function(intermediate_data_file=NULL, intermediate_data=NULL, chr, save_dir, bin_size)
{	
    time0 = Sys.time()
    log_file = paste0(save_dir, '/chr', chr, '_sub_domains_log.txt')

   	cat('\r', '>>>> Begin compute sub-domains at:', as.character(Sys.time()))
   	cat('>>>> Begin compute sub-domains at:', as.character(Sys.time()), '\n', file=log_file, append=FALSE)

	if(is.null(intermediate_data)) intermediate_data = readRDS(intermediate_data_file)
	{

	    if(intermediate_data$chr!=chr) stop('intermediate_data$chr!=chr; check your input parameters\n') 
	    if( !setequal(rownames(intermediate_data$mat), intermediate_data$bin_names) ) stop('!setequal(rownames(intermediate_data$mat), intermediate_data$bin_names) \n')     
	    compartment_segs = generate_compartment_segs( intermediate_data$initial_clusters )

		cat('\r', '>>>> Begin compute sub-domains within each compartment domain at:', as.character(Sys.time()))   			
		cat('>>>> Begin compute sub-domains within each compartment domain at:', as.character(Sys.time()), '\n', file=log_file, append=TRUE)

		sub_domains_raw = HRG_zigzag_compartment_domain_main_fun(intermediate_data$mat, './', compartment_segs, min_n_bins=2)   

	    no_output = post_process_sub_domains(chr, sub_domains_raw, ncores=1, out_dir=save_dir, bin_size=bin_size)

	    cat('>>>> Finish compute sub-domains within each compartment domain at:', as.character(Sys.time()), '\n', file=log_file, append=TRUE)
	    cat('\r', '>>>> Finish compute sub-domains within each compartment domain at:', as.character(Sys.time()), '\n')

	   	time1 = Sys.time()
        # delta_time  = gsub('Time difference of', 'Total time used for computing compartment domains and their hierachy:', print(time1 - time0))
        delta_time <- time1 - time0
		timediff <- format(round(delta_time, 2), nsmall = 2)

        cat('\n\n', 'Total time used for computing sub-domains:', timediff, '\n', file=log_file, append=TRUE)
	}
	# return(NULL)
}



############################################################

create_compartment_bed_v4 = function(chr_bin_domain, bin_size)
{
	# for( chr in chrs )
	{
		v = chr_bin_domain
		# v$intra_domain = as.character(6 - (as.numeric(v$intra_domain))) ## invert the labeling
		# v$intra_domain = names(cols)[(as.numeric(v$intra_domain))]
		v = v[order(v$bin_index), ]


		borders_non_consecutive = which(diff(v$bin_index)!=1)
		borders_domain = cumsum(rle(v$comp)$lengths)
		borders = sort(union(borders_non_consecutive, borders_domain))
		bins = v$bin_index
		to_id = as.numeric(bins[borders])
		from_id = as.numeric(bins[c(1, head(borders, length(borders)-1)+1)])

		pos_start = (from_id-1)*bin_size + 1
		pos_end = to_id*bin_size
		# chr = as.numeric( gsub('chr', '', v$chr) )
		chr = gsub('chr', '', v$chr) ## no need for as.numeric, also makes it compatible with chrX

		compartment_info_tab = data.frame(chr=rep(unique(chr), length(pos_start)), pos_start=pos_start, pos_end=pos_end, domain=v$comp[borders])
	}
	return(compartment_info_tab)
}

############################################################

generate_hierachy_bed = function(chr, res, save_dir, bin_size)
{
	chr_name = paste0('chr', chr)
	# res = reses[[chr_name]][[CELL_LINE]]
	hc = res$CALDER_hc

	hc_k_labels_full = try(get_cluser_levels(hc, k_clusters=Inf, balanced_4_clusters=FALSE)$cluster_labels)
	bin_comp = data.table::data.table(chr=chr, bin_index=as.numeric(res$bin_names), comp=rep(hc_k_labels_full, sapply(res$initial_clusters, length)))
	chr_bin_domain = bin_comp
	chr_bin_domain$chr = paste0('chr', chr_bin_domain$chr)

	# chr_bin_domain = chr_bin_domain[order(bin_index)]

	compartment_info_tab = create_compartment_bed_v4(chr_bin_domain, bin_size=bin_size)

	boundaries = unname(sapply(res$initial_clusters, max))
	boundaries_ori = as.numeric(res$bin_names[boundaries])*bin_size

	compartment_info_tab$is_boundary = 'gap'
	compartment_info_tab[compartment_info_tab$pos_end %in% boundaries_ori, 'is_boundary'] = 'boundary'
	
	colnames(compartment_info_tab)[4] = 'compartment_label'
	compartments_tsv_file = paste0(save_dir, '/chr', chr, '_domain_hierachy.tsv')
	compartments_bed_file = paste0(save_dir, '/chr', chr, '_sub_compartments.bed')
	boundary_bed_file = paste0(save_dir, '/chr', chr, '_domain_boundaries.bed')

	options(scipen=999)
	write.table( compartment_info_tab, file=compartments_tsv_file, quote=FALSE, sep='\t', row.names=FALSE )


	comp_cols = c("#FF0000", "#FF4848", "#FF9191", "#FFDADA", "#DADAFF", "#9191FF", "#4848FF", "#0000FF") 	
	names(comp_cols) = c('A.1.1', 'A.1.2', 'A.2.1', 'A.2.2', 'B.1.1', 'B.1.2', 'B.2.1', 'B.2.2')
	comp_val = (8:1)/8
	names(comp_val) = names(comp_cols)

	comp_8 = substr(compartment_info_tab$compartment_label, 1, 5)

	compartment_bed = data.frame(chr=paste0('chr', compartment_info_tab$chr), compartment_info_tab[, 2:4], comp_val[comp_8], '.', compartment_info_tab[, 2:3], comp_cols[comp_8])
	write.table( compartment_bed, file=compartments_bed_file, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE )

	bounday_bed_raw = subset(compartment_info_tab, is_boundary=='boundary')
	bounday_bed = data.frame(chr=paste0('chr', compartment_info_tab$chr), compartment_info_tab[,3], compartment_info_tab[,3], '', '.', compartment_info_tab[,3], compartment_info_tab[,3], '#000000')
	write.table( bounday_bed, file=boundary_bed_file, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE )
}



project_to_major_axis <- function(PC_12_atanh)
{
	Data = data.frame(x=PC_12_atanh[,1], y=PC_12_atanh[,2])
	Data = Data[order(Data$x),]
	loess_fit <- loess(y ~ x, Data)

	more_x = seq(min(PC_12_atanh[,1]), max(PC_12_atanh[,1]), len=10*length(PC_12_atanh[,1]))
	major_axis = cbind(x=more_x, y=predict(loess_fit, newdata=more_x))
	new_x_axis = cumsum(c(0, sqrt(diff(major_axis[,1])^2 + diff(major_axis[,2])^2))) ## the new xaxis on the curved major_axis

	dis = fields::rdist(PC_12_atanh[, 1:2], major_axis)
	projected_x = new_x_axis[apply(dis, 1, which.min)]
	names(projected_x) = rownames(PC_12_atanh)
	# projected_x = major_axis[apply(dis, 1, which.min)]
	project_info = list(x_pro=projected_x, major_axis=major_axis)
	return(project_info)
}


get_best_reorder <- function(hc_hybrid_x_pro, x_pro)
{
	n = length(x_pro)
	reordered_names_x_pro_list = list()

	reordered_names_x_pro_list[[1]] = reorder_dendro(hc_hybrid_x_pro, (x_pro), aggregateFun=mean) ## here the clusters are assigned into A.1 A.2 B.1 B.2

	best_index = which.max(sapply(reordered_names_x_pro_list, function(v) cor(1:n, unname(rank(x_pro, ties.method='first')[v]))))
	return(reordered_names_x_pro_list[[1]])
}

