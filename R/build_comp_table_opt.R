	
	## function to obtain the optimal compartment calling from various resolutions (bin_sizes)

	build_comp_table_opt = function(save_dir, chrs, bin_sizes, with_ref) 
	{

    	`%dopar%` <- foreach::`%dopar%`
		`%do%` <- foreach::`%do%`

		get_opt_index = function(consist_tab)
		{
			if(ncol(consist_tab)==2) ## if only one bin_size value
			{
				opt_index = rep(1, nrow(consist_tab))
				return(opt_index)			
			}

			mins = apply(consist_tab[,-1], 1, max) - 0.05 ## allow at most 0.05 smaller than the maximum
			opt_index = apply(1*((consist_tab[,-1] >= mins)==1), 1, function(v) min(which(v==1)))
			return(opt_index)
		}
		
		################################

		identities = 'bulk'

		{
			consist_tab_li = foreach::foreach(identity=identities) %do%
			{
				# consist_tab = foreach::foreach(bin_size=bin_sizes, .combine=merge) %do%
				# {
				# 	bin_size_kb = sprintf("%skb", bin_size/1E3)
				# 	save_dir_binsize = file.path(save_dir, bin_size_kb)
				#     cor_log_file = paste0(save_dir_binsize, '/cor_with_ref.txt')
				#     log_file = paste0(save_dir_binsize, '/chr', chr, '_log.txt')

				#     as.numeric(strsplit(readLines(log_file)[5], 'this chr is:')[[1]][2])

				#     cor_tab = data.table::fread(cor_log_file)
				#     colnames(cor_tab) = c("chr", bin_size_kb)
				#     cor_tab
				# }


				consist_tab = foreach::foreach(bin_size=bin_sizes, .combine=merge) %do%
				{
					bin_size_kb = sprintf("%skb", bin_size/1E3)
            		save_dir_binsize = file.path(save_dir, 'intermediate_data/sub_compartments', bin_size_kb)
					consist_tab_tmp = data.table::data.table(chr=paste0('chr', chrs), val=0)
					consist_tab_tmp$val = foreach::foreach(chr=chrs, .combine=c) %do%
					{
					    log_file = paste0(save_dir_binsize, '/chr', chr, '_log.txt')
					    # print(log_file)
					    cor_val = as.numeric(strsplit(readLines(log_file)[5], 'this chr is:')[[1]][2])
					    # print(cor_val)
					    cor_val
					}
					colnames(consist_tab_tmp)[2] = bin_size_kb
					consist_tab_tmp
				}

				s = gsub('chr', '', consist_tab[['chr']])

				x <- suppressWarnings(as.numeric(s)) ## https://stackoverflow.com/questions/70080294/sort-column-in-r-strings-first-alphabetically-then-numbers-numerically
				consist_tab = consist_tab[order(x, 'is.na<-'(s, !is.na(x))), ]
				# print((consist_tab))


			}

			names(consist_tab_li) = identities

			min_consist = (sapply(consist_tab_li, function(v) min(v[,2])))

			################################ Choose the best bin size (as small as possible, and save the chosen comp to the opt dir)

			# silent_out = foreach::foreach(identity=identities, .combine=merge) %do% ## do not use dopar
			{

				save_dir_opt = file.path(save_dir, "sub_compartments")
				dir.create(save_dir_opt, recursive=TRUE, showWarnings=FALSE)
				
				consist_tab = consist_tab_li[[identity]]		
				opt_index = get_opt_index(consist_tab)
				names(opt_index) = gsub(":", "", consist_tab$chr)
				opt_bin_tab = data.table::data.table(chr=names(opt_index),  opt_binsize=(colnames(consist_tab)[-1])[opt_index])

				consist_tab_tmp = cbind( consist_tab, opt_bin_tab[, 'opt_binsize'])
				colnames(consist_tab_tmp)[ncol(consist_tab_tmp)] = 'opt_binsize'

				# consist_tab_tmp$chr_num = gsub(':', '', gsub("chr", "", consist_tab_tmp$chr))
				# consist_tab_tmp[chr=="chrX:"]$chr_num = 23
				# consist_tab_tmp = consist_tab_tmp[order(chr_num)]

				cor_all_file = file.path(save_dir_opt, "cor_with_ref.ALL.txt")
				data.table::fwrite(consist_tab_tmp, file=cor_all_file, col.names=TRUE, sep="\t")



				cor_with_ref = foreach::foreach(j=1:nrow(opt_bin_tab), .combine=c) %do%
				{
					chr_query = opt_bin_tab$chr[j]
					opt_binsize = opt_bin_tab$opt_binsize[j]

					row_index = which(consist_tab$chr == chr_query)
					col_index = which(colnames(consist_tab)==opt_binsize)
					as.data.frame(consist_tab)[row_index, col_index] ## some wired behavior when using data.table. Convert to data.frame
				}

				
				cor_with_ref_tab = data.table::data.table(chr=consist_tab$chr, cor=cor_with_ref, chr_num=gsub("chr", "", opt_bin_tab$chr))

				s = cor_with_ref_tab$chr_num
				x <- suppressWarnings(as.numeric(s)) ## https://stackoverflow.com/questions/70080294/sort-column-in-r-strings-first-alphabetically-then-numbers-numerically
				cor_with_ref_tab = cor_with_ref_tab[order(x, 'is.na<-'(s, !is.na(x))), ]


				cor_opt_file = file.path(save_dir_opt, "cor_with_ref.txt")
				data.table::fwrite(cor_with_ref_tab[, 1:2], file=cor_opt_file, col.names=TRUE, sep="\t")


				################################ make plot

				cor_with_ref_tab$index = 1:nrow(cor_with_ref_tab)
				cor_opt_plot_file = file.path(save_dir_opt, "cor_with_ref.pdf")


			    pdf(cor_opt_plot_file, width=8, height=6)
			    p = ggplot2::ggplot(data=cor_with_ref_tab, ggplot2::aes(x=index, y=cor, label = format(round(cor_with_ref_tab$cor, 2), nsmall = 2))) + ggplot2::geom_line(color="red") + ggplot2::geom_point()
			    p = p + ggplot2::geom_label() + ggplot2::scale_x_continuous(labels=cor_with_ref_tab$chr_num, breaks=1:nrow(cor_with_ref_tab))
			    p = p + ggplot2::geom_hline(yintercept=ifelse(with_ref, 0.4, 0.15), linetype=2, color='darkblue')
			    if(with_ref) p = p + ggplot2::xlab('chr') + ggplot2::ylab('Rho') + ggplot2::ggtitle('Correlation of compartment rank with reference\nRho < 0.4 indicates imprecise compartment calling')
			    if(!with_ref) p = p + ggplot2::xlab('chr') + ggplot2::ylab('Rho') + ggplot2::ggtitle('Correlation of compartment rank with reference\nRho < 0.15 indicates imprecise compartment calling')
			    print(p)    
			    dev.off()			

			    ################################

				opt_sub_comp_beds = foreach::foreach(i=1:nrow(opt_bin_tab)) %do%
				{
					opt_bin_size_kb = opt_bin_tab[[2]][i]
					chr_name = opt_bin_tab[[1]][i]
					opt_sub_comp_bed_chr = sprintf("%s/%s_sub_compartments.bed", file.path(save_dir, 'intermediate_data/sub_compartments', opt_bin_size_kb), chr_name)
					opt_sub_comp_bed_chr
				}

				# save_opt_file = file.path(save_dir_opt, "all_sub_compartments.bed")
				save_opt_file = file.path(save_dir_opt, sprintf("all_sub_compartments.bed"))

				cmd_opt = sprintf("cat %s > %s", paste0(opt_sub_comp_beds, collapse=" "), save_opt_file)
				system(cmd_opt)

				# cmd_replaceX = sprintf("sed -i 's/chrNA/chrX/g' %s", save_opt_file) ## the bed file contains NA because of chrX
				# system(cmd_replaceX)

				comp_raw = data.table::fread(save_opt_file)
				comp_tab = make_comp_tab(comp_raw, bin_size=10E3)


				cols2keep = c('chr', 'pos_start', 'pos_end', 'comp_name', 'comp_rank', 'continous_rank')
				comp_tab = comp_tab[, c('chr', 'pos_start', 'pos_end', 'comp_name', 'comp_rank', 'continous_rank')]
				# comp_tab = comp_tab[, mget(cols2keep)] ## mget does not work for some reasons

				save_opt_tab_file = file.path(save_dir_opt, sprintf("all_sub_compartments.tsv"))
				data.table::fwrite(comp_tab, file=save_opt_tab_file, sep='\t')
			}
		}
	}

	
	