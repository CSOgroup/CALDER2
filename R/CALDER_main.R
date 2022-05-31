
    CALDER = function(contact_tab_dump=NULL, contact_file_dump=NULL, contact_file_hic=NULL, chrs, bin_size, save_dir, save_intermediate_data=FALSE, swap_AB=0, genome='others', feature_track=NULL, black_list_bins=NULL, n_cores=1, sub_domains=FALSE, single_binsize_only=FALSE)
    {
    		########################### https://stackoverflow.com/questions/30216613/how-to-use-dopar-when-only-import-foreach-in-description-of-a-package

    		`%dopar%` <- foreach::`%dopar%`
			`%do%` <- foreach::`%do%`

    		###########################

            chrs = as.character(chrs)
            n_chrs = length(chrs)

            if(n_chrs > 1) ## when computing for multiple chrs, contact_tab_dump or contact_file_dump should be a list with elements matching to that of chrs
            {
            	if(!is.null(contact_tab_dump))
            	{
            		if(class(contact_tab_dump)!='list' | length(contact_tab_dump)!=n_chrs) stop('contact_tab_dump should be a list of contact table with its length equal to the length of chrs you specified\n')
                    names(contact_tab_dump) = chrs
            	}

            	if(!is.null(contact_file_dump))
            	{
            		if(class(contact_file_dump)!='list' | length(contact_file_dump)!=n_chrs) stop('contact_file_dump should be a list of contact table with its length equal to the length of chrs you specified\n')
                    names(contact_file_dump) = chrs
            	}
            }

            ########################### try multiple bin_size parameters and choose the one that generates reliable compartments

            if(bin_size==5E3) bin_sizes = c(5E3, 10E3, 50E3, 100E3)
            if(bin_size==10E3) bin_sizes = c(10E3, 50E3, 100E3)
            if(bin_size==20E3) bin_sizes = c(20E3, 40E3, 100E3)
            if(bin_size==25E3) bin_sizes = c(25E3, 50E3, 100E3)
            if(bin_size==40E3) bin_sizes = c(40E3, 80E3)
            if(bin_size==50E3) bin_sizes = c(50E3, 100E3)
            if(!(bin_size %in% c(5E3, 10E3, 20E3, 25E3, 40E3, 50E3))) bin_sizes = bin_size

            if(genome=='others') single_binsize_only = TRUE
            if(single_binsize_only) bin_sizes = bin_size ## do not try multiple bin_sizes

            ## choose the already computed good compartment as reference to decide correctly the A/B compartment direction
            if(genome!='others') ref_compartment_file = system.file("extdata", sprintf('%s_all_sub_compartments.bed', genome), package = 'CALDER')
            if(genome=='others') ref_compartment_file = NULL

            ###########################

            if(is.null(ref_compartment_file) & is.null(feature_track)) stop('Should either specify genome in one of hg19, hg38, mm9, mm10, or provide a feature_track\n')

            ###########################

            doParallel::registerDoParallel(cores=n_cores)

            for(bin_size2look in bin_sizes)
            {
            	bin_size_kb = sprintf("%skb", bin_size2look/1E3)
            	save_dir_binsize = file.path(save_dir, 'intermediate_data/sub_compartments', bin_size_kb)
            	dir.create(save_dir_binsize, recursive=TRUE, showWarnings=FALSE)
            }

            para_tab = data.table::data.table(expand.grid(bin_sizes, chrs))
            for(i in 1:ncol(para_tab)) para_tab[[i]] = as.vector(para_tab[[i]])
            colnames(para_tab) = c('bin_size', 'chr')

            ########################### compute compartment for each bin_size and chr

            silent_out =  foreach::foreach(i=1:nrow(para_tab)) %dopar%
	        {	
	        	bin_size2look = para_tab$bin_size[i]
	           	chr = para_tab$chr[i]
           		bin_size_kb = sprintf("%skb", bin_size2look/1E3)
            	save_dir_binsize = file.path(save_dir, 'intermediate_data/sub_compartments', bin_size_kb)

            	save_intermediate_data_tmp = save_intermediate_data*(bin_size2look==bin_size)

	            CALDER_CD_hierarchy_v2(contact_tab_dump=contact_tab_dump[[chr]], contact_file_dump=contact_file_dump[[chr]], contact_file_hic=contact_file_hic, chr=chr, bin_size_input=bin_size, bin_size2look=bin_size2look, save_dir=save_dir_binsize, save_intermediate_data=save_intermediate_data_tmp, swap_AB=swap_AB, ref_compartment_file=ref_compartment_file, feature_track=feature_track, black_list_bins=black_list_bins)
	        }

            ########################### build optimal compartment from multiple resolutions

            try(build_comp_table_opt(save_dir=save_dir, chrs=chrs, bin_sizes=bin_sizes, with_ref=!is.null(genome)))

            ########################### if computing sub-domains

            if(sub_domains==TRUE) 
            {
            	silence_out = foreach::foreach(chr=chrs) %do% ## use one core instead of multi-cores
            	{
            		chr_num = gsub('chr', '', chr, ignore.case=TRUE)
            		intermediate_data_file = sprintf('%s/%skb/chr%s_intermediate_data.Rds', file.path(save_dir, 'intermediate_data/sub_compartments') , bin_size/1E3, chr_num) 
            		save_dir_sub_domains = file.path(save_dir, 'intermediate_data/sub_domains')
					dir.create(save_dir_sub_domains, showWarnings = FALSE)
            		try(CALDER_sub_domains(intermediate_data_file=intermediate_data_file, chr=chr, save_dir=save_dir_sub_domains, bin_size=bin_size))
            	}

				save_dir_opt = file.path(save_dir, "sub_domains")
				dir.create(save_dir_opt, recursive=TRUE, showWarnings=FALSE)
				save_opt_file = file.path(save_dir_opt, sprintf("all_nested_boundaries.bed"))

				opt_sub_domain_beds = sprintf('%s/chr%s_nested_boundaries.bed', file.path(save_dir, 'intermediate_data/sub_domains') , gsub('chr', '', chrs, ignore.case=TRUE)) 
				cmd_opt = sprintf("cat %s > %s", paste0(opt_sub_domain_beds, collapse=" "), save_opt_file)
				system(cmd_opt)
            }
    }

    
