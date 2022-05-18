	## Yuanlong LIU, 09-08-2018

	generate_compartments_bed <- function(input_mat, p_thresh, out_file_name, chr, window.sizes=3, stat_window_size=NULL, bin_size)
	{	
		input_mat_extended = data.frame(chr=paste0('chr', chr), pos_start=0:(nrow(input_mat)-1), pos_end=1:nrow(input_mat), mat=input_mat)
		res_input_mat = TopDom_v2(input_mat_extended, window.size=NULL, NULL, T, p_thresh=p_thresh, window.sizes=window.sizes, stat_window_size=stat_window_size, domain_size_min=NULL)
		# return(res_input_mat)
		to_id = as.numeric(rownames(input_mat)[res_input_mat$domain$to.id])
		from_id = as.numeric(rownames(input_mat)[res_input_mat$domain$from.id])
		
		start_poses = (from_id-1)*bin_size + 1
		# end_poses = start_poses + n2one*bin_size
		end_poses = start_poses + bin_size
		
		input_mat_compartments_bed = data.frame(paste0('chr', chr), as.character(start_poses), as.character(end_poses) )
		if(!is.null(out_file_name)) write.table( input_mat_compartments_bed, file=out_file_name, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=' ' )
		return( res_input_mat )
	}
	