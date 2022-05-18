
	get_PCs = function( mat, which ) ## fast way to compute PCs
	{
		PC_mat = crossprod(mat)
		res_eigs_sym = rARPACK::eigs_sym( PC_mat, k=max(which), which = "LM" )
		if(any(res_eigs_sym$values <0)) stop('Non-positive eigenvalues for A^2')
 		PCs = mat%*%(res_eigs_sym$vectors)
		return( PCs[, which] )
	}

	
	PC_compartment <- function(A_oe) ## compute PC and define A/B compartment. Not used currently
	{
		cA_oe = fast_cor(A_oe)
		PC_mat = crossprod(cA_oe)
		res_eigs_sym = rARPACK::eigs_sym( PC_mat, k=2, which = "LM" )
		PC1 = cA_oe%*%res_eigs_sym$vectors[,1]
		PC2 = cA_oe%*%res_eigs_sym$vectors[,2]

		borders = which(diff(1*(PC1 > 0))!=0)
		to_id = as.numeric(rownames(A_oe)[borders])
		from_id = as.numeric(rownames(A_oe)[c(1, head(borders, length(borders)-1)+1)])
			
		start_poses = (from_id-1)*bin_size + 1
		end_poses = to_id*bin_size

		compartment_AB = data.frame(chr=paste0('chr', chr), start_poses=start_poses, end_poses=end_poses, A_or_B=NA, zero=0, dot='.', start_poses_2=start_poses, end_poses_2=end_poses, col=NA)
		compartment_AB[which((1:nrow(compartment_AB))%%2==0), 'A_or_B'] = 'A'
		compartment_AB[which((1:nrow(compartment_AB))%%2==1), 'A_or_B'] = 'B'
			
		compartment_AB[which((1:nrow(compartment_AB))%%2==0), 'col'] = '112,128,144'
		compartment_AB[which((1:nrow(compartment_AB))%%2==1), 'col'] = '255,255,0'

		compartments_bed_files = paste0(sub_folder, '/chr', chr, '_compartments_PCA_AB_2', '.bed')
		write.table( compartment_AB, file=compartments_bed_files, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=' ' )
	}

	PC_compartment_slim <- function(A_oe, downsratio=NULL) ## compute and save the PC values only. Not used currently
	{
		cA_oe = fast_cor(A_oe)
		PC_mat = crossprod(cA_oe)
		res_eigs_sym = rARPACK::eigs_sym( PC_mat, k=2, which = "LM" )
		PC1 = cA_oe%*%res_eigs_sym$vectors[,1]
		PC2 = cA_oe%*%res_eigs_sym$vectors[,2]

		compartment_AB = data.frame(PC1=PC1, bin_names=rownames(A_oe))

		if( is.null(downsratio) ) compartments_AB_file = paste0(sub_folder, '/chr', chr, '_compartments_PCA_AB.Rdata')
		if( !is.null(downsratio) ) compartments_AB_file = paste0(sub_folder, '/chr', chr, '_compartments_PCA_AB_downsratio_', downsratio, '.Rdata')
		
		save( compartment_AB, file=compartments_AB_file )
	}


