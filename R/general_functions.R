    get_ccA_atanh = function(contact_mat_file, compress_size, bin_size)
    {
        contact_mat_raw = data.table::fread(contact_mat_file)
        contact_mat_raw = subset(contact_mat_raw, !is.na(V3))
        contact_mat_raw[,1] = contact_mat_raw[,1]/bin_size
        contact_mat_raw[,2] = contact_mat_raw[,2]/bin_size
        contact_mat = contact_mat_raw

        colnames(contact_mat) = c('pos_1', 'pos_2', 'val')  
        if(!all(contact_mat[[2]] >= contact_mat[[1]])) stop('\nYour provided matrix does not represent an upper triangular matrix!\n\n')

        n_bins_whole = max(max(contact_mat[[1]]), max(contact_mat[[2]])) + 1 ## should +1 because contact_mat index starts from 0 (bin 0 represents: 0-10E3, checked by looking at the juicebox map, 2018-11-19)
        contact_mat_sparse = Matrix(0, nrow=n_bins_whole, ncol=n_bins_whole)
        contact_mat_sparse[cbind(contact_mat[[1]]+1, contact_mat[[2]]+1)] <- contact_mat[[3]]

        rownames(contact_mat_sparse) = colnames(contact_mat_sparse) = as.character( 1:nrow(contact_mat_sparse) )
        
        contact_mat_sparse = Matrix::forceSymmetric(contact_mat_sparse, uplo='U')
        # if(bin_size!=bin_size_initial) contact_mat_sparse = mat_10to40kb( contact_mat_sparse, bin_size, bin_size_initial )
        A = remove_blank_cols(contact_mat_sparse, sparse=TRUE, ratio=zero_ratio) ## has the same rows/cols as A
        if(nrow(A) < 100) A = remove_blank_cols(contact_mat_sparse, sparse=TRUE, ratio=0) ## when all are dense
        while(min(apply(A, 1, sd))==0) ## sometimes after removing the cols / rows, the remained part will all be 0
        {
            A = remove_blank_cols(A, sparse=TRUE, ratio=1E-7) ## has the same rows/cols as A
            if(nrow(A) < 1) stop('ERROR IN GENERATING MEANINGFUL A at the data generating step')
        }

        # #########################################################

        n_bins_A = nrow(A) - nrow(A)%%compress_size
        A_2_compress = A[, 1:n_bins_A]

        bin_names = rownames(A)

        A_compressed = compress_mat_fast( as.matrix(A_2_compress), compress_size=compress_size )
        colnames(A_compressed) = bin_names
        rm(A_2_compress); gc()
        
        A_compressed_sparse = A_compressed
        A_compressed = as.matrix(A_compressed)
        A_compressed_log = log2(A_compressed + 1)
        
        # #########################################################
        cat('Compute correlation matrix\n')

        cA_compressed_log = fast_cor(A_compressed_log)
        ccA_compressed_log = fast_cor(cA_compressed_log)

        cat('Finished computing correlation matrix\n')

        # #########################################################
        # # ccA_compressed_atanh = atanh(ccA_compressed - 1E-7)
        ccA_atanh = atanh(ccA_compressed_log / (1+1E-7))
        A = as.matrix(A)
        rm(A_compressed, A_compressed_sparse, cA_compressed_log, A_compressed_log); gc()
        # #########################################################
        cat('Finished data generation\n')
        return(list(A=A, ccA_atanh=ccA_atanh))
    }


    
    get_detailed_full_tree <- function(res_info)
    {
        # attach(res_info)
        for( i in 1:nrow(res_info$segmentss) )
        {
            index_left2adj = which(igraph::V(res_info$hi_tree)$left_rescaled == res_info$segmentss_nadj[i,1])
            index_right2adj = which(igraph::V(res_info$hi_tree)$right_rescaled == res_info$segmentss_nadj[i,2])
            igraph::V(res_info$hi_tree)[index_left2adj]$left_rescaled = res_info$segmentss[i,1]
            igraph::V(res_info$hi_tree)[index_right2adj]$right_rescaled = res_info$segmentss[i,2]
            igraph::V(res_info$hi_tree)$width_rescaled = igraph::V(res_info$hi_tree)$right_rescaled - igraph::V(res_info$hi_tree)$left_rescaled + 1
            igraph::V(res_info$hi_tree)$name = paste('(',igraph::V(res_info$hi_tree)$left_rescaled, ',', igraph::V(res_info$hi_tree)$right_rescaled, ')', sep='')
        }
        
        ## xenocfraf:
        # res_inner = res_inner[sapply(res_inner, length) > 0]
        ## xenocfraf:
        # res_inner = res_inner[sapply(res_inner, length) > 0]
        branches = lapply( res_info$res_inner, get_tree_v0 )
        for( i in 1:length(branches) ) branches[[i]] = update_branch_name(branches[[i]], root_start=res_info$segmentss[i,1])
        trunk = res_info$hi_tree
        
        full_tree = xenocraft( trunk, branches )
        names_tmp = do.call(rbind, strsplit(igraph::V(full_tree)$name, ','))
        igraph::V(full_tree)$left = substring(names_tmp[,1], 2)
        igraph::V(full_tree)$right = substring(names_tmp[,2], 1, nchar(names_tmp[,2])-1)
        if(!is_binary_tree(full_tree)) stop("Trunk + branches do not produce a binary tree")
        
        detailed_full_tree = tree_germination(full_tree)
        # detach(res_info)
        return(detailed_full_tree)
    }
    
    

    insulation_score_fun = function(A, size=3)
    {
        insulation_score_fun_helper = function(bin_start)
        {
            bin_mid = bin_start + size -1
            bin_end = bin_start + 2*size - 1
            up = sum(A[bin_start:bin_mid, bin_start:bin_mid])
            down = sum(A[(bin_mid+1):bin_end, (bin_mid+1):bin_end])
            inter = 2*sum(A[bin_start:bin_mid, (bin_mid+1):bin_end])
            insulation = 1 - inter / (up + down)
            return(insulation)
        }

        bin_starts = 1:(nrow(A) - 2*size + 1)
        insulations = sapply(bin_starts, function(bin_start) insulation_score_fun_helper(bin_start))
        names(insulations) = rownames(A)[size:(nrow(A) - size)]
        return(insulations)
    }

    

    ## generating the compartment_segs based on A and initial_clusters

    generate_compartment_segs <- function( initial_clusters )
    {
        bins_seq = unname(unlist(initial_clusters))
        if(!all(diff(bins_seq)==1)) stop('check generate_compartment_segs in generate_compartment_segs.R')
        if(bins_seq[1]!=1) stop('check generate_compartment_segs in generate_compartment_segs.R')
        compartment_segs = do.call(rbind, lapply(initial_clusters, function(v) v[c(1, length(v))]))
        compartment_segs = data.frame( start_pos=compartment_segs[,1], end_pos=compartment_segs[,2] )
        return( compartment_segs )
    }

    ## This function gets the nodes at level k
    ## Yuanlong LIU
    ## 02_07_2018
    
    get_level_k_nodes <- function(tree, k)
    {
        nodes = igraph::ego(tree, order=k, nodes=1, mode = "out", mindist = 1)[[1]]
        return(nodes)
    }
    
    

    ## This function gets the TADs at level k
    ## Yuanlong LIU
    ## 18_05_2018
    
    get_level_k_TADs <- function(branches, k)
    {
        tad_sizes_ind = lapply( branches, function(branch) 
        {
            branch_level_k = igraph::induced.subgraph( branch, igraph::V(branch)$depth <= k )
            tad_size_ind = get_leaves(branch_level_k, 'igraph')$width
        })
        
        tad_sizes = unlist(tad_sizes_ind)
        end_pos = cumsum(tad_sizes)
        start_pos = c(1, 1 + end_pos[-length(end_pos)])
        tads = data.frame(start_pos=start_pos, end_pos=end_pos)
        return( tads )      
    }
    
    


    get_least_residue_matrix <- function(pA_sym, max_nbins, allowed_shifts)
    {
        nbins = nrow( pA_sym )
        max_nbins_bp = max_nbins
        # remainders = nbins %% (max_nbins + (-10:10))
        remainders = nbins %% (max_nbins + allowed_shifts)
            
        max_nbins = (max_nbins + allowed_shifts)[ which.min(remainders) ]
        n2one = floor(nbins/max_nbins)
        remainder = nbins %% max_nbins
        if(max_nbins!=max_nbins_bp) warning('Value of max_nbins has been adjusted from ', max_nbins_bp, ' to ', max_nbins, '. This results in ', remainder, ' rows/columns excluded for the analysis.')

        if(remainder!=0) A_final = pA_sym[-(tail(1:nbins, remainder)), -(tail(1:nbins, remainder))]
        if(remainder==0) A_final = pA_sym
        res = list(A_final=A_final, n2one=n2one, max_nbins_new=max_nbins)
    }



    ## merge clusters based on their corr_mat
    ## Yuanlong LIU
    ## 01-07-2018
    
    sim.fun <- function( corr_mat, mod1, mod2, method='mean' )
    {
        if(method=='mean') return(mean(corr_mat[mod1, mod2]))
        if(method=='median') return(median(corr_mat[mod1, mod2]))
        if(method=='max') return(max(corr_mat[mod1, mod2]))
        
    }
    

    ############################################################    
        
    ## this function removes all-zero columns / rows
    rm_zeros <- function(A) ## mat is a symmetric matrix
    {
        zero_indices = apply( A, 2, function(v) all(v==0) )
        if( all(zero_indices==0) ) return(A)
        A = A[!zero_indices, !zero_indices]
        return(A)
    }
    
    my_dist <- function(m) {mtm <- Matrix::tcrossprod(m); sq <- rowSums(m*m);  res = suppressWarnings(sqrt(outer(sq,sq,"+") - 2*mtm)); res[is.nan(res)]=0; diag(res)=0; return(res)}
    
    
    MoC = function(P, Q)
    {
        len_p = length(P)
        len_q = length(Q)
        if((len_p==len_q) & (len_p==1)) return(1)
        grids = expand.grid(1:len_p, 1:len_q)
        vecs = apply( grids, 1, function(x) {u=x[1]; v=x[2]; length(intersect( P[[u]], Q[[v]] ))^2/length(P[[u]])/length(Q[[v]])} )
        res = 1/( sqrt(len_p*len_q) - 1)*(sum(vecs) - 1)
        return(res)
    }

    get_stair_vecs <- function(mat, from, to)
    {
        n = nrow(mat)
        solve <- function(mid) c(mat[1:mid,(mid+1):n])
        lapply(from:to, solve)
    }
    
    correct_A_fast_divide_by_mean <- function(A, remove_zero=TRUE, divide_or_substract='divide', mean_or_median='mean')
    {
        f = function(i, d) d*(d-1)/2+1+(2*d+i)*(i-1)/2
        
        n = nrow(A)
        upper_tri_values = A[upper.tri(A, diag=TRUE)]
        means_mat = indices = A*0
        indices[upper.tri(indices, diag=TRUE)] = 1:((n+1)*n/2)
        
        diag_values = orders = rep(list(list()), n)

        for(d in seq_along(orders)) orders[[d]] = f(i=seq(1, n-d+1), d=d)
        for(d in seq_along(orders)) diag_values[[d]] = upper_tri_values[f(i=seq(1, n-d+1), d=d)]

        if(mean_or_median=='mean'){
            if(remove_zero) means = unlist(lapply( diag_values, function(v) {m=mean(v[v!=0]); v=v*0+m; v[is.na(v)]=0; return(v) } ))
            if(!remove_zero) means = unlist(lapply( diag_values, function(v) {m=mean(v); v=v*0+m; return(v) } ))            
        }

        if(mean_or_median=='median'){
            if(remove_zero) means = unlist(lapply( diag_values, function(v) {m=median(v[v!=0]); v=v*0+m; v[is.na(v)]=0; return(v) } ))
            if(!remove_zero) means = unlist(lapply( diag_values, function(v) {m=median(v); v=v*0+m; return(v) } ))          
        }

        means[means==0] = 1
        means_mat[upper.tri(means_mat, diag=TRUE)] = means[order(unlist(orders))]
        if(divide_or_substract=='divide') A_corrected = A / means_mat
        if(divide_or_substract=='substract') A_corrected = A - means_mat
        A_corrected = as.matrix( Matrix::forceSymmetric(A_corrected) )
        rownames( A_corrected ) = colnames( A_corrected ) = rownames( A )
        return( A_corrected )
    }


    correct_A_fast_equal_mean <- function(A, remove_zero=TRUE, divide_or_substract='divide', mean_or_median='mean')
    {
        ## may also try log normal, since the off-diagnal values follow well a log-normal distribution
        ## commented by Yuanlong LIU, 2018-07-26
        f = function(i, d) d*(d-1)/2+1+(2*d+i)*(i-1)/2
        
        n = nrow(A)
        upper_tri_values = A[upper.tri(A, diag=TRUE)]
        means_mat = indices = A*0
        indices[upper.tri(indices, diag=TRUE)] = 1:((n+1)*n/2)

        diag_values = orders = rep(list(list()), n)

        for(d in seq_along(orders)) orders[[d]] = f(i=seq(1, n-d+1), d=d)
        for(d in seq_along(orders)) diag_values[[d]] = upper_tri_values[f(i=seq(1, n-d+1), d=d)]

        if(mean_or_median=='mean'){
            if(remove_zero) means = unlist(lapply( diag_values, function(v) {m=mean(v[v!=0]); v=v*0+m; v[is.na(v)]=0; return(v) } ))
            if(!remove_zero) means = unlist(lapply( diag_values, function(v) {m=mean(v); v=v*0+m; return(v) } ))            
        }

        if(mean_or_median=='median'){
            if(remove_zero) means = unlist(lapply( diag_values, function(v) {m=median(v[v!=0]); v=v*0+m; v[is.na(v)]=0; return(v) } ))
            if(!remove_zero) means = unlist(lapply( diag_values, function(v) {m=median(v); v=v*0+m; return(v) } ))          
        }

        means[means==0] = 1
        means_mat[upper.tri(means_mat, diag=TRUE)] = means[order(unlist(orders))]
        if(divide_or_substract=='divide') A_corrected = A / means_mat
        if(divide_or_substract=='substract') A_corrected = A - means_mat
        A_corrected = as.matrix( Matrix::forceSymmetric(A_corrected) )
        rownames( A_corrected ) = colnames( A_corrected ) = rownames( A )
        return( A_corrected )
    }

    creat_phylo_object <- function(tree)
    {
        creat_phylo_object_inner <- function(tree)
        {
            twins = igraph::ego(tree, order=1, node=1, mode='out', mindist=1)[[1]]
            branches = igraph::decompose(tree - 1)
            left_branch = branches[[1]]
            right_branch = branches[[2]]
            left_branch_flag = igraph::vcount( left_branch ) == 1
            right_branch_flag = igraph::vcount( right_branch ) == 1
            short_walk_dist = 1
            igraph::diameters = sapply(branches, igraph::diameter)
            long_walk_dist = abs(diff(igraph::diameters)) + 1
            
            if( left_branch_flag ) left_branch_newick = igraph::V(left_branch)$name
            if( !left_branch_flag ) left_branch_newick = creat_phylo_object_inner(left_branch)
            
            if( right_branch_flag ) right_branch_newick = igraph::V(right_branch)$name
            if( !right_branch_flag ) right_branch_newick = creat_phylo_object_inner(right_branch)   

            if( igraph::diameters[1] > igraph::diameters[2] ) tree_newick = paste( '(', left_branch_newick, ':', short_walk_dist, ',', right_branch_newick, ':', long_walk_dist, ')', sep='' )
            if( igraph::diameters[1] <= igraph::diameters[2] ) tree_newick = paste( '(', left_branch_newick, ':', long_walk_dist, ',', right_branch_newick, ':', short_walk_dist, ')', sep='' )
            return( tree_newick )
        }
        tree_newick = creat_phylo_object_inner(tree)
        final_newick = paste(tree_newick, ';', sep='')
        tree = read.tree(text=final_newick)
        return( tree )
    }
    

    dist_max_min <- function( L_diff )
    {
        N = dim(L_diff)[1]
        L_diff_dist = list()
        for( dist in 1:(N-1) )
        {
            maxLs = numeric( N-dist )
            for( row in 1:(N-dist) )
            {
                i = row
                j = row + dist
                maxLs[i] = L_diff[i, j]
            }
            
            L_diff_dist[[dist]] = maxLs
        }
        return(L_diff_dist)
    }
    

    divide_into_groups <- function(A, n_group)
    {
        N = nrow(A)
        dists = 1:(N-1)
        counts = sapply(dists, function(v) n_cells2compute( A, v ))
        average = (ceiling(sum(counts) / n_group)) ## the average computational complexity
        groups <- cumsum(counts) %/% average + 1
        borders = c(0, which(diff(groups)!=0))
        binsizes = c()
        for( i in 1:(length(borders)-1) )
        {
            binsizes = c(binsizes, sum(counts[ (borders[i]+1):borders[i+1] ]))
        }
        info = list(borders=borders, binsizes=binsizes)
        return(info)
    }
    

    get_leaves <- function(tree, type='name')
    {
        leaf_indices = which(igraph::degree( tree, mode='out' ) == 0)
        if( type=='name' ) return(igraph::V(tree)[leaf_indices]$name)
        if( type=='igraph_node' ) return(igraph::V(tree)[leaf_indices])
        if( type=='igraph' ) return(igraph::V(tree)[leaf_indices])
        if( type=='index' ) return(leaf_indices)
    }

    get_segments <- function(hi_tree, binsize_thresh, return_segmentss_tree=FALSE)
    {
        leaf_widths =  get_leaves(hi_tree, type='igraph')$width_rescaled
        if( binsize_thresh <= min(leaf_widths) )
        {
            warning( 'Your input binsize_thresh has a value: ', binsize_thresh, ', which is smaller than the minimum of the leaf width of hi_tree: ', min(leaf_widths) )
        }
        
        nodes2rm = numeric()
        for( i in 2:igraph::vcount(hi_tree) ) ## root node not take into account
        {
            node = igraph::V(hi_tree)[i]
            width = node$right_rescaled - node$left_rescaled + 1
            if(width <= binsize_thresh)
            {
                parent = igraph::ego(hi_tree, order=1, nodes=node, mode = "in", mindist = 1)[[1]]
                width_parent = parent$right_rescaled - parent$left_rescaled + 1
                if(width_parent <= binsize_thresh) nodes2rm = c(nodes2rm, i)    
            }
        }
        
        trimmed_tree = hi_tree - nodes2rm
        leaves = get_leaves(trimmed_tree, type='igraph_node')
        segmentss = cbind(leaves$left_rescaled, leaves$right_rescaled)
        if(return_segmentss_tree==TRUE) return( trimmed_tree )

        return( segmentss )
    }
    

get_tree_decoration = function( single_res_info, decoration=TRUE, distr, n_parameters, imputation_num=1E2 )
    {
        cA = single_res_info$cA
        # cat(dim(cA), '\n')
        if( !is.null(single_res_info$full_tree) ) tree=single_res_info$full_tree
        if( is.null(single_res_info$full_tree) ) tree = get_tree_v0(single_res_info)
        
        if(length(tree)==1) if( class(tree)=="character" & tree=='bad_tree') return( tree )

        ## more decoration
        if( decoration==FALSE ) return( tree )

        leaf_indices = get_leaves( tree, type='index' )
        igraph::V(tree)$mean_diff = 0
        zero_ratios = numeric()
        for(i in setdiff(1:igraph::vcount(tree), leaf_indices))
        {
            node = igraph::V(tree)[i]
            # igraph::V(tree)[i]$L = L[node$left, node$right]

            # A_union = cA[node$left:node$right, node$left:node$right]
            # igraph::V(tree)[i]$L_union = get_prob_nb( A_union[upper.tri(A_union, diag=TRUE)] )
            
            # the p-value of: H0: unioned model; H1: hierarchical model
            twins = igraph::ego(tree, 1, node, mode='out', mindist=1)[[1]]          
            mid = sort(c(as.numeric(twins[1]$left), as.numeric(twins[1]$right), as.numeric(twins[2]$left), as.numeric(twins[2]$right)))[2]
            
            # if(distr=='nb') test_info = p_likelihood_ratio_nb( cA, head=node$left, mid=mid, tail=node$right )
            # if(distr=='norm') test_info = p_likelihood_ratio_norm( cA, head=node$left, mid=mid, tail=node$right )
            # if(distr=='gamma') 
            # {
                # test_info = p_likelihood_ratio_gamma( cA, head=node$left, mid=mid, tail=node$right, n_parameters, imputation=FALSE )
                # test_info_imputation = p_likelihood_ratio_gamma( cA, head=node$left, mid=mid, tail=node$right, n_parameters, imputation=TRUE )
                # igraph::V(tree)[i]$imp_p = test_info_imputation$p
                # igraph::V(tree)[i]$imp_Lambda = test_info_imputation$Lambda
            # }
            
            if(distr=='lnorm')
            {
                if(!is.null(imputation_num))
                {
                    test_info_imp  = p_likelihood_ratio_lnorm( cA, head=as.numeric(node$left), mid=mid, tail=as.numeric(node$right), n_parameters=n_parameters, imputation= TRUE, imputation_num=imputation_num )       
                    igraph::V(tree)[i]$imp_p = test_info_imp$p
                    igraph::V(tree)[i]$mean_diff = test_info_imp$mean_diff                  
                }

                test_info_nimp = p_likelihood_ratio_lnorm( cA, head=as.numeric(node$left), mid=mid, tail=as.numeric(node$right), n_parameters=n_parameters, imputation=FALSE, imputation_num=imputation_num )       

                igraph::V(tree)[i]$nimp_p = test_info_nimp$p
                igraph::V(tree)[i]$mean_diff = test_info_imp$mean_diff
                # zero_ratios = c(zero_ratios, sum(cA[node$left:node$right, node$left:node$right]==0) / (length(node$left:node$right))^2)
            }

            if(distr=='wilcox')
            {
                if(i==1) test_info = p_wilcox_test( is_CD=TRUE, cA, head=as.numeric(node$left), mid=mid, tail=as.numeric(node$right), alternative='less' )      
                if(i!=1) test_info = p_wilcox_test( is_CD=FALSE, cA, head=as.numeric(node$left), mid=mid, tail=as.numeric(node$right), alternative='less' )     

                # if(i==1) test_info = p_wilcox_test_nested( is_CD=TRUE, cA, head=as.numeric(node$left), mid=mid, tail=as.numeric(node$right), alternative='less' )     
                # if(i!=1) test_info = p_wilcox_test_nested( is_CD=FALSE, cA, head=as.numeric(node$left), mid=mid, tail=as.numeric(node$right), alternative='less' )        


                igraph::V(tree)[i]$nimp_p = igraph::V(tree)[i]$imp_p = igraph::V(tree)[i]$wilcox_p = test_info$p
                igraph::V(tree)[i]$mean_diff = test_info$mean_diff
                # zero_ratios = c(zero_ratios, sum(cA[node$left:node$right, node$left:node$right]==0) / (length(node$left:node$right))^2)
            }
            
            test_info_mean_diff = lognormal_mean_test( cA, head=as.numeric(node$left), mid=mid, tail=as.numeric(node$right) )       

            igraph::V(tree)[i]$aa_p = test_info_mean_diff$p_Aa
            igraph::V(tree)[i]$bb_p = test_info_mean_diff$p_Ab
        }
        
        igraph::V(tree)$L_diff = 0
        # igraph::V(tree)$L_diff = igraph::V(tree)$L - igraph::V(tree)$L_union
        if(!is.null(imputation_num)) igraph::V(tree)[leaf_indices]$imp_p = 1
        igraph::V(tree)[leaf_indices]$nimp_p = igraph::V(tree)[leaf_indices]$wilcox_p = igraph::V(tree)[leaf_indices]$imp_p = 1     
        igraph::V(tree)[leaf_indices]$aa_p = 1
        igraph::V(tree)[leaf_indices]$bb_p = 1
        return(tree)
    }

  
    add_boundary_binsignal_to_decrated_branches <- function( decorated_branches, bin_signals_5_10 )
    {
        tree_widths = sapply(decorated_branches, function(v) igraph::V(v)[1]$width)
        if(sum(tree_widths) != length( bin_signals_5_10 )) stop('Check add_boundary_binsignal_to_decrated_branches')
        for(k in 1:length(decorated_branches))
        {
            tree = decorated_branches[[k]]
            leaf_indices = get_leaves( tree, type='index' )

            for(i in setdiff(1:igraph::vcount(tree), leaf_indices))
            {
                node = igraph::V(tree)[i]
                # igraph::V(tree)[i]$L = L[node$left, node$right]

                # A_union = cA[node$left:node$right, node$left:node$right]
                # igraph::V(tree)[i]$L_union = get_prob_nb( A_union[upper.tri(A_union, diag=TRUE)] )
                
                # the p-value of: H0: unioned model; H1: hierarchical model
                twins = igraph::ego(tree, 1, node, mode='out', mindist=1)[[1]]          
                mid = sort(c(as.numeric(twins[1]$left), as.numeric(twins[1]$right), as.numeric(twins[2]$left), as.numeric(twins[2]$right)))[2]
                
                absolute_mid = mid + sum(tree_widths[(1:k)-1])
                igraph::V(tree)[i]$binsignal = bin_signals_5_10[absolute_mid]
            }
            igraph::V(tree)[leaf_indices]$binsignal = 0
            igraph::V(tree)[which(is.na(igraph::V(tree)$binsignal))]$binsignal = 0
            decorated_branches[[k]] = tree
        }
        return(decorated_branches)
    }
    
    get_tree_v0 <- function( single_res_info )
    {
        ancestors = single_res_info$ancestors
        L = single_res_info$L
        N = dim(ancestors)[1]
        current_node = ancestors[1, N]

        recursive <- function(current_node)
        {
            # cat(current_node, '\n')
            seqs = as.numeric(strsplit(current_node,'-')[[1]])
            # cat( current_node, '\n' )

            left_node = ancestors[seqs[1],seqs[2]]
            right_node = ancestors[seqs[3],seqs[4]]

            flag_left = (seqs[1]!=seqs[2]) & (left_node!="")
            flag_right = (seqs[3]!=seqs[4]) & (right_node!="")

            if(is.na(flag_left)) return( 'bad_tree' ) ## added: 30-04-2020
            
            left_leaf = paste( seqs[1],seqs[2], sep='-' )
            right_leaf = paste( seqs[3],seqs[4], sep='-' )
            
            if( flag_left ) left_tree = recursive(left_node)
            if( !flag_left ) left_tree = left_leaf

            if( flag_right ) right_tree = recursive(right_node)
            if( !flag_right ) right_tree = right_leaf

            if( !flag_left ) left_node = left_leaf
            if( !flag_right ) right_node = right_leaf
            
            tree_raw = igraph::graph.empty() + current_node + left_tree + right_tree
            tree = igraph::add_edges(tree_raw, c(current_node, left_node, current_node, right_node))

            return( tree )
        }
        
        tree = recursive(current_node)
        # stop('I stop here')
        
        # if(tree=='bad_tree') return('bad_tree') ## added: 30-04-2020
        if(class(tree)!='igraph') return('bad_tree')
 
        ## this part tries to decorate the tree by adding various node attributes
        
        igraph::V(tree)$left = sapply( igraph::V(tree)$name, function(v) {tmp=strsplit(v,'-')[[1]]; return(as.numeric(head(tmp,1)))} )
        igraph::V(tree)$right = sapply( igraph::V(tree)$name, function(v) {tmp=strsplit(v,'-')[[1]]; return(as.numeric(tail(tmp,1)))} )
        igraph::V(tree)$name = sapply( igraph::V(tree)$name, function(v) {tmp=strsplit(v,'-')[[1]]; paste( '(', head(tmp,1), ',', tail(tmp,1), ')', sep='' )} )
        igraph::V(tree)$width = igraph::V(tree)$right - igraph::V(tree)$left + 1
        
        return(tree)
    }


    is_binary_tree <-function(tree)
    {
        leaves = get_leaves( tree, 'index' )
        not_leaves = setdiff(1:igraph::vcount(tree), leaves)
        degrees = igraph::degree(tree, v = not_leaves, mode = 'out')
        if(!igraph::is.connected( tree )) return(FALSE)
        if(all(degrees==2)) return(TRUE)
        return(FALSE)
    }
    
    join_left_or_right <- function(tree)
    {
        leaves = get_leaves( tree, 'igraph' )

        left_or_right = numeric()
        left_or_right[1] = 0
        left_or_right[length(leaves) - 1] = 0

        for( i in 2:(length(leaves) - 1) )
        {
            leave = leaves[i]
            parent = igraph::ego(tree, nodes=leave, order=1, mindist=1, mode='in')[[1]]
            if( parent$left < leave$left ) left_or_right[i] = -1 ## left
            if( parent$right > leave$right ) left_or_right[i] = 1 ## right      
        }
        return( left_or_right )
    }

    
    long_slices <- function(left_or_right, thresh=3)
    {
        dup_lens = rle( left_or_right )$lengths
        dup_values = rle( left_or_right )$values
        long_dups = which( dup_lens >= thresh )
        long_dup_lens = dup_lens[long_dups]
        long_dup_start_pos = cumsum( dup_lens )[long_dups-1] + 1 # 1: shift to the right by one
        long_dup_direction = left_or_right[ long_dup_start_pos ]
        indices_of_slice = as.vector(unlist(mapply(function(u,v) {seq(from=u, length=v, by=1)}, long_dup_start_pos, long_dup_lens)))
        dup_info = list( long_slices_lens=long_dup_lens, long_slices_start_pos=long_dup_start_pos, long_dup_direction=long_dup_direction, indices_of_slice=indices_of_slice )
        return(dup_info)
    }
    
    
    n_cells2compute <- function( A, dist, min_n_bins=1 )
    {
        N = nrow(A)
        if(min_n_bins==1)
        {
            count = dist*(1+dist)*(2+dist)*(N-dist)/6
            return(count)
        }
        
        count = (dist - 2*min_n_bins + 2)*(dist^2 + 2*dist*min_n_bins + dist - 2*(min_n_bins - 2)*min_n_bins)*(N - dist) / 6
        return(count)
    }
    
    
    
    plot_tree <- function(tree, seg_col='blue', which_part='upper', indices_of_slice=NULL, ...)
    {
        root_node = igraph::V(tree)[1]
        leaves = igraph::V(tree)[which(igraph::degree( tree, mode='out' ) == 0)]
        plot_inner <- function( node ){
            left = node$left
            right = node$right
            if(right - left == 1)
            {
                return()
            }
            
            x0 = left - 0.5
            x1 = right + 0.5
            if(which_part=='upper')
            {
                segments(x0, x0, x0, x1, col=seg_col, ...)
                segments(x0, x1, x1, x1, col=seg_col, ...)      
            }

            if(which_part=='lower')
            {
                segments(x0, x0, x1, x0, col='black', ...)
                segments(x1, x1, x1, x0, col='black', ...)          
            }

            if(which_part=='both')
            {
                segments(x0, x0, x0, x1, col=seg_col, ...)
                segments(x0, x1, x1, x1, col=seg_col, ...)      
                segments(x0, x0, x1, x0, col='black', ...)
                segments(x1, x1, x1, x0, col='black', ...)  
            }
            
            if( !is.null(indices_of_slice) )
            {
                if(left %in% indices_of_slice)
                segments(x0, x0, x0, x1, col='yellow', ...)

                if(right %in% indices_of_slice)
                segments(x0, x1, x1, x1, col='black', ...)  
            }
            
            if( node %in% leaves ) return()
            twins = igraph::ego(tree, node=node, order=1, mindist=1, mode='out')[[1]]
            left_node = twins[1]
            plot_inner(left_node)
            if( length(twins) > 1 )
            {
                # cat('hello', '\n')
                right_node = twins[2]
                plot_inner(right_node)
            }
        }
        plot_inner(root_node)
        return()
    }

    remove_blank_cols <- function( mat, sparse=FALSE, row_or_col='both', ratio=0.05 )
    {
        if(sparse==TRUE)
        {
            if(row_or_col=='row')
            {
                non_zero_indices = (Matrix::rowSums(mat != 0) / ncol(mat) >  ratio)
                new_mat = mat[non_zero_indices, ]
                return(new_mat)
            }

            if(row_or_col=='col')
            {
                non_zero_indices = (Matrix::colSums(mat != 0) / nrow(mat) >  ratio)
                new_mat = mat[, non_zero_indices]
                return(new_mat)
            }

            if( ratio > 0 ) ## remove col/rows with non-zero number smaller than 5% percentile
            {
                positive_num = unname(Matrix::rowSums(mat != 0))
                positive_num_thresh = quantile(positive_num[positive_num > 0], ratio)
                non_zero_indices = (Matrix::rowSums(mat != 0) >  positive_num_thresh)
                new_mat = mat[non_zero_indices, non_zero_indices]
                return(new_mat)
            }


            non_zero_indices = (Matrix::rowSums(mat != 0) / ncol(mat) >  ratio)
            new_mat = mat[non_zero_indices, non_zero_indices]
            return(new_mat)
        }

        col_sum = apply( mat, 2, sum )
        blank_col_indices = which( col_sum==0 )
        new_mat = mat[ -blank_col_indices, -blank_col_indices ]
        return(new_mat)
    }

    
    tree_germination <- function(tree)
    {
        leaves = get_leaves( tree )
        for( leaf in leaves )
        {
            cat(leaf, '\n')
            if((igraph::vcount(tree) - ecount(tree)) !=1) break
            bins = as.character(igraph::V(tree)[leaf]$left:igraph::V(tree)[leaf]$right)
            if(length(bins) == 1) stop('Some leave node contains only one bin, please check')
            if(length(bins) == 2) ## this is an exception. added on 11/06/2018
            {
                nodes_supp = as.character(bins)
                edges2add = c(leaf, bins[1], leaf, bins[2])
                tree = tree %>% igraph::add_vertices(length(nodes_supp), name=nodes_supp) %>% igraph::add_edges(edges2add)
                next
            }
            
            joints_supp = paste('j', bins[1:(length(bins) - 2)], sep='')
            joints = c(leaf, joints_supp, bins[length(bins)])
            edges2add = character()
            for(i in 1:(length(joints)-1))
            {
                edges2add = c(edges2add, c( joints[i], bins[i], joints[i], joints[i+1] ))
            }
            nodes_supp = c(joints_supp, bins)
            tree = tree %>% igraph::add_vertices(length(nodes_supp), name=nodes_supp) %>% igraph::add_edges(edges2add)
        }
        return( tree )
    }

    
    trim_tree_adaptive <- function( tree, L_diff_thresh=-Inf, max_imp_p=Inf, max_nimp_p=Inf, width_thresh=-Inf, boundary_signal_thresh=Inf, peak_thresh=-Inf )
    {
        if(igraph::vcount(tree) ==1) return(tree) ## so this can be used!

        width_thresh = min(width_thresh, igraph::V(tree)[1]$width-1) ## cannot merge if tree is already very small
        
        nodes2rm_Ldiff = igraph::V(tree)[which( igraph::V(tree)$L_diff <= L_diff_thresh )]$name
        
        
        nodes2rm_imp_p = igraph::V(tree)[which( igraph::V(tree)$imp_p >= max_imp_p )]$name
        nodes2rm_nimp_p = igraph::V(tree)[which( igraph::V(tree)$nimp_p >= max_nimp_p )]$name
        nodes2rm_boundary_signal = igraph::V(tree)[which( igraph::V(tree)$binsignal <= boundary_signal_thresh )]$name
        nodes2rm_not_peak = igraph::V(tree)[which( igraph::V(tree)$is_peak <= peak_thresh )]$name
        
        nodes2rm = unique( c( nodes2rm_Ldiff, nodes2rm_imp_p, nodes2rm_nimp_p, nodes2rm_boundary_signal, nodes2rm_not_peak ) )
        nodes_not_2rm = setdiff( igraph::V(tree)$name, nodes2rm )

        flags = sapply( nodes2rm, function(v) { children = names(unlist( igraph::ego( tree, nodes = v, order=igraph::diameter(tree), mindist=1, mode='out' ))); all(children %in% nodes2rm) } )
        nodes2rm_final = names(unlist( igraph::ego(tree, nodes = nodes2rm[which(flags==TRUE)], order=igraph::diameter(tree), mindist=1, mode='out' ) ) )
        
        if(!is.null(nodes2rm_final)) tree = tree - nodes2rm_final
        if(is.null(nodes2rm_final)) cat('No nodes are removed at the given p_thresh\n')
    
        if( width_thresh==-Inf ) return( tree )
        
        tree = trim_tree_adaptive_width_thresh(tree, width_thresh)
        return( tree )
    }
    
    trim_tree_adaptive_width_thresh <- function( tree, width_thresh )
    {
        ## first remove all nodes of a parent if both are < width_thresh
        nodes2rm_width = igraph::V(tree)[which( igraph::V(tree)$width <= width_thresh )]$name
        parentOfnodes2rm_width = names(unlist(igraph::ego(tree, node=nodes2rm_width, order=1, mindist=1, mode='in')))
        ## if two nodes have the same parent, remove them
        index2rm = c(which(duplicated(parentOfnodes2rm_width)), which(duplicated(fromLast=TRUE, parentOfnodes2rm_width)))
        nodes2rm = nodes2rm_width[index2rm]
        tree = tree - nodes2rm
        if(!is_binary_tree(tree)) stop('Error in trim_tree_adaptive_width_thresh')

        ## merge TADs that are too small    
        nodes2rm_width = igraph::V(tree)[which( igraph::V(tree)$width <= width_thresh )]$name

        di = igraph::diameter(tree)
        while( length(nodes2rm_width) > 0 )
        {
            node = nodes2rm_width[1]
            parent = igraph::ego( tree, nodes = node, order=1, mindist=1, mode='in' )[[1]]
            nodesOfSameParent = igraph::ego( tree, nodes = parent, order=di, mindist=1, mode='out' )[[1]]
            left_siblings = nodesOfSameParent[which(nodesOfSameParent$right < igraph::V(tree)[node]$left)]
            right_siblings = nodesOfSameParent[which(nodesOfSameParent$left > igraph::V(tree)[node]$right)]
            
            left_siblings2change_right = left_siblings[left_siblings$right == (igraph::V(tree)[node]$left-1)]
            right_siblings2change_left = right_siblings[right_siblings$left == (igraph::V(tree)[node]$right+1)]
            
            igraph::V(tree)[ left_siblings2change_right ]$right = igraph::V(tree)[node]$right
            igraph::V(tree)[ right_siblings2change_left ]$left = igraph::V(tree)[node]$left
            igraph::V(tree)$width = igraph::V(tree)$right - igraph::V(tree)$left + 1
            tree = tree - node
            nodes2rm_width = igraph::V(tree)[which( igraph::V(tree)$width <= width_thresh )]$name
        }
        igraph::V(tree)$name = paste('(',igraph::V(tree)$left, ',', igraph::V(tree)$right, ')', sep='')
        return(tree)
    }
    
    ## should be different from trim_tree_adaptive_width_thresh because binsignal is not monotonic
    trim_tree_adaptive_binsig_thresh <- function( tree, boundary_signal_thresh )
    {

        ## merge TADs that are too small    
        parent_of_nodes2rm = igraph::V(tree)[which( igraph::V(tree)$binsignal <= boundary_signal_thresh )]$name
        leaves = get_leaves( tree )
        nodes2rm = intersect(leaves,  names(unlist(igraph::ego( tree, nodes = parent_of_nodes2rm, order=1, mindist=1, mode='out' ))))

        while( length(nodes2rm) > 0 )
        {
            di = igraph::diameter(tree)
            node = nodes2rm[1]
            parent = igraph::ego( tree, nodes = node, order=1, mindist=1, mode='in' )[[1]]
            sibling = setdiff(igraph::ego( tree, nodes = parent, order=1, mindist=1, mode='out' )[[1]]$name, node)
            twins_of_sibling = igraph::ego( tree, nodes = sibling, order=1, mindist=1, mode='out' )[[1]]$name
            nodesOfSameParent = igraph::ego( tree, nodes = parent, order=di, mindist=1, mode='out' )[[1]]
            left_siblings = nodesOfSameParent[which(nodesOfSameParent$right < igraph::V(tree)[node]$left)]
            right_siblings = nodesOfSameParent[which(nodesOfSameParent$left > igraph::V(tree)[node]$right)]
            
            left_siblings2change_right = left_siblings[left_siblings$right == (igraph::V(tree)[node]$left-1)]
            right_siblings2change_left = right_siblings[right_siblings$left == (igraph::V(tree)[node]$right+1)]
            
            igraph::V(tree)[ left_siblings2change_right ]$right = igraph::V(tree)[node]$right
            igraph::V(tree)[ right_siblings2change_left ]$left = igraph::V(tree)[node]$left
            
            
            if(length(twins_of_sibling) > 0) tree = igraph::add_edges(tree, c( parent$name, twins_of_sibling[1], parent$name, twins_of_sibling[2]) )
            tree = tree - node - sibling
            
            igraph::V(tree)$width = igraph::V(tree)$right - igraph::V(tree)$left + 1
            
            parent_of_nodes2rm = igraph::V(tree)[which( igraph::V(tree)$binsignal <= boundary_signal_thresh )]$name
            leaves = get_leaves( tree )
            nodes2rm = intersect(leaves,  names(unlist(igraph::ego( tree, nodes = parent_of_nodes2rm, order=1, mindist=1, mode='out' ))))
        }
        igraph::V(tree)$name = paste('(',igraph::V(tree)$left, ',', igraph::V(tree)$right, ')', sep='')
        return(tree)
    }

    visualize_left_or_right <- function(left_or_right, ...)
    {
        len = length( left_or_right )
        plot(1, type="n", xlab="", ylab="", xlim=c(0, len), ylim=c(-2, 2))
        for( i in 1:len )
        {
            if( left_or_right[i]==1 ) segments(i-0.1, 1, i+1, 1, col='red', ...)
            if( left_or_right[i]==-1 ) segments(i-0.1, -1, i+1, -1, col='blue', ...)
        }
    }

    xenocraft <- function( trunk, branches )
    {
        
        xenocraft_nodes = sapply(branches, function(v) igraph::V(v)[1]$name)
        if(!all(xenocraft_nodes %in% igraph::V(trunk)$name)) stop("Check xenocraft function")
        
        children_xenocraft_nodes = unique(unlist(sapply(igraph::ego(trunk, order=igraph::diameter(trunk), node=xenocraft_nodes, mode='out', mindist=1), function(v) v$name)))
        prunned_trunck = trunk - children_xenocraft_nodes
        full_tree = Reduce( union, c(list(prunned_trunck), branches) )
        att2delete = setdiff(vertex_attr_names(full_tree), 'name')
        for( att in att2delete ) full_tree = delete_vertex_attr(full_tree, att) 
        return( full_tree )
    }

    fast_cor <- function(mat)
    {
        res = 1/( NROW( mat ) -1)*crossprod ( scale( mat , TRUE , TRUE ) )
        # scaled_mat = scale( mat , TRUE , TRUE )
        # res = 1/( NROW( mat ) -1)*matrix_multiplication_sym_cpp( scaled_mat )

        return(res)
    }

    fast_cor_cor <- function(mat, k)
    {
        scaled_mat = scale(mat , TRUE , TRUE )
        coeff = 1/( NROW( mat ) -1) ## there is no need to multiply this coeff when computing Pearson coeff

        res_begin = coeff*crossprod( scaled_mat[, 1:k], scaled_mat )

        for(j in (k+1):nrow(mat))
        {
            res_slice = res_begin[-1,]
            res_next = coeff*crossprod( scaled_mat[, j], scaled_mat )
            res_slice = rbind(res_slice, res_next)
            new_values = cor(t(res_slice), t(res_next))
            cat(j, '\n')        
        }
    }


    get_bin_singals_CHiP <- function(chr, hc_ordered, res, ks, ChiP_NAME, CELL_LINE, ChiP_data_already_loaded=FALSE)
    {

        #################
        n_bins_of_CD = sapply(res$initial_clusters, length)
        pos_start_end = lapply(res$initial_clusters, function(v)
        {
            bins_ori = as.numeric(res$bin_names[v])
            from_id = min(bins_ori)
            to_id = max(bins_ori)
            start_pos = (from_id-1)*bin_size + 1
            end_pos = to_id*bin_size
            return(c(start_pos, end_pos))
        })

        pos_start_end = do.call(rbind, pos_start_end)
        #################

        ordered_bin_signals_ALL_k_list = lapply(ks, function(k){

            hc_k_labels = get_cluser_levels(hc_ordered, k_clusters=k, balanced_4_clusters=FALSE)$cluster_labels

            if(is.null(ChiP_NAME))
            {
                CD_index = labels(hc_ordered)
                pos_start_end = pos_start_end[match( CD_index, rownames(pos_start_end) ), ]
                ordered_bin_signals_ALL = data.frame(CD_index=CD_index, n_bins=n_bins_of_CD[CD_index], pos_start=pos_start_end[,1], pos_end=pos_start_end[,2],  compartment=hc_k_labels[CD_index])
                return(ordered_bin_signals_ALL)
            }

            cluster_signal = get_names_by_H3k4me1(chr, res, ChiP_NAME, kmeans_cluster_assigment=hc_k_labels, CELL_LINE=CELL_LINE, ChiP_data_already_loaded=ChiP_data_already_loaded)
            # cluster_signal = rbind(cluster_signal, sorted_coverage)
            ordered_bin_signals_ALL = apply(cluster_signal, 1, function(v)
            {
                bin_signals = v[hc_k_labels]
                names(bin_signals) = names(hc_k_labels)
                ordered_bin_signals = bin_signals[labels(hc_ordered)] ## order by dendro topology
                return(ordered_bin_signals)
            })

            ordered_bin_signals_ALL = as.data.frame(ordered_bin_signals_ALL)
            # CD_index = rownames(ordered_bin_signals_ALL)
            CD_index = labels(hc_ordered)

            pos_start_end = pos_start_end[match( CD_index, rownames(pos_start_end) ), ]

            ordered_bin_signals_ALL = data.frame(CD_index=CD_index, n_bins=n_bins_of_CD[CD_index], pos_start=pos_start_end[,1], pos_end=pos_start_end[,2],  compartment=hc_k_labels[CD_index], ordered_bin_signals_ALL)
            return(ordered_bin_signals_ALL)
            })
        names(ordered_bin_signals_ALL_k_list) = paste0('clusters_', ks)
        names(ordered_bin_signals_ALL_k_list)[which(names(ordered_bin_signals_ALL_k_list)=="clusters_Inf")] = 'compartment_domains'
        return(ordered_bin_signals_ALL_k_list)
    }



        densplot = function(x,y,points = FALSE, pch=19, cex=1, xlim=c(min(x),max(x)), ylim=c(min(y),max(y)), ...){
        df = data.frame(x,y)
        d = densCols(x,y, colramp=colorRampPalette(c("black", "white")))
        df$dens = col2rgb(d)[1,] + 1L
        
        cols = colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
        
        df$col = cols[df$dens]
        df=df[order(df$dens),]
        if(points)
                points(df$x,df$y, pch=pch, col=df$col, cex=cex, ...)
        else
                plot(df$x,df$y, pch=pch, col=df$col, cex=cex, xlim=xlim, ylim=ylim, ...)
        }



    my_merge = function(...) merge(..., all=TRUE)
    
    
    build_chr_bin_domain_fun <- function( CELL_LINE, chrs, cluster_level, p_thresh, ob_oe, downsratio=NULL, compress_size=10 )
    {
        chrs = as.character(chrs)
        chr_bin_domain_tmp = lapply(chrs, function(chr) get_clusters_bins_xy(CELL_LINE, chr, cluster_level, p_thresh, ob_oe, downsratio=downsratio, compress_size=compress_size))
        names( chr_bin_domain_tmp ) = chrs ## can use mapply
        chr_bin_domain = lapply(chrs, function(v) 
        {
            chr_bin_domain_ind = data.frame( chr=paste0('chr', v), bin_index=as.numeric(unlist(chr_bin_domain_tmp[[v]])), intra_domain=rep(names(chr_bin_domain_tmp[[v]]), sapply(chr_bin_domain_tmp[[v]], length)) )
            chr_bin_domain_ind = chr_bin_domain_ind[order(chr_bin_domain_ind$bin_index), ]
            return(chr_bin_domain_ind)
        })
        names(chr_bin_domain) = chrs
        res = do.call(rbind, chr_bin_domain)
        rownames(res) = NULL
        return( res )
    }

    build_chr_bin_domain_fun_direct <- function( chr, initial_clusters, cluster_vec, bin_names ) ## directly from R workingspace instead of loading
    {
        chrs = as.character(chr) ## for historic reason, chr -- chrs
        chr_bin_domain_tmp = lapply(chrs, function(chr) get_cluster_bin_names(initial_clusters, cluster_vec, bin_names))
        names( chr_bin_domain_tmp ) = chrs ## can use mapply
        chr_bin_domain = lapply(chrs, function(v) data.frame( chr=paste0('chr', v), bin_index=as.numeric(unlist(chr_bin_domain_tmp[[v]])), intra_domain=rep(names(chr_bin_domain_tmp[[v]]), sapply(chr_bin_domain_tmp[[v]], length)) ))
        names(chr_bin_domain) = chrs
        return( do.call(rbind, chr_bin_domain) )
    }

    get_clusters_bins_xy <- function(CELL_LINE, chr, cluster_level, p_thresh, ob_oe='oe', downsratio, compress_size, sort=TRUE)
    {
        if(ob_oe=='oe') sub_folder = paste0('./', CELL_LINE, '/oe_chr_', chr, '_', bin_size/1E3, 'kb_', compress_size, 'to1_', p_thresh)
        if(ob_oe=='ob') sub_folder = paste0('./', CELL_LINE, '/ob_chr_', chr, '_', bin_size/1E3, 'kb_', compress_size, 'to1_', p_thresh)

        if(!is.null(downsratio)) compartments_Rdata_file = paste0(sub_folder, '/chr', chr, '_compartments_atanh_log_AB', ws, 'downsratio_', downsratio, '.Rdata')

        if(is.null(downsratio)) compartments_Rdata_file = paste0(sub_folder, '/chr', chr, '_compartments_atanh_log_AB3_3.Rdata')
        
        load(compartments_Rdata_file)
        clusters_bins = get_cluster_bin_names(sort=sort, res$initial_clusters, res$clusters[, cluster_level], res$bin_names)
        rm( res )
        return(clusters_bins)
    }


    get_cluster_bin_names <- function(initial_clusters, cluster_vec, bin_names, sort=TRUE)
    {
        if(sort==TRUE) cluster_indices = sort(unique(cluster_vec), decreasing=TRUE)
        if(sort==FALSE) cluster_indices = unique(cluster_vec)

        cluster_bins = lapply(cluster_indices, function(v)
        {
            indices = which(cluster_vec==v)
            bin_names[unlist(initial_clusters[indices])]
        })
        names(cluster_bins) = cluster_indices
        return(cluster_bins)
    }

    # ave_cor <- function(mat, seg_len) ## seg_len is the length of a segment
    # {
    #   d = 1:nrow(mat)
    #   seq_index = split(d, ceiling(seq_along(d)/10))
    #   tmp = simplify2array( lapply(seq_index[1:22], function(v) fast_cor( mat[v, ] )/10*length(v)) )
    #   ave_cor_val = rowMeans(tmp, dims = 2)
    # }

    # if( kmeans_cluster_assigment ): compute the cluster_level signals
   

    reorder_dendro <- function(hc_object, named_weights, return_g=FALSE, aggregateFun=mean)
    {
        get_children <- function(g, root_node)
        {
            leaves = igraph::V(g)[igraph::degree(g)==1]$name
            children = intersect(leaves, igraph::ego(g, mode='out', root_node, order=igraph::diameter(g))[[1]]$name)
            return(children)
        }

        # hc_dendro = as.dendrogram(hc_object)
        # cat('hello I am here\n')
        # print(class(hc_object))
        g = igraph::as.igraph(ape::as.phylo(hc_object))
        # stop('I am here')
        leaves = igraph::V(g)[igraph::degree(g)==1]$name

        igraph::V(g)$weight = sapply(1:igraph::vcount(g), function(v) aggregateFun(named_weights[get_children(g, v)]))

        if(return_g==TRUE) return(g)

        swap_branches <- function(g, root_node)
        {
            twins = igraph::ego(g, mode='out', root_node, order=1, mindist=1)[[1]]$name
            twins_weight = igraph::V(g)[twins]$weight

            if(twins_weight[1] > twins_weight[2])
            {
                children_of_twin_A = get_children(g, twins[1])
                children_of_twin_B = get_children(g, twins[2])
                leaves = swap_names( leaves, children_of_twin_A, children_of_twin_B )
            }
            return(leaves)
        }

        swap_names <- function( leaves, names2swap_A, names2swap_B )
        {
            names2swap_indices = match(c(names2swap_A, names2swap_B), leaves)
            leaves[ names2swap_indices ] = c(names2swap_B, names2swap_A)
            return(leaves)
        }

        for( root_node in setdiff(igraph::V(g)$name, leaves) )
        {
            leaves = swap_branches(g, root_node)
            # cat(leaves, '\n')
        }
            
        return(leaves)
    }


    get_cluser_assignment = function(hc, k_clusters, leaves_hclust_pc)
    {
        clusters_raw = cutree(hc, k_clusters)
        clusters = tapply(as.numeric(names(clusters_raw)), clusters_raw, function(v) list(v))

        od = rank(sapply( clusters, function(v) sort(match(v, as.numeric(leaves_hclust_pc)))[1]))

        cluster_assignment = numeric(sum(sapply(clusters, length)))
        for(j in 1:length(clusters)) cluster_assignment[clusters[[j]]] = od[j]
        return( cluster_assignment )
    }



    get_cluster_boudaries <- function(hc, k_clusters, named_weights)
    {
        len = length(labels(hc))
        clusters_raw = cutree(hc, k_clusters)
        clusters = tapply(as.numeric(names(clusters_raw)), clusters_raw, function(v) list(v))
        clusters_pc_rank = lapply(clusters, function(v) unname(sort(rank(named_weights)[v])))
        bondaries = sapply(clusters_pc_rank, function(v) tail(v,1))
        return(setdiff(bondaries, c(1, len)))   
    }

    # get_cluser_vector <- function(hc, k_clusters, named_weights) ## get cluster assignment of 1,2,3,4,5... (A1, A2, ...)
    # {
    #   clusters_raw = cutree(hc, k_clusters)
    #     clusters = tapply(as.numeric(names(clusters_raw)), clusters_raw, function(v) list(v))
    #     clusters_pc_rank = unname(rank(sapply(clusters, function(v) unname(sort(rank(named_weights)[v]))[1])))
    #     cluster_vector = numeric()
    #     for( i in 1:length(clusters) ) cluster_vector[ clusters[[i]] ] = clusters_pc_rank[i]
    #     return( cluster_vector )  
    # }

    ## if return binary tree, A1, A2, B1, B2 are forced to be returned
    get_cluser_levels <- function(hc_ordered, k_clusters, balanced_4_clusters=FALSE) ## get the detailed A1, A2, ..., B1, B2...
    {
        assign_twins_name <- function(graph, node)
        {
            twins = igraph::ego(graph, node, mode='out', order=1, mindist=1)[[1]]$name
            igraph::V(graph)[twins]$level_name = paste0(igraph::V(graph)[node]$level_name, '.', c(2,1))
            if(node==igraph::V(graph)[1]$name) igraph::V(graph)[twins]$level_name = c('B', 'A')
            return(graph)
        }

        if(k_clusters==Inf) k_clusters = length(labels(as.dendrogram(hc_ordered)))
        #################################################
        # cat('haha I am here\n')
        # print(class(hc_ordered))
        graph = igraph::as.igraph(ape::as.phylo(hc_ordered))
        leave_names = get_leaves(graph)

        bfs_names = igraph::bfs(graph, 1)$order$name
        dfs_names = igraph::dfs(graph, 1)$order$name

        igraph::V(graph)[1]$level_name = ''
        for( node in bfs_names )
        {
            graph = assign_twins_name(graph, node)
            # if( !any(is.na(igraph::V(graph)[common_father_name]$level_name)) ) break ## when all common_father_name have level_name
            if( !any(is.na(igraph::V(graph)$level_name)) ) break ## when all common_father_name have level_name
        }
        
        if( balanced_4_clusters==TRUE ) ##A1, A2, B1, B2
        {
            branch_root_name = c('A.1', 'A.2', 'B.1', 'B.2')
            branch_root = match(branch_root_name, igraph::V(graph)$level_name)

            children = igraph::ego(graph, order = igraph::diameter(graph), nodes = branch_root, mode = 'out', mindist = 0) ## mindist==0, when itself is a branch
            tmp = lapply( children, function(v) intersect(leave_names, v$name) )
            cluster_labels = rep(branch_root_name, sapply(tmp, length))
            names(cluster_labels) = unlist( tmp )
            cluster_labels = cluster_labels[ as.character(1:length(cluster_labels)) ]
            return(cluster_labels)
        }

        #################################################
        clusters_raw = cutree(hc_ordered, k_clusters)[labels(hc_ordered)] ## labels are ordered according to pc 
        clusters = tapply(as.numeric(names(clusters_raw)), clusters_raw, function(v) list(v))[unique(clusters_raw)]
        names(clusters) = 1:k_clusters # named by pc order 
        
        ## get the vector. Named from 1 to 5 from left leaves to right leaves
        cluster_vector = rep( 1:k_clusters, sapply(clusters, length) )
        names(cluster_vector) = labels(hc_ordered)
        cluster_vector = cluster_vector[ as.character( sort(as.numeric(names(cluster_vector))) ) ]

        #################################################

        ## which node is the common father of all nodes in a cluster
        common_father_index = sapply( clusters, function(u) 
            {
                if(length(u)==1) return( which(u==igraph::V(graph)$name) )
                return(max(which(sapply(igraph::ego(graph, order = igraph::diameter(graph), nodes = igraph::V(graph), mode = 'out', mindist = 1), function(v) all(u %in% v$name) ))))
            })

        common_father_name = igraph::V(graph)[common_father_index]$name





        ## whether common_father_name are ordered
        if(!all(diff(match( common_father_name, dfs_names )) >= 0)) stop('!!!!!!!!check get_all_children in header_funs.R')
        


        level_names = igraph::V(graph)[common_father_name]$level_name ## label names are ordered
        if(any(sort(level_names, decreasing=TRUE)!=level_names)) warning('!!!!!!!!Need to check get_all_children in header_funs.R')

        cluster_labels = rep( level_names, sapply(clusters, length) )
        names(cluster_labels) = labels(hc_ordered)
        cluster_labels = cluster_labels[ as.character( sort(as.numeric(names(cluster_labels))) ) ]

        return( list(cluster_vector=cluster_vector, cluster_labels=cluster_labels))
    }



     
    get_hkmeans_cluser_levels <- function(hk_cluster_centers, PC1)
    {
        hc = hclust(as.dist(my_dist(hk_cluster_centers)), method='com')
        reordered_names = as.character(rank(PC1))
        hc_ordered = dendextend::rotate(hc, reordered_names) 
        hk_clust_labels = get_cluser_levels(hc_ordered, k_clusters=length(PC1))$cluster_labels
        return( hk_clust_labels )
    }   


    
    ## names_A_final is the rownames of the A_final
    ## names_A_final is matched to the names of the starting contact matrix
    
    get_original_tad_indices <- function(names_A_final, TADs, bin_size)
    {
        start_pos = as.numeric(names_A_final[TADs$start_pos])
        end_pos = as.numeric(names_A_final[TADs$end_pos])
        start_pos_ori = (start_pos - 1)*bin_size + 1
        end_pos_ori = end_pos*bin_size
        TADs = data.frame( start_pos_ori=start_pos_ori, end_pos_ori=end_pos_ori )
        return( TADs )
    }
    



    ## This code updates the branch name 
    ## Branches obtained from branches = lapply( res_inner, get_tree_v0 )
    ## The original branch name start from 1
    update_branch_name <- function(branch, root_start)
    {
        igraph::V(branch)$left = igraph::V(branch)$left + root_start - 1
        igraph::V(branch)$right = igraph::V(branch)$right + root_start - 1
        igraph::V(branch)$name = paste('(',igraph::V(branch)$left, ',', igraph::V(branch)$right, ')', sep='')
        return(branch)
    }

