
	
 	## k-means with replicatable seeds
 	my_kmeans = function(iter.max=1E3, nstart=50, ...)
 	{
 		set.seed(1)
 		res = kmeans(iter.max=iter.max, nstart=nstart, ...)
 		return(res)
 	}

 	## This function tries to adjust the height of each split, in order to generate a valid hclust object and with balanced compartments A.1 A.2 B.1 B.2
 	## Clusters with more nodes will get bigger height in case of same height
 	adjust_hs <- function(l_r_h)
 	{
 		hs = sapply(l_r_h, function(v) v$h)
 		all_names = sapply(l_r_h, function(v) paste0(collapse='_', sort(c(v$l, v$r))))
  		r_names = sapply(l_r_h, function(v) paste0(collapse='_', sort(c(v$r))))
		
	 	sizes = sapply(l_r_h, function(v) length(v$l) + length(v$r)) ##

	 	################ This part deals with duplicated heights
	 	hs = hs + sizes*1E-7
	 	################ This part tries to make the top-level left and right branch to have similar height, such that to make balanced A.1, A.2, B.1, B.2 compartments
	 	## Find the index of second branch, whose number of nodes is n_total - n_left: sizes[1] - sizes[2]
	 	l_b = 2 ## left sub-branch
	 	# r_b = which(sizes==(sizes[1] - sizes[2]))[1] ## right sub-branch
	 	r_b = which(r_names[1]==all_names) ## right sub-branch

	 	l_h = hs[l_b]
	 	r_h = hs[r_b]
	 	max_h = max(l_h, r_h) ## the maximum height of the two branches
	 	hs_new = mean(sort(hs, decreasing=TRUE)[2:3]) ## hs_new is the 3rd largest height
		hs[l_b] = ifelse(l_h > r_h, max_h, hs_new)
		hs[r_b] = ifelse(r_h > l_h, max_h, hs_new)

	 	if(any(duplicated(hs))) stop('ERROR: DUPLICATED HEIGHTS exist in bisecting_kmeans')
	 	return( hs )
 	}

 	bisecting_kmeans <- function(data)
 	{
 		dist_mat = as.matrix(stats::dist(data))
 		indices = 1:nrow(data)
 		l_r_h <<- list()

	 	get_h <- function(l_indices, r_indices)
	 	{
			combined_indices = c(l_indices, r_indices)
			idx <- as.matrix(expand.grid(combined_indices, combined_indices))
			max(dist_mat[idx]) ## diameter		
	 	}

	 	get_sub_tree <- function( indices )
	 	{
		 	n_nodes = length(indices)

		 	if(n_nodes==1) ## if only two nodes
		 	{
		 		h = NULL
		 		# tree = list(h=h, leaf=indices)
		 		return()
		 	}

		 	############# if more than two nodes
		 	if(n_nodes==2) cluster=c(1,2) else cluster = my_kmeans(x=data[indices, ], centers=2)$cluster
		 	l_indices = indices[cluster==1]
		 	r_indices = indices[cluster==2]
		 	h = get_h(l_indices, r_indices)
		 	l_r_h <<- c(l_r_h, list(list(l=l_indices, r=r_indices, h=h)))

		 	# cat(h, '\n')
			l_branch = get_sub_tree( l_indices )
			r_branch = get_sub_tree( r_indices )
		 	# tree = list(h=h, l_branch=l_branch, r_branch=r_branch, l_indices=l_indices, r_indices=r_indices)
		 	# return(tree)
	 	}

	 	get_sub_tree(indices)

	 	hs = adjust_hs(l_r_h)

	 	r_hs = rank(hs)
	 	for( i in 1:length(l_r_h) ) {name=r_hs[i]; names(name)=paste0(collapse='_', sort(c(l_r_h[[i]]$l, l_r_h[[i]]$r))); l_r_h[[i]]$name=name}
	 	pos_names = sapply(l_r_h, function(v) v$name)
	 	neg_names = -(1:length(indices)); names(neg_names) = 1:length(indices); all_names = c(pos_names, neg_names)
	 	for( i in 1:length(l_r_h) ) {l_r_h[[i]]$l_name=unname(all_names[paste0(l_r_h[[i]]$l, collapse='_')]); l_r_h[[i]]$r_name=unname(all_names[paste0(l_r_h[[i]]$r, collapse='_')]) }
		
		merge_height = data.frame(l=sapply(l_r_h, function(v) v$l_name), r=sapply(l_r_h, function(v) v$r_name), h=hs)
		merge_height = merge_height[order(merge_height$h), ]
		rownames(merge_height) = NULL
		
		data_tmp = cbind(c(0,0,1,1), c(0,1,1,0))
		hc = hclust(stats::dist(data_tmp), "com")
		hc$merge = as.matrix(unname(merge_height[,1:2]))
		hc$height = merge_height$h
		# hc$order = unname(unlist(res, recursive=TRUE)[grepl('leaf', names(unlist(res, recursive=TRUE)))])
		# hc$order = 1:length(indices)
		hc$labels = 1:length(indices)
		den <- as.dendrogram(hc)
		hc_r <- as.hclust(reorder(den, 1:length(indices)))
		hc_r$method = "complete"
		hc_r$dist.method = "euclidean"
		l_r_h <<- list()
		rm(l_r_h)
		return(hc_r)
 	}
 
