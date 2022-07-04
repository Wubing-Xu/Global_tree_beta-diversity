##A function to calculate beta diversity among a focal cell and its neighboring cells.
# The first two rows in community data should be the coordinates (position).

library(betapart)
library(picante)
library(ape)
library("phytools")
library(spdep)

get_beta <- function(comdata=NA, tree=NA, trait=NA, method="sorensen", root=TRUE, null=FALSE, 
                     beta.multi=FALSE, alpha=FALSE, n=4, d2=300, ncell_resample=NA){
  
  # drop species without distributions in the phylogeny
  if(any(!is.na(tree))){
    tree <- keep.tip(tree, colnames(comdata)[-c(1:2)])
    tree$root.edge <- 0
  }
  
  # drop species without distributions in the traits (rownames are species names)
  if(any(!is.na(trait))){
    trait <- trait[rownames(trait) %in% colnames(comdata)[-c(1:2)], ]
  }

  # global null model: randomize species across all species in the phylogeny and trait data
	if(null == "global"){
		if(any(!is.na(tree)))  tree$tip.label <- sample(tree$tip.label)
		if(any(!is.na(trait))) rownames(trait) <- sample(rownames(trait))
	}
	
	# add a root in the phylogeny
	if(any(!is.na(tree)) & root){
		tree <- ladderize(tree)
		tree.tips <-  tree$tip.label
		root.node <- getMRCA(tree,c(tree.tips[1], tree.tips[length(tree.tips)]))
		tree <- bind.tip(tree, "root", where=root.node)
	}
	
	# get neighboring grid cells of the focal cells
	nbs <- dnearneigh(x=comdata[, 1:2], d1=0, d2=d2, row.names=rownames(comdata), longlat=FALSE)
	nbsmat <- nb2mat(neighbours=nbs, style="B", zero.policy=TRUE)

	# the grid cells that have more than n neighboring grid cells
	rows.beta <- which(rowSums(nbsmat) >= n)
	
	
	# a sub-function to calculate beta for each focal cell
	get_beta.cell <- function(i) {
	  # each focal cell and its neighboring cells form a region; get a subset of community for the region
	  cols.beta <- which(nbsmat[rows.beta[i], ] == 1)
	  comdata.sub <- comdata[c(rows.beta[i], cols.beta), -c(1:2)]
	  id <- colSums(comdata.sub)>0
	  comdata.sub <- comdata.sub[,id]
	  
	  # get number of neighboring cells; richness in the focal cell, richness in the region
	  cell_info <- c(nrow(comdata.sub)-1, sum(comdata.sub[1, ]), ncol(comdata.sub))
	  names(cell_info) <- c("n.neigh", "alpha.focal", "gamma")
	  
	  # simplify the phylogeny to reduce calculation
	  if(method %in% c("phylosor","unifrac")){
	    if(!root) {
	      tree.sub <- drop.tip(tree, setdiff(colnames(comdata)[-c(1:2)], colnames(comdata.sub)), subtree=TRUE, trim.internal=TRUE)
	    }
	    if(root) {
	      tree.sub <- keep.tip(tree, c(colnames(comdata.sub), "root"))
	    }
	    tree.sub$root.edge <- 0
	  }
	  
	  if(method %in% c("Dpw","Dnn") & any(!is.na(tree))){
	    tree.sub <- keep.tip(tree, colnames(comdata.sub))
	  }
	  
	  if(any(!is.na(trait))){
	    trait.sub <- trait[rownames(trait) %in% colnames(comdata.sub), ]
	  }
	  
	  # local null model: randomize the species of the region  in the phylogeny and trait data
	  if(null=="local"){
	    if(any(!is.na(tree))) {
	      id <- tree.sub$tip.label %in% colnames(comdata.sub)
	      tree.sub$tip.label[id] <- sample(tree.sub$tip.label[id])
	    }
	    if(any(!is.na(trait))) rownames(trait.sub) <- sample(rownames(trait.sub))
	  }
	  
	  # get the phylogentic or traits distance for calculate Dpw or Dnn if applicable
	  if(method %in% c("Dpw","Dnn")){
	    if(any(!is.na(tree))) {trait.sub.dist <- cophenetic(tree.sub)}
	    if(any(!is.na(trait))){trait.sub.dist <- as.matrix(dist(trait.sub))}
	  }
	  
	  # if want to keep same number of neighboring cells for all focal cells, resample neighboring cells
	  if(!is.na(ncell_resample) & ncell_resample < nrow(comdata.sub)-1){
	    id <- sample(2:nrow(comdata.sub), ncell_resample)
	    comdata.sub <- comdata.sub[c(1, id), ]
	  }
	  
	  beta.cell <- vector()
	  
	  # calculate averaged pairwise beta and components between focal cell and its neighboring cells
	  if(!beta.multi){
	    if(!alpha){
	      for(j in 2:nrow(comdata.sub)){
	        if(method=="sorensen") beta.1pair <- unlist(beta.pair(comdata.sub[c(1,j), ], index.family="sorensen"))
	        if(method=="jaccard") beta.1pair <- unlist(beta.pair(comdata.sub[c(1,j),], index.family="jaccard"))
	        if(method=="phylosor") beta.1pair <- unlist(phylo.beta.pair(comdata.sub[c(1,j),], tree=tree.sub, index.family="sorensen"))
	        if(method=="unifrac") beta.1pair <- unlist(phylo.beta.pair(comdata.sub[c(1,j),], tree=tree.sub, index.family="jaccard"))
	        if(method=="funcsor") {
	          beta.1pair <- try(unlist(functional.beta.pair(comdata.sub[c(1,j), ],traits=trait.sub, index.family="sorensen")), silent=TRUE)
	          if(inherits(beta.1pair, "try-error")) beta.1pair <- rep(NA, 3)
	        }
	        if(method=="funcjac") {
	          beta.1pair <- try(unlist(functional.beta.pair(comdata.sub[c(1,j), ], traits=trait.sub, index.family="jaccard")), silent=TRUE)
	          if(inherits(beta.1pair, "try-error")) beta.1pair <- rep(NA,3)
	        }
	        if(method=="Dpw") {
	          beta.1pair <- as.vector(comdist(comdata.sub[c(1,j),], dis=trait.sub.dist))
	          names(beta.1pair) <- "Dpw"
	        }
	        if(method=="Dnn") {
	          beta.1pair <- as.vector(comdistnt(comdata.sub[c(1,j),], dis=trait.sub.dist))
	          names(beta.1pair) <- "Dnn"
	        }
	        beta.cell <- rbind(beta.cell, beta.1pair)
	      }
	    }
	    
	    # besides pairwise beta, also get shared richness, alpha1, alpha2, max.notshared and min.notshared
	    if(alpha){
	      for(j in 2:nrow(comdata.sub)){
	        if(method=="sorensen"){
	          comdata.sub.core <- betapart.core(comdata.sub[c(1,j), ])
	          beta.1pair <-  unlist(beta.pair(comdata.sub.core, index.family="sorensen"))
	          beta.1pair <- c(beta.1pair, c(comdata.sub.core[[5]][c(2,1,4)], comdata.sub.core[[8]][2], comdata.sub.core[[9]][2]))
	        }
	        if(method=="jaccard"){
	          comdata.sub.core <- betapart.core(comdata.sub[c(1,j), ])
	          beta.1pair <-  unlist(beta.pair(comdata.sub.core, index.family="jaccard"))
	          beta.1pair <- c(beta.1pair, c(comdata.sub.core[[5]][c(2,1,4)], comdata.sub.core[[8]][2], comdata.sub.core[[9]][2]))
	        }
	        if(method=="phylosor"){
	          comdata.sub.core <- phylo.betapart.core(comdata.sub[c(1,j),], tree=tree.sub)
	          comdata.sub.pd <- pd(comdata.sub[c(1,j+1),], tree=tree.sub)
	          beta.1pair <-  unlist(phylo.beta.pair(comdata.sub.core, index.family="sorensen"))
	          beta.1pair <- c(beta.1pair, c(comdata.sub.core[[3]][1], comdata.sub.pd[1:2,1], comdata.sub.core[[5]][1], comdata.sub.core[[6]][1]))
	        }
	        if(method=="unifrac"){
	          comdata.sub.core <- phylo.betapart.core(comdata.sub[c(1,j+1),], tree=tree.sub)
	          comdata.sub.pd <- pd(comdata.sub[c(1,j+1),], tree=tree.sub)
	          beta.1pair <-  unlist(phylo.beta.pair(comdata.sub.core, index.family="jaccard"))
	          beta.1pair <- c(beta.1pair, c(comdata.sub.core[[3]][1], comdata.sub.pd[1:2,1], comdata.sub.core[[5]][1], comdata.sub.core[[6]][1]))
	        }
	        beta.cell <- rbind(beta.cell, beta.1pair)
	      }
	      colnames(beta.cell)[4:8] <-  c("shared", "alpha1", "alpha2", "max.notshared", "min.notshared")
	    }
	    
	    beta.cell <- colMeans(beta.cell, na.rm=TRUE)
	  }
	  
	  # calculate multiple-site beta-diversity
	  if(beta.multi){
	    if(!alpha){
	      if(method=="sorensen") beta.cell <- unlist(beta.multi(comdata.sub, index.family="sorensen"))
	      if(method=="jaccard") beta.cell <- unlist(beta.multi(comdata.sub, index.family="jaccard"))
	      if(method=="phylosor") beta.cell <- unlist(phylo.beta.multi(comdata.sub, tree=tree.sub, index.family="sorensen"))
	      if(method=="unifrac") beta.cell <- unlist(phylo.beta.multi(comdata.sub, tree=tree.sub, index.family="jaccard"))		  
	    }
	    
	    # besides multiple-site beta, also get averaged shared richness, alpha, max.notshared and min.notshared
	    if(alpha){
	      if(method=="sorensen"){
	        comdata.sub.core <- betapart.core(comdata.sub)
	        beta.com <- c(mean(as.dist(comdata.sub.core[[5]])), mean(diag(comdata.sub.core[[5]])), 
	                      mean(as.dist(comdata.sub.core[[8]])), mean(as.dist(comdata.sub.core[[9]])))
	        beta.cell <- unlist(beta.multi(comdata.sub.core, index.family="sorensen"))
	        beta.cell <- c(beta.cell, beta.com)
	      }
	      if(method=="jaccard"){
	        comdata.sub.core <- betapart.core(comdata.sub)
	        beta.com <- c(mean(as.dist(comdata.sub.core[[5]])), mean(diag(comdata.sub.core[[5]])),
	                      mean(as.dist(comdata.sub.core[[8]])), mean(as.dist(comdata.sub.core[[9]])))
	        beta.cell <- unlist(beta.multi(comdata.sub.core, index.family="jaccard"))
	        beta.cell <- c(beta.cell, beta.com)
	      }
	      if(method=="phylosor"){
	        comdata.sub.core <- phylo.betapart.core(comdata.sub, tree=tree.sub)
	        beta.com <- c(mean(comdata.sub.core[[3]]), comdata.sub.core[[1]]/nrow(comdata.sub), 
	                      mean(comdata.sub.core[[5]]), mean(comdata.sub.core[[6]]))
	        beta.cell <- unlist(phylo.beta.multi(comdata.sub.core, index.family="sorensen"))
	        beta.cell <- c(beta.cell, beta.com)
	      }
	      if(method=="unifrac"){
	        comdata.sub.core <- phylo.betapart.core(comdata.sub, tree=tree.sub)
	        beta.com <- c(mean(comdata.sub.core[[3]]), comdata.sub.core[[1]]/nrow(comdata.sub),
	                      mean(comdata.sub.core[[5]]), mean(comdata.sub.core[[6]]))
	        beta.cell <- unlist(phylo.beta.multi(comdata.sub.core, index.family="jaccard"))
	        beta.cell <- c(beta.cell, beta.com)
	      }
	      names(beta.cell)[4:7] <-  c("shared", "alpha", "max.notshared", "min.notshared")
	    }
	  }
	  
	  # combine cell information with beta values
	  beta.cell <-  c(cell_info, round(beta.cell, 4))
	  
	  return(beta.cell)
	}
	
	# calculate beta for all cells
	beta.res <- vector()
	for(i in 1:length(rows.beta)) {
	  beta.cell <- try(get_beta.cell(i), silent = TRUE)
	  if(inherits(beta.cell, "try-error")) beta.cell <- NA
	  beta.res <- rbind(beta.res, beta.cell)
	}

	rownames(beta.res) <- rows.beta
	beta.res <- data.frame("ID"= as.numeric(rownames(comdata)[rows.beta]), beta.res)
	return(beta.res)
}
