#
# Use shrunken CCA, perhaps supervised, to pick loci for exploration via heatmap
#
# Arguments: 
#
# object: an eSet of some sort (this function was constructed for MethyLumiSets)
# design: a design matrix (optional -- can be specified instead of picking covs)
# covs: alternatively, specify the pData column names to include in the model
# loci: restrict to certain loci for expression/methylation (default: use all)
# samples: which samples should be used for the analysis? (default: all of them)
# y: an optional outcome with which we would like to associate canonical vectors
# grp.x: optional smoothing factor for ordering/partitioning the covariates in x
# grp.z: optional smoothing factor for grouping the specified loci (or all loci)
# pen.x: the lasso or fused lasso penalty to apply against x while computing CCA
# pen.z: the lasso or fused lasso penalty to apply against z while computing CCA
# niter: maximum iterations to be performed while computing CCA (default is 100)
# method: what clustering method to use? (options are modal, ward, and rpmm)
# standardize: should the columns of X and Z be centered and scaled into Z-vals?
# title: specify a title for the resulting heatmap here, if that's what you want
# colorbar: optional colorcoded matrix with colnames to display covariate values
# locbar: optional color-coded matrix with colnames to display locus codings
# parallel: should we attempt to run in parallel if possible? (needs multicore)
# plot: should a heatmap be generated? (as might be expected, default is TRUE)
# forced: any columns in the colorbar to force into the plot regardless of CCA
# annot.CpG: try to automatically annotate CpG islands? (only for methylation)
# annot.promoter: try to automatically annotate promoters? (for methylation)
# x.dendro: plot a dendrogram for the covariates in X? (default: true)
# z.dendro: plot a dendrogram for the covariates in Z? (default: false)
# colorcode.x: should the program attempt to assign X colors by correlation?(T)
# rotate: transpose the heatmap? (default: samples=rows, loci=columns)
#
# Value:
# 
# The results of the CCA run and (passively, only if default plot=T) a heatmap.
#
# Warning:
#
# You don't actually want the fused lasso, you want the grouped lasso here...
# by and large, attempts to use grp.x and grp.z will not make you happy.
#
autoheatmap <- function(object,design=NULL,covs=NULL,loci=NULL,samples=NULL,y=NULL,grp.x=NULL,grp.z=NULL,pen.x=0.3,pen.z=0.2,niter=100,method="ward",standardize=T,title=NULL,colorbar=NULL,locbar=NULL,parallel=F,plot=T,forced=NULL,annot.CpG=F,annot.promoter=F,x.dendro=T,z.dendro=F,rotate=F,colorcode.x=T,how='betas'){ 

  require(impute) # eh.
  require(grplasso) 
  require(flashClust)
  require(heatmap.plus)
  
  # must have design matrix if not an eSet
  if(!is(object, 'eSet')) {
    stopifnot(class(object) == 'matrix')
    stopifnot(class(design) == 'matrix')
    IDs <- colnames(object)
  } else {
    if(is.element('ID', varLabels(object))) {
      IDs <- pData(object)[['ID']]
    } else {
      IDs <- sampleNames(object)
    }
    # fix this to use an arbitrary coordinate-based annotation instead...
    if(annot.CpG||annot.promoter) {
      annot <- annotation(object)
      require(paste(annot,'db',sep='.'), character.only=T)
    }
  }

  # just in case we have extras... (!)
  if(is.null(covs) && is(object, 'eSet')) {
    covs <- intersect(covs, varLabels(object)) 
  }

  # any subsetting to be done on the columns?
  if(is.character(samples)) samples <- match(samples, sampleNames(object))
  if(is.null(samples)) {
    if(is(object, 'eSet')) s <- 1:dim(object)[2]
    else s <- ncol(object)
  } else {
    s <- samples
  }

  # generate an appropriate matrix for X
  if(!is.null(design)) {
    X <- design
    standardize <- F
  } else {
    X <- data.matrix(pData(object)[s,covs])
  }
  X <- impute.knn(X)$data

  # generate an appropriate matrix for Z
  if(is.null(loci) && is(object, 'eSet')) {
    loci <- featureNames(object)
  } else if(is.null(loci) && is(object, 'matrix')) {
    loci <- rownames(object)
  } else if(!is.null(loci) && is(object, 'eSet')) {
    loci <- intersect(loci, featureNames(object))
  } else {
    loci <- intersect(loci, rownames(object))
  }
  if(is(object, 'eSet')) {
    if(how=='betas') z <- betas(object)[loci,s]
    if(how=='exprs') z <- exprs(object)[loci,s]
    if(how=='logit') z <- beta.logit(betas(object)[loci,s], n=0.9999999)
    if(how=='probit') z <- beta.probit(betas(object)[loci,s], n=0.9999999)
    if(how=='arcsin') z <- arcsinize(betas(object)[loci,s])
    else z <- exprs(object)[loci,s]
  } else {
    z <- object[loci,s]
  }
  colnames(z) <- IDs
  Z <- impute.knn(t(z))$data

  # summarize the results as appropriate
  CCA.result <- CCA(x=X, z=Z, y=y, chromz=grp.z, chromx=grp.x, niter=niter,
                    typex=ifelse(is.null(grp.x), 'standard', 'ordered'),
                    typez=ifelse(is.null(grp.z), 'standard', 'ordered'),
                    penaltyx=pen.x, penaltyz=pen.z, standardize=standardize)

  # positive/negative correlated loadings
  pos.loci <- which(CCA.result$v > 0)
  neg.loci <- which(CCA.result$v < 0)
  if(!is.null(locbar)) locbar <- locbar[ c(pos.loci, neg.loci), ] # original
  correlations <- round( 100 * (CCA.result$v[ c(pos.loci, neg.loci) ] + 1)/2 )
  lcor.pal <- colorRampPalette(c("darkred","red","white","green","darkgreen"), 
                               space="Lab")(100)
  lcor.bar <- lcor.pal[ correlations ]

  # annotation
  if(annot.CpG){
    orig.cols <- colnames(locbar)
    # is.cpg <- unlist(mget(loci[c(pos.loci,neg.loci)],
    #                  get(paste(annot,'ISCPGISLAND',sep=''))))
    is.cpg <- fData(object)$CGIS_IRIZARRY_2009[c(pos.loci,neg.loci)]
    locbar <- cbind(locbar, matrix(c('#FFFFFF','#0000FF')[is.cpg+1], ncol=1))
    colnames(locbar) <- c(orig.cols, 'CpGI?')
  }
  if(annot.promoter){
    warning("Promoter/exon annotation is not yet implemented (feel free...)")
  }
  if(!is.null(locbar)) {
    orig.cols <- colnames(locbar)
    locbar <- cbind( lcor.bar, locbar )  
    colnames(locbar) <- c('Strength', orig.cols)
  } else {
    locbar <- cbind(lcor.bar, lcor.bar)
    if(rotate) 
      colnames(locbar) <- c('','')
    else 
      colnames(locbar) <- c('(canonical)','Association')
  }

  # try & dig up the symbols
  if( is(object, 'eSet') ) { 
    teh.loci <- intersect(loci,featureNames(object))[ c(pos.loci, neg.loci) ]
  } else { 
    teh.loci <- c(pos.loci, neg.loci) 
  }  
  Z.subset <- Z[,teh.loci]
  if( is(object, 'eSet') ) {
    symbol.env <- paste(annotation(object), 'SYMBOL', sep='')
    colnames(Z.subset) <- unlist(mget(teh.loci, get(symbol.env)))
  }

  # positive/negative correlated loadings
  pos.covs <- which(CCA.result$u > 0)
  neg.covs <- which(CCA.result$u < 0)
  if(is.character(forced)) forced <- grep(forced, colnames(colorbar))
  if(is.null(design) && !is.null(colorbar)) {
    keep <- unique(c(pos.covs, neg.covs, forced))
    if(length(keep)<2) {
      colorbar <- cbind(colorbar, rep('#FFFFFF', nrow(colorbar)))
      colnames(colorbar)[ncol(colorbar)] <- ' '
      keep <- c(keep, ncol(colorbar))
    }
  } else if(!is.null(colorbar)) {
    keep <- 1:ncol(colorbar)
  }

  # colorcode canonical associations?
  if(colorcode.x) {
    for(i in pos.covs) colorbar[X[,i]!=0, i] <- '#AA0000'
    for(i in neg.covs) colorbar[X[,i]!=0, i] <- '#00AA00'
  }

  # prepare the heatmap (default is 'methylation colors' blue->black->yellow)
  pal <- colorRampPalette(c("blue","black","yellow"), space="Lab")(255)
  hward <- function(d) hclust(d, method='ward')
  if(!z.dendro) {
    Z.subset <- Z.subset[, hward(dist(t(Z.subset), 'euclidean'))$order]
    zdend <- NA
  } else {
    zdend <- NULL
  }
  if(!x.dendro) {
    Z.subset <- Z.subset[hward(dist(Z.subset, 'euclidean'))$order, ]
    xdend <- NA
  } else {
    xdend <- NULL
  }
  if( plot ) {
    if(!is.null(colorbar)) {
      if( rotate ) {
        heatmap.plus(t(Z.subset), RowSideColors=locbar, Rowv=zdend,
                     ColSideColors=colorbar[,keep], Colv=xdend, hclustfun=hward,
                     scale="none", col=pal, main=paste("\n\n",title))
      } else {
        heatmap.plus(Z.subset, ColSideColors=locbar, Colv=zdend,
                     RowSideColors=colorbar[,keep], Rowv=xdend, hclustfun=hward,
                     scale="none", col=pal, main=paste("\n\n",title))
      }
    } else {
      if( rotate ) {
        heatmap.plus(t(Z.subset),RowSideColors=locbar,Rowv=zdend,Colv=xdend,
                     hclustfun=hward, scale="none", col=pal, 
                     main=paste("\n\n",title))
      } else {
        heatmap.plus(Z.subset, ColSideColors=locbar,Colv=zdend,Rowv=xdend,
                     hclustfun=hward, scale="none", col=pal, 
                     main=paste("\n\n",title))
      }
    }
  }

  return(CCA.result)
}
