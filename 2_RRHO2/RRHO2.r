library(dplyr)
require(parallel)
require(VennDiagram)
require(RRHO)

## Compute the overlaps between two *numeric* lists:
numericListOverlap<- function(sample1, sample2, stepsize, method="hyper", tol = 0.5, offset = 1 , mcores = ncores){
  n<- length(sample1)
  
  overlap_hyper <- function(a,b) {
    count<-as.integer(sum(as.numeric(sample1[1:a] %in% sample2[1:b])))    
    signs<- 1L
    log.pval<- -phyper(q=count-1, m=a, n=n-a+1, k=b, lower.tail=FALSE, log.p=TRUE)    
    #log.pval[is.na(log.pval)]<-0
    
    return(c(counts=count, 
             log.pval=as.numeric(log.pval),
             signs=as.integer(signs)
    ))    
  }
  
  overlap_fisher <- function(a,b) {
    s1 <- sample1[1:a]
    s2 <- sample2[1:b]
    lenA <- as.integer(sum(as.numeric(s1 %in% s2))) 
    lenB <- length(s1)
    lenC <- length(s2)
    Odds<-((lenA+offset)*(n-lenB-lenC+lenA+offset))/((lenC-lenA+offset)*(lenB-lenA+offset))
    logOdds <- log(abs(Odds))*sign(Odds)
    signs<- 1L
    
    return(c(counts=lenA, 
             log.pval=as.numeric(logOdds),
             signs=as.integer(signs)
    ))    
  }
  
  indexes<- expand.grid(i=seq(1,n,by=stepsize), j=seq(1,n,by=stepsize))
  if(method=="hyper"){
    #overlaps<- apply(indexes, 1, function(x) overlap_hyper(x['i'], x['j']))
    overlaps<- parallel::mcmapply(overlap_hyper , a = indexes[ ,1] , b = indexes[ ,2] , mc.cores = mcores)
  } else if(method=="fisher"){
    #overlaps<- apply(indexes, 1, function(x) overlap_fisher(x['i'], x['j']))
    overlaps<- parallel::mcmapply(overlap_fisher , a = indexes[ ,1] , b = indexes[ ,2] , mc.cores = mcores)
  }
  
  nrows<- sqrt(ncol(overlaps))
  matrix.counts<- matrix(overlaps['counts',], ncol=nrows)  
  matrix.log.pvals<- matrix(overlaps['log.pval',], ncol=nrows)  
  matrix.signs<- matrix(overlaps['signs',], ncol=nrows)  
  
  return(list(counts=matrix.counts, log.pval=matrix.log.pvals, signs = matrix.signs))
}


RRHO2_initialize <- function (list1, list2, labels = NULL, log10.ind = FALSE, multipleTesting="none", boundary = 0.1)
{
  stepsize = RRHO:::defaultStepSize(list1, list2)
  if (any(duplicated(list1[,1])))
    stop("Non-unique gene identifier found in list1")
  if (any(duplicated(list2[,1])))
    stop("Non-unique gene identifier found in list2")
  if(!is.null(labels)){
    stopifnot(length(labels) == 2)
  }
  if(any(is.na(list1))){
    stop("NA value exists in list1, please remove them.")
  }
  if(any(is.na(list2))){
    stop("NA value exists in list2, please remove them.")
  }
  if(!all(list1[, 1] %in% list2[, 1]) | !all(list2[, 1] %in% list1[, 1])){
    stop("The gene names of the two lists must be identical.")
  }
  
  list1 <- list1[order(list1[, 2], decreasing = TRUE), ]
  list2 <- list2[order(list2[, 2], decreasing = TRUE), ]
  nlist1 <- length(list1[, 1])
  nlist2 <- length(list2[, 1])
  N <- max(nlist1, nlist2)
  
  .hypermat_normal<- numericListOverlap(list1[, 1], list2[, 1], stepsize, method='hyper' , mcores = 1)
  hypermat_normal<- .hypermat_normal$log.pval
  
  .hypermat_flipX <- numericListOverlap(rev(list1[, 1]), list2[, 1], stepsize, method='hyper', mcores = 1)
  hypermat_flipX <- .hypermat_flipX$log.pval

  if(multipleTesting=="none"){
    ;
  } else if(multipleTesting=="BH"){
    hypermatvec  <- matrix(hypermat_normal,nrow=nrow(hypermat_normal)*ncol(hypermat_normal),ncol=1)
    hypermat.bhvec  <- p.adjust(exp(-hypermatvec),method="BH")
    hypermat_normal <- matrix(-log(hypermat.bhvec),nrow=nrow(hypermat_normal),ncol=ncol(hypermat_normal))       	  	
    
    hypermatvec  <- matrix(hypermat_flipX,nrow=nrow(hypermat_flipX)*ncol(hypermat_flipX),ncol=1)
    hypermat.bhvec  <- p.adjust(exp(-hypermatvec),method="BH")
    hypermat_flipX <- matrix(-log(hypermat.bhvec),nrow=nrow(hypermat_normal),ncol=ncol(hypermat_normal))       	  	
  } else if(multipleTesting=="BY"){
    hypermatvec  <- matrix(hypermat_normal,nrow=nrow(hypermat_normal)*ncol(hypermat_normal),ncol=1)
    hypermat.byvec  <- p.adjust(exp(-hypermatvec),method="BY")
    hypermat_normal <- matrix(-log(hypermat.byvec),nrow=nrow(hypermat_normal),ncol=ncol(hypermat_normal))       	  	
    
    hypermatvec  <- matrix(hypermat_flipX,nrow=nrow(hypermat_flipX)*ncol(hypermat_flipX),ncol=1)
    hypermat.byvec  <- p.adjust(exp(-hypermatvec),method="BY")
    hypermat_flipX <- matrix(-log(hypermat.byvec),nrow=nrow(hypermat_normal),ncol=ncol(hypermat_normal))       	  	
  } else {
    stop("no such multiple testing procedure, please check the multipleTesting argument in the help file for RRHO2_initialize")
  }
  
  
  stepList1 <- seq(1, nlist1, stepsize)
  stepList2 <- seq(1, nlist2, stepsize)
  
  len1 <- length(stepList1)
  len2 <- length(stepList2)
  
  lenStrip1 <- round(len1*boundary)
  lenStrip2 <- round(len2*boundary)
  
  boundary1 <- sum(list1[stepList1,2] > 0)
  boundary2 <- sum(list2[stepList2,2] > 0)
  
  hypermat <- matrix(NA, nrow = nrow(hypermat_normal) + lenStrip1,
                     ncol = ncol(hypermat_normal) + lenStrip2)
  
  ## d1d2, quadrant I
  hypermat[lenStrip1 + (boundary1+1):len1, lenStrip2 + (boundary2+1):len2] <- hypermat_normal[(boundary1+1):len1, (boundary2+1):len2]

  ## u1d2, quadrant II
  hypermat[1:boundary1, lenStrip2 + (boundary2+1):len2] <- hypermat_flipX[len1:(len1 - boundary1 + 1),(boundary2+1):len2]
  
  ## u1u2, quadrant III
  hypermat[1:boundary1, 1:boundary2] <- hypermat_normal[1:boundary1,1:boundary2]
  
  ## u1d2, quadrant IV
  hypermat[lenStrip1 + (boundary1+1):len1, 1:boundary2] <- hypermat_flipX[(len1 - boundary1):1,1:boundary2]
  
  if(any(is.infinite(hypermat[!is.na(hypermat)]))	){
    warning("Inf was generated because of the multiple testing procedure. I.e., the multiple testing procedure cannot handle extreme small p-values. Suggest to use raw p-value (multipleTesting='none')")
  }
  
  if (log10.ind){
    hypermat <- hypermat * log10(exp(1))
  }
  
  #### dd: down in 1 and down in 2
  maxind.dd <- which(max(hypermat[lenStrip1 + (boundary1+1):len1, lenStrip2 + (boundary2+1):len2],
                         na.rm = TRUE) == hypermat, arr.ind = TRUE)
  maxind.dd <- maxind.dd[maxind.dd[,1]>=lenStrip1 + (boundary1+1) & maxind.dd[,1]<=lenStrip1 +len1 & 
                           maxind.dd[,2]>=lenStrip2 + (boundary2+1) & maxind.dd[,2]<=lenStrip2 + len2,]
  #
  if(!is.null(dim(maxind.dd))){
    maxind.dd <- maxind.dd[1, ]
  }
  
  indlist1.dd <- seq(1, nlist1, stepsize)[maxind.dd[1] - lenStrip1]
  indlist2.dd <- seq(1, nlist2, stepsize)[maxind.dd[2] - lenStrip2]
  gene_list1_dd <- list1[indlist1.dd:nlist1, 1]
  gene_list2_dd <- list2[indlist2.dd:nlist2, 1]
  gene_list_overlap_dd <- intersect(gene_list1_dd,
                                    gene_list2_dd)
  genelist_dd <- list(gene_list1_dd=gene_list1_dd, 
                      gene_list2_dd=gene_list2_dd,
                      gene_list_overlap_dd=gene_list_overlap_dd
  )
  #### uu: up in 1 and up in 2
  maxind.uu <- which(max(hypermat[1:boundary1, 1:boundary2],
                         na.rm = TRUE) == hypermat, arr.ind = TRUE)
  maxind.uu <- maxind.uu[maxind.uu[,1]>=1 & maxind.uu[,1]<=boundary1 & maxind.uu[,2]>=1 & maxind.uu[,2]<=boundary2,]
  if(!is.null(dim(maxind.uu))){
    maxind.uu <- maxind.uu[1, ]
  }
  
  indlist1.uu <- seq(1, nlist1, stepsize)[maxind.uu[1]]
  indlist2.uu <- seq(1, nlist2, stepsize)[maxind.uu[2]]
  gene_list1_uu <- list1[1:indlist1.uu, 1]
  gene_list2_uu <- list2[1:indlist2.uu, 1]
  gene_list_overlap_uu <- intersect(gene_list1_uu,
                                    gene_list2_uu)
  genelist_uu <- list(gene_list1_uu=gene_list1_uu, 
                      gene_list2_uu=gene_list2_uu,
                      gene_list_overlap_uu=gene_list_overlap_uu
  )
  #### ud: up in 1 and down in 2
  maxind.ud <- which(max(hypermat[1:boundary1, lenStrip2 + (boundary2+1):len2],
                         na.rm = TRUE) == hypermat, arr.ind = TRUE)
  #
  maxind.ud <- maxind.ud[maxind.ud[,1]>=1 & maxind.ud[,1]<=boundary1 & maxind.ud[,2]>= lenStrip2 + (boundary2+1) & maxind.ud[,2]<=lenStrip2 + len2,]
  if(!is.null(dim(maxind.ud))){
    maxind.ud <- maxind.ud[1, ]
  }
  
  indlist1.ud <- seq(1, nlist1, stepsize)[maxind.ud[1]]
  indlist2.ud <- seq(1, nlist2, stepsize)[maxind.ud[2] - lenStrip2]
  gene_list1_ud <- list1[1:indlist1.ud, 1]
  gene_list2_ud <- list2[indlist2.ud:nlist2, 1]
  gene_list_overlap_ud <- intersect(gene_list1_ud,
                                    gene_list2_ud)
  genelist_ud <- list(gene_list1_ud=gene_list1_ud, 
                      gene_list2_ud=gene_list2_ud,
                      gene_list_overlap_ud=gene_list_overlap_ud
  )
  
  #### du: down in 1 and up in 2
  maxind.du <- which(max(hypermat[lenStrip1 + (boundary1+1):len1, 1:boundary2],
                         na.rm = TRUE) == hypermat, arr.ind = TRUE)
  #
  maxind.du <- maxind.du[maxind.du[,1]>=lenStrip1 + (boundary1+1) & maxind.du[,1]<=lenStrip1 + len1 & maxind.du[,2]>=1 & maxind.du[,2]<=boundary2,]
  if(!is.null(dim(maxind.du))){
    maxind.du <- maxind.du[1, ]
  }
  
  indlist1.du <- seq(1, nlist1, stepsize)[maxind.du[1] - lenStrip1]
  indlist2.du <- seq(1, nlist2, stepsize)[maxind.du[2]]
  gene_list1_du <- list1[indlist1.du:nlist1, 1]
  gene_list2_du <- list2[1:indlist2.du, 1]
  gene_list_overlap_du <- intersect(gene_list1_du,
                                    gene_list2_du)
  genelist_du <- list(gene_list1_du=gene_list1_du, 
                      gene_list2_du=gene_list2_du,
                      gene_list_overlap_du=gene_list_overlap_du
  )
  
  
  result <- list(hypermat=hypermat, 
                 labels=labels,
                 log10.ind=log10.ind,
                 method='hyper',
                 genelist_uu=genelist_uu,
                 genelist_dd=genelist_dd,
                 genelist_ud=genelist_ud, 
                 genelist_du=genelist_du
  )
  return(result)
}

### get best log10(p-val) of each quadrent 
best_rrho2_values = function (rrho2Obj) {
  
  max_na_row = which (is.na(rrho2Obj$hypermat [1,])) %>% max 
  max_na_col = which (is.na(rrho2Obj$hypermat [,1])) %>% max 
  min_na_row = which (is.na(rrho2Obj$hypermat [1,])) %>% min 
  min_na_col = which (is.na(rrho2Obj$hypermat [,1])) %>% min 
  
  return (
    list (
      best_UU = (as.numeric(unlist(rrho2Obj$hypermat [1:min_na_row, 1:min_na_col] )) %>% na.omit %>% max ),
      best_DD = (as.numeric(unlist(rrho2Obj$hypermat [max_na_row:nrow(rrho2Obj$hypermat) ,max_na_col:ncol(rrho2Obj$hypermat)] )) %>% na.omit %>% max ),
      best_UD = (as.numeric(unlist(rrho2Obj$hypermat [max_na_row:nrow(rrho2Obj$hypermat), 1:min_na_col] )) %>% na.omit %>% max ) ,
      best_DU = (as.numeric(unlist(rrho2Obj$hypermat [1:min_na_row , max_na_col:ncol(rrho2Obj$hypermat)] )) %>% na.omit %>% max ) )
  )
}

RRHO2_heatmap <- function(RRHO_obj, maximum=NULL, minimum=NULL, colorGradient=NULL, labels=NULL, ...)
{
  
  hypermat <- RRHO_obj$hypermat
  method <- RRHO_obj$method
  
  if(is.null(labels)){
    labels <- RRHO_obj$labels	
  }
  
  if(!is.null(maximum)){
    hypermat[hypermat>maximum] <- maximum
  } else {
    maximum <- max(hypermat,na.rm=TRUE)
  }
  
  if(!is.null(minimum)){
    hypermat[hypermat<minimum] <- minimum
  } else {
    minimum <- min(hypermat,na.rm=TRUE)
  }
  
  if(minimum > maximum){
    stop("minimum > maximum, please check these function arguments!")
  }
  
  color.bar <- function(lut, min, max=-min, 
                        nticks=11, 
                        ticks=seq(min, max, len=nticks), 
                        title='') {
    scale  <- (length(lut)-1)/(max-min)
    plot(c(0,10), c(min,max), type='n', bty='n', 
         xaxt='n', xlab='', yaxt='n', ylab='')
    mtext(title,2,2.3, cex=0.8)
    axis(2, round(ticks,0), las=1,cex.lab=0.8)
    for (i in 1:(length(lut)-1)) {
      y  <- (i-1)/scale + min
      rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
  }
  
  if(is.null(colorGradient)){
    jet.colors  <- colorRampPalette(
      c("#00007F", "blue", "#007FFF", "cyan", 
        "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    colorGradient <- jet.colors(101)
  }
  layout(matrix(c(rep(1, 6), 2), 1, 7, byrow = TRUE))
  
  breaks <- seq(minimum,maximum,length.out = length(colorGradient) + 1)
  image(hypermat, col = colorGradient,breaks=breaks,
        axes = FALSE, ...)
  
  if(!is.null(labels)){
    mtext(labels[2],2,0.5)
    mtext(labels[1],1,0.5)
  }
  
  if(method == "hyper"){
    atitle <- ifelse(RRHO_obj$log10.ind, "-log10(P-value)", "-log(P-value)")
    color.bar(colorGradient, min = minimum, max = maximum, nticks = 6, title = atitle)
  } else if (method == "fisher"){
    atitle <- "log Odds"
    color.bar(colorGradient, min = minimum, max = maximum, nticks = 6, title = atitle)
  } else {
    stop("internal error (1), please report this error to https://github.com/RRHO2/RRHO2/issues")
  }
  invisible(hypermat)
}


#################
###   RRHO2   ###
#################

#Data
data = readRDS('res_single_mRNA.rds')

female = data$diff_female
female$DDE = -log10(female$pvalue) * ifelse(female$log2FoldChange > 0, 1, -1)
female$Genes = rownames(female)
female = female[, c('Genes', 'DDE')]

male = data$diff_male
male$DDE = -log10(male$pvalue) * ifelse(male$log2FoldChange > 0, 1, -1)
male$Genes = rownames(male)
male = male[, c('Genes', 'DDE')]

male = male[male$Genes %in% female$Genes,]
female = female[female$Genes %in% male$Genes,]

rm(data)

#Computing
RRHO_obj <-  RRHO2_initialize(male, female, labels = c("male", "female"), log10.ind = T)

best_rrho2_values(RRHO_obj)

RRHO2_heatmap(RRHO_obj)
