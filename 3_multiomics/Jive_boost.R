#!/usr/local/bin/Rscript
#.libPaths( c( .libPaths(), "/home/amazigh.mokhtari/NeuroDev_ADD/R/r_3.6.0/") )

library (argparser)
library (r.jive)
library (dplyr)
library(corrplot)


args = arg_parser("Input Arguments")
args = add_argument( args, "--input", , help = "Path to RDS input containing the multiOmic Data")
args = add_argument( args, "--nfact", help = "Number of Factors" , type = "integer")
args = add_argument( args, "--outputDir", help = "Path to output result")


argv = parse_args(args)


#########################
#########################
#########################

#### Function to extract Jive Results
Jive_extract = function (JIVE_dry) {
    ### Extract Jive Results 
    rankJV   = JIVE_dry$rankJ
    rankIV.v = JIVE_dry$rankA   

    J  = numeric(0)
    ng = 0

    factors_jive = features_w_jive = list()

    for(j in 1:length(JIVE_dry$data)){

        J  = rbind(J,JIVE_dry$joint[[j]]);
        ng = c(ng,dim(JIVE_dry$joint[[j]])[1])
      }

    svd.o = svd(J);
    jV    = svd.o$v %*% diag(svd.o$d);

    for(j in 1:length(JIVE_dry$data)){
        features_w_jive[[j]]  =  svd.o$u[(1+sum(ng[1:j])):sum(ng[1:j+1]),1:rankJV]; ###error in dimension
        rownames(features_w_jive[[j]]) = colnames(omics [[names(JIVE_dry$data[j] )]])   #rownames( t(JIVE_dry$data[[j]]))
        colnames(features_w_jive[[j]]) = paste0("JIVE_" , 1:num_factors )
      }

    factors_jive = jV[,1:rankJV]
    colnames(factors_jive) = paste0("JIVE_" , 1:10 )


    return (list( factors_jive = factors_jive, features_w_jive = features_w_jive )) }

########################

cor.mtest = function(mat, ...) {
    mat = as.matrix(mat)
    n = ncol(mat)
    p.mat = matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp = cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] = p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) = rownames(p.mat) <- colnames(mat)
  p.mat
}



cor_mat = function (x){
    M = cor(x)
p.mat = cor.mtest(x)
col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
plot = corrplot(M, method="color", col=col(200),  
         type="upper", 
         addCoef.col = "black", # Ajout du coefficient de corrélation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # Cacher les coefficients de corrélation sur la diagonale
         diag=FALSE 
         )
# get top fact cor to Group 
    pval_group = p.mat [ 1:nrow(p.mat)-1, ncol(p.mat)] 
    
    name_fact_ = names(which.min(pval_group [ pval_group < 0.05]))
    if (length(name_fact_) != 0 ) {
         top.fact =  name_fact_   
    }
    
    return( list (cor.mat = M , p.mat = p.mat , plot = plot , top.fact = top.fact))
}


#### Get top features 


get_top_feature = function (features_w , top_cor_ ) {
    ### remove the Cov table if exist    
    lapply (1:length(features_w) , function (x) {   X =  features_w [[x]] %>% 
                         .[order ( - abs(.[,  top_cor_  ])  ) ,] ; # order the features according to the top factor that corr the most with phenotype
                                if ( grepl("miR" , rownames(X) [1] ) ){ ### if omic == miRNA %>% get the top 2% (~ 14 miRNA)
                                    X %>% head (n = round (nrow(.) *0.02 )  ) %>% 
                                    rownames
                                }else{
                                    X %>% head (n = round (nrow(.) *0.1 )  ) %>%  ### if omic is else  %>% get the top 10% 
                                    rownames
                                                                                                              }
                                                                                                              } )  

}
            

#########################
###### Variables ########
#########################




omics = readRDS (argv$input )
num_factors = argv$nfact


                                       
factorizations_jive = jive( lapply (omics,t) , 
                            rankJ   = num_factors, 
                            rankA   = rep ( num_factors, length( omics ) ), 
                            method  = "given", 
                            conv    = "default", 
                            maxiter = 1000, 
                            showProgress = FALSE)


JIVE_Fact   =  Jive_extract(factorizations_jive) [["factors_jive"]] 



features_w_jive =   Jive_extract(factorizations_jive) [["features_w_jive"]]   
features_w_jive = features_w_jive[1:(length(features_w_jive)-1) ]

cor_meth_jive = lapply ( JIVE_Fact , function (y) cor.mtest (  cbind(y , Group )) [1:10 , 11] )
top_cor_jive =  lapply ( 1:length(cor_meth_jive) ,  function (y) names(cor_meth_jive [[y]]) [which.min(cor_meth_jive [[y]])]   )
                                             
        
                        
                        
Top_feature = get_top_feature (features_w_jive , top_cor_jive )