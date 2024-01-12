library (argparser)

args = arg_parser("Input Arguments")
args = add_argument( args, "--input" , help = "Path to RDS input containing the multiOmic Data")
args = add_argument( args, "--nfact", help = "Number of Factors" , type = "integer")
args = add_argument( args, "--outputDir", help = "Path to output result")
args = add_argument( args, "--pathRgcca", help = "Path to RGCCA factory R files") 

argv = parse_args(args)

omics = readRDS (argv$input )
num_factors = argv$nfact



####################
## Load functions ##
####################
fun_dir = dirname(getwd())
# RGCCA
devtools::load_all ( file.path( argv$pathRgcca , "/RGCCA_factory/RGCCA/."))
source ( file.path ( argv$pathRgcca , "/additional_functions/compute_SGCCA_multi-Copy1.R"))
source ( file.path ( argv$pathRgcca , "/additional_functions/scale_test.R"))



blocks = omics

result_sgcca = rgcca(blocks  = blocks  , 
                     ncomp = rep(num_factors, length(blocks)), 
                     scheme = "centroid", 
                     scale = TRUE, 
                     init = "svd",
                     bias = TRUE, 
                     tol = 1e-08, 
                     verbose = T,
                     superblock = T) 


saveRDS ( result_sgcca , file = paste0 (argv$outputDir , "/06_factorizations_SGCCA_dry.RDS" ))
