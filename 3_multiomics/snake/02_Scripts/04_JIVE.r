library ( r.jive )
library (argparser)

args = arg_parser("Input Arguments")
args = add_argument( args, "--input", help = "Path to RDS input containing the multiOmic Data")
args = add_argument( args, "--nfact", help = "Number of Factors" , type = "integer")
args = add_argument( args, "--outputDir", help = "Path to output result")

argv = parse_args(args)

omics = readRDS (argv$input )
num_factors = argv$nfact

### JIVE



                                       
factorizations_jive = jive( lapply (omics,t) , 
                            rankJ   = num_factors, 
                            rankA   = rep ( num_factors, length( omics ) ), 
                            method  = "given", 
                            conv    = "default", 
                            maxiter = 1000, 
                            showProgress = FALSE)

                                            
saveRDS ( factorizations_jive , file = file.path (argv$outputDir, "04_factorizations_jive_dry.RDS" ))
