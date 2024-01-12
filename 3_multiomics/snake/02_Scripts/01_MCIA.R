library (argparser)
library (omicade4)


args = arg_parser("Input Arguments")
args = add_argument( args, "--input", help = "Path to RDS input containing the multiOmic Data")
args = add_argument( args, "--nfact", help = "Number of Factors" , type = "integer")
args = add_argument( args, "--outputDir", help = "Path to output result")


argv = parse_args(args)


omics = readRDS (argv$input )
num_factors = argv$nfact



factorizations_mcia = mcia( omics  , cia.nf = num_factors)

saveRDS ( factorizations_mcia , file = file.path (argv$outputDir, "01_factorizations_mcia_dry.RDS" ))
### Get the factors Matrix
factors_mcia = as.matrix(factorizations_mcia$mcoa$SynVar)
colnames (factors_mcia) = paste0 ( "MCIA_", 1:10)
saveRDS(factors_mcia , file = file.path  (argv$outputDir, "factors_mcia.RDS") )
