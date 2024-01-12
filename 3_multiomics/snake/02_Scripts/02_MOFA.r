library (argparser)
library(MOFA2)
library("rhdf5")
library(reticulate)

args = arg_parser("Input Arguments")
args = add_argument( args, "--input" , help = "Path to RDS input containing the multiOmic Data")
args = add_argument( args, "--nfact", help = "Number of Factors" , type = "integer")
args = add_argument( args, "--outputDir", help = "Path to output result")


argv = parse_args(args)

omics = readRDS (argv$input )

names (omics) = c("miRNA","mRNA","DNAm","Hamilton")
num_factors = argv$nfact


# Creating MOFA object from a list of matrice
MOFAobject             = create_mofa( lapply (omics , t ) )
DataOptions            = get_default_data_options(MOFAobject)
model_opts             = get_default_model_options(MOFAobject)
model_opts$num_factors = num_factors 
TrainOptions           = get_default_training_options(MOFAobject)


# Prepare the MOFA Object 
MOFAobject =  prepare_mofa( object = MOFAobject,
                      data_options = DataOptions,
                     model_options = model_opts,
                  training_options = TrainOptions )

# Configure python for MOFA2
use_basilisk = FALSE
reticulate::use_python("/opt/miniconda3/bin/python3.7", required=TRUE)
outfile = file.path(getwd(),"model2.hdf5")


MOFAobject_trained = run_mofa(MOFAobject, outfile)

saveRDS ( MOFAobject_trained , file = file.path (argv$outputDir, "02_MOFAobject_trained_dry.RDS" ))

factors_mofa = MOFAobject_trained@expectations$Z$group1
colnames (factors_mofa) =  paste0 ( "MOFA_" , 1:10)

saveRDS ( factors_mofa , file = file.path (argv$outputDir, "factors_s_mofa.RDS" ))

#Delete temp file if it exists

if (file.exists(outfile)) {
  file.remove(outfile)
}