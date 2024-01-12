library (IntNMF)
library (argparser)

args = arg_parser("Input Arguments")
args = add_argument( args, "--input", help = "Path to RDS input containing the multiOmic Data")
args = add_argument( args, "--nfact", help = "Number of Factors" , type = "integer")
args = add_argument( args, "--outputDir", help = "Path to output result")

argv = parse_args(args)

omics = readRDS (argv$input )
num_factors = argv$nfact

### IntNMF

omics_pos<-list()
omics_pos<-list()
  for(j in 1:length(omics[1:length(omics)])){
    if(min(omics[[j]])<0){
      omics_pos[[j]]<-omics[[j]]+abs(min(omics[[j]]))
    }else{
      omics_pos[[j]]<-omics[[j]]
    }
    omics_pos[[j]]<-omics_pos[[j]]/max(omics_pos[[j]])
  }



factorizations_intnmf = nmf.mnnals(dat = lapply (omics_pos, as.matrix ) ,  k  = num_factors)


saveRDS ( factorizations_intnmf , file = file.path (argv$outputDir, "03_factorizations_intnmf_dry.RDS" ))
