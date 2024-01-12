library (argparser)

args = arg_parser("Input Arguments")
args = add_argument( args, "--input", help = "Path to RDS input containing the multiOmic Data")
args = add_argument( args, "--nfact", help = "Number of Factors" , type = "integer")
args = add_argument( args, "--outputDir", help = "Path to output result")

argv = parse_args(args)

omics = readRDS (argv$input )
num_factors = argv$nfact


temp.folder <- paste0(argv$outputDir ,'/scikit_MDD/')
dir.create( temp.folder, showWarnings = FALSE)


omics_ = lapply (omics , t)
  files<-""
  for(j in 1:length(omics_)){
    write.table(omics_[[j]],paste(temp.folder,"/omics",j,".txt",sep=""),sep=",",col.names=T, row.names=T)
    files<-paste(files,paste("omics",j,".txt",sep=""),sep=" ")
  }