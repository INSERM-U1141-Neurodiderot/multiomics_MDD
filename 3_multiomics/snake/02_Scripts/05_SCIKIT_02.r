library (argparser)

args = arg_parser("Input Arguments")
args = add_argument( args, "--nfact", help = "Number of Factors" , type = "integer")
args = add_argument( args, "--outputDir", help = "Path to output result")

argv = parse_args(args)

num_factors = argv$nfact


temp.folder <- paste0(argv$outputDir ,'/scikit_MDD/')

metagenes_scikit = list()  
factors_scikit<-as.matrix(read.table(paste(temp.folder,"signals.txt",sep=""),sep="\t",header=F))
  colnames(factors_scikit)<-1:num.factors
  rownames(factors_scikit)<-colnames(omics_[[1]])
  metagenes_scikit<-list()
  for(j in 1:(length(omics_))){
    metagenes_scikit[[j]]<-as.matrix(read.table(paste(temp.folder,"projomics",j,".txt",sep=""),sep="\t",header=F))
    rownames(metagenes_scikit[[j]])<-rownames(omics_[[j]])
    colnames(metagenes_scikit[[j]])<-1:num.factors
  }
  factorizations_skit<-list(factors_scikit,metagenes_scikit)

factors_scikit<-as.matrix(read.table(paste(temp.folder,"signals.txt",sep=""),sep="\t",header=F))


saveRDS ( factors_scikit , file = file.path (argv$outputDir, "05_factorizations_scikit_dry.RDS" ))
