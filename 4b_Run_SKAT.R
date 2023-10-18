### STEP 4b: Run SKAT

workdir <- "path_here"
outdir<-paste(workdir,"out/",sep="")
tempdir<-paste(workdir,"temp/",sep="")
sh_tempdir<-paste(tempdir,"sh/",sep="")
SSD_Dir <- "path_here"
logdir<-paste(workdir,"log/",sep="")

rundate<-gsub("-","",Sys.Date())
userID<-"ahangarim"

###Where is R script to run SKAT
Rscript<-paste(workdir,"4a_SKAT_Script.R",sep="")
###This R scripts requires 6 arguments
#mafCut <- args[1]     ###Max minor allele frequency
#estMAF <- args[2]     ###Which method to estimate MAF
#geneNameW <- args[3]  ###Gene / Interval Name
#weightName<-args[4]   ###Gene / Intervals can have multiple weight files
#phenoName<-args[5]    ###Names of phenotype being testing. Need for pulling from phenofile in the pheno directory
#Chr <-args[6]         ### chromosome for the interval/gene


###Define for arguments
phenoDir<-"path_here"
phenoFile<-"file.pheno"
PhenoName<-"pheno_name"
mafCut<-0.01
estMAF<-2

WeightName<-"weight_name"

#get the intervals
SSD_FileList  <- list.files(SSD_Dir,recursive = TRUE) #list of the SSD files
SSDs <- SSD_FileList[grepl(".SSD",SSD_FileList)] #get the SSD lists 
SSDs <- SSDs[!grepl(".txt",SSDs)] #remove files with .txt suffix in that folder
SSDs <- SSDs[!grepl("MakeSSD.R",SSDs)] #remove MakeSSD.R in that folder
SSDs <- SSDs[!grepl("MakeSSD.log",SSDs)] #remove the log file in that folder 
SSDs <- SSDs[!grepl("submit_SSD.sh",SSDs)] #remove the bash submission file in that folder
SSDsLength <- length(SSDs) #length of the intervals testing.
nGenes<-SSDsLength
#nGenes <- 15

##Pick Queue
QueueList<-c("queue1", "queue2")
QueueName<- QueueList[1]

###How many threads
nnodes<-1
nthreads<-2
JobName<-"R_SKAT"

###What is maximum number of jobs
maxAll<-35

##################################
###Start ii loop to submit jobs
for (ii in 1:nGenes){
  geneNameW<-SSDs[ii]
  Chr <- gsub("/.*","",geneNameW)
  Chr <- gsub("chr","",Chr)
  geneNameW <- gsub(".*/","",geneNameW)
  geneNameW
  geneNameW<- gsub(".SSD","",geneNameW)
  weightNameW<-WeightName
  
  nJobsAll<-as.numeric(system(
    paste("qstat | grep ",userID," | wc -l",sep=""),
    intern=TRUE))
  while(nJobsAll>maxAll) {
    print(Sys.time())
    print(nJobsAll)
    Sys.sleep(10)
    nJobsAll<-as.numeric(system(
      paste("qstat | grep ",userID," | wc -l",sep=""),
      intern=TRUE))
  }
  
  ###Create header for shell script 
  startSH<- paste(
    "#!/bin/bash \
#PBS -q ",QueueName," \
#PBS -N ",JobName," \
### Pass all current environment variables to the job \
#PBS -V \
#PBS -S /bin/bash \
###Path for output and error \
###Request resources ppn=processors /node ncpus=number of cpus \
#PBS -l nodes=",nnodes,":ppn=",nthreads," \
###Where to direct output \
#PBS -o ",logdir,JobName,"_Chr",Chr,"_",geneNameW,"_",PhenoName,"_",weightNameW,"_model3log_weight","_",rundate,".log \
#PBS -e ",logdir,JobName,"_Chr",Chr,"_",geneNameW,"_",PhenoName,"_",weightNameW,"_model3log_weight","_",rundate,".err \
 \
cd $PBS_O_WORKDIR \
date \
", sep="")
  
  ### 4a_SKAT_Script.R requires 10 arguments
  RCom1<-paste("Rscript",Rscript,
               mafCut,
               estMAF,
               PhenoName,
               WeightName,
               geneNameW,
               Chr,
               sep=" \\\n")
  
  comAll<-paste(startSH,RCom1,"date",sep="\n")
  #shellName<-paste(tempdir,JobName,"_",geneNameW,"_",rundate,".sh",sep="")
  shellName<-paste(sh_tempdir,JobName,"_",geneNameW,"_",PhenoName,"_",weightNameW,"_",mafCut,"_",rundate,".sh",sep="")
  write(comAll, file=shellName)
  
  ###Submit shell script
  system(paste ("qsub ",shellName, sep=""))
  print(Sys.time())
  print(ii)
  Sys.sleep(1)
}
