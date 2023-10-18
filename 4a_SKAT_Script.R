### Step 4a: Create script for SKAT
######### SKAT and SKAT-O Script

args <- commandArgs(TRUE) #allow arguments 
mafCut<-args[1] #maf cutoff
estMAF<-args[2] #estMAF method
PhenoName <- args[3] #start range for the intervals to be tested 
WeightName <-  args[4] #end range for the intervals to be tested
GeneName <- args[5] #phenotype 
Chr <- args[6]

t1 <- Sys.time()
library(SKAT) #load skat
library(data.table) #load data.table 
rundate<-gsub("-","",Sys.Date()) #get the rundate

#set the results directory
ResultsDir <- "path_here"
ResultPath <- ResultsDir #result directory for skat tests
setwd(ResultsDir) #set the working directory to that folder

PhenoPath <- "path_here" #path to the phenotype file 
PhenoFile <- "file.pheno" #phenotype file 
CovFile <- "file.cov" #covariate file

#load pheno and cov files
Pheno <- as.data.frame(fread(paste(PhenoPath,PhenoFile,sep=""))) #read in the pheno file
Cov <- as.data.frame(fread(paste(PhenoPath,CovFile,sep=""))) #read in the covariate file

#path to the weights
WeightsDir <- "path_here"
WeightsDir <- gsub("xx",WeightName,WeightsDir)
WeightsPath <- WeightsDir #weight directory for custom weights
WeightRDS <- "EmpWeights_UKBexome_chrxx_GRCh38_xx.rds" #name of the weight files

#load the empirical weights
EmpWeight <- gsub("xx",Chr,WeightRDS) # get the right name for the weight file
load(paste(WeightsPath,EmpWeight,sep="")) #load it into R
mapOut$model3log <- log(mapOut$model3) #log transform model 3
mapOut$model3log[mapOut$model3log<0] <- 0 #minus values should be zero. rarely happens (less than 5 usually, due to transformation)


#load fam file to prepare phenotype and covariate objects
ExomePath <- "path_here" #path to the exomes
Fam<-"file.fam" #fam file. xx will be replaced with the chromosome fed through the argument file

#load the fam file to get the phenotype and covariate files in order and ready for SKAT
FamName <- gsub("xx",Chr,Fam) #get the right fam file
File.Fam <- paste(ExomePath,FamName,sep="") #get the full path to the fam file
FAM<-Read_Plink_FAM(File.Fam, Is.binary=FALSE) #read it using SKAT
phenoName<-PhenoName #get the phenotype 
Pheno<-Pheno[,c("iid",phenoName)] #subset for that phenotype
Pheno<-na.omit(Pheno) #remove NAs
Pheno<-Pheno[order(Pheno$iid),] #order them
y <- Pheno[,phenoName] #put the outcome in a vector format
Cov <- Cov[,-c(1,3,4)] #get the covariates we need 
Cov <- Cov[order(Cov$iid),] #order them 
Cov <- as.matrix(Cov[,-c(1)]) #remove the ids and put it in a matrix format

#path to the SSD files
SSDPath <- "path_here"
SSDFileIn <- paste(SSDPath,"chr",Chr,"/",GeneName,".SSD",sep="")
SSDInfoIn <- paste(SSDPath,"chr",Chr,"/",GeneName,".info",sep="")
length(SSDFileIn)==length(SSDInfoIn) #should be TRUE

obj<-SKAT_Null_Model(y ~ Cov, out_type="C") #run the null + covars
File.SSD <- SSDFileIn #get the full path to that SSD file
File.Info <- SSDInfoIn #get the full path to that info file
SSD.INFO<-Open_SSD(File.SSD, File.Info) #open the SSD file
Z <- Get_Genotypes_SSD(SSD.INFO, 1) #get the genotypes
Z <- as.data.frame(Z) #put it into a dataframe format for manipulation
tmp <- cbind(Z,FAM[,2]) #bind it to the fam file IDs
length1 <- length(colnames(tmp)) #get ids
colnames(tmp)[length1] <- "iid" #call it iids
tmp$iid <- as.numeric(tmp$iid) #make it numeric
tmp <- tmp[order(tmp$iid),] #order it to match the phenotypes
tmp2 <- merge(tmp,Pheno,by="iid") #merge it with the pheno file
length2 <- length(colnames(tmp2)) #get the length 
Z <- tmp2[,-c(1,length2)] #create the new ordered and merged with phenotype genotypes
Z[Z == 9] <- NA #mark missing (9) as NA
Z <- as.matrix(Z) #put into matrix format
Variants <- colnames(Z) #get the variant IDs
WeightFileIn <- subset(mapOut, ID %in% Variants) #merge the weights with the variants.
WeightFileIn <- WeightFileIn[,"model3log"] #get model 3 log transformed

out.skat <- tryCatch(SKAT(Z, obj, max_maf=mafCut, estimate_MAF=estMAF, method="SKAT"), error=function(e) NULL )  #run skat default, catch any errors
out.skato <- tryCatch(SKAT(Z, obj,max_maf=mafCut, estimate_MAF=estMAF, method="SKATO"),error=function(e) NULL ) #run skato default, catch any errors
Custom.out.skat <- tryCatch(SKAT(Z, obj, max_maf=mafCut, estimate_MAF=estMAF, weights=WeightFileIn, method="SKAT"), error=function(e) NULL) #run skat custom, catch any errors
Custom.out.skato <- tryCatch(SKAT(Z, obj,max_maf=mafCut,estimate_MAF=estMAF, weights=WeightFileIn, method="SKATO"), error=function(e) NULL) #run skato custom, catch any errors

tmp.out.skat <- NULL #create tmp files to put the results in there for skat
tmp.out.skato <- NULL #create tmp files to put the results in there for skato
tmp.Custom.out.skat <- NULL #create tmp files to put the results in there for skat custom
tmp.Custom.out.skato <- NULL #create tmp files to put the results in there for skato custom


# run if else to catch the possible errors. if a job was not completed sucessfully due to not having enough markers after maf filtering,
# or cusotm weights all being zero, simply put NA in their place so the job does not fail and you at least get an empty output.  

if(is.null(out.skat)){ print(paste("ERROR! Default weight SKAT analysis of ",GeneName," failed. Did not have enough markers after MAF cutoff",sep=""))}
if(is.null(out.skat)){ tmp.out.skat <- as.data.frame(cbind(tmp.out.skat$Test.Type <- "NA",
                                                           tmp.out.skat$Q <-"NA",
                                                           tmp.out.skat$p.value <- "NA",
                                                           tmp.out.skat$param$liu_pval <- "NA",
                                                           tmp.out.skat$param$Is_Converged <- "NA",
                                                           tmp.out.skat$param$n.marker <- "NA",
                                                           tmp.out.skat$param$n.marker.test <- "NA"
))
} else {tmp.out.skat <- cbind(tmp.out.skat$Test.Type <- as.character(out.skat$Test.Type),
                              tmp.out.skat$Q <- as.numeric(out.skat$Q), #Q value for default skat
                              tmp.out.skat$p.value <- as.numeric(out.skat$p.value), #defaul skat pvalue
                              tmp.out.skat$liu_pval <- as.numeric(out.skat$param$liu_pval), #liu pvalue for defaul skat
                              tmp.out.skat$Is_Converged <- as.numeric(out.skat$param$Is_Converged), #was it converged
                              tmp.out.skat$n.marker <- as.numeric(out.skat$param$n.marker), #number of markers in total in the interval
                              tmp.out.skat$n.marker.test <- as.numeric(out.skat$param$n.marker.test))
}
if(length(tmp.out.skat)<7) {print("length less than 7! all SNPs were filtered out because of MAF filtering")}
if(length(tmp.out.skat)<7) {tmp.out.skat <- NULL}
if(length(tmp.out.skat)<7){tmp.out.skat <- as.data.frame(cbind(tmp.out.skat$Test.Type <- "NA",
                                                               tmp.out.skat$Q <-"NA",
                                                               tmp.out.skat$p.value <- "NA",
                                                               tmp.out.skat$param$liu_pval <- "NA",
                                                               tmp.out.skat$param$Is_Converged <- "NA",
                                                               tmp.out.skat$param$n.marker <- "NA",
                                                               tmp.out.skat$param$n.marker.test <- "NA"
))}	




if(length(out.skato$param$rho_est)>1){print(paste("NOTE! Default weight SKATO analysis of ",GeneName," did not converge!",sep=""))}
if(length(out.skato$param$rho_est)>1){out.skato="NULL"}
if(is.null(out.skato)){ print(paste("ERROR! Default weight SKATO analysis of ",GeneName," failed. Did not have enough markers after MAF cutoff",sep=""))}
if(is.null(out.skato)){ tmp.out.skat <- as.data.frame(cbind(tmp.out.skato$p.value <-"NA",
                                                            tmp.out.skato$param$minp <- "NA",
                                                            tmp.out.skato$param$rho_est <- "NA"
))
} else { tmp.out.skato <- as.data.frame(cbind(tmp.out.skato$p.value <- as.numeric(out.skato$p.value),
                                              tmp.out.skato$param$minp <- as.numeric(out.skato$param$minp),
                                              tmp.out.skato$param$rho_est <- as.numeric(out.skato$param$rho_est)))
}
if(length(tmp.out.skato)<3) {print("length less than 3! all SNPs were filtered out because of MAF filtering")}
if(length(tmp.out.skato)<3) {tmp.out.skato <- NULL}
if(length(tmp.out.skato)<3) {tmp.out.skato <- as.data.frame(cbind(tmp.out.skato$p.value <-"NA",
                                                                  tmp.out.skato$param$minp <- "NA",
                                                                  tmp.out.skato$param$rho_est <- "NA"
))}



if(is.null(Custom.out.skat)){ print(paste("ERROR! Default weight SKAT analysis of ",GeneName," failed. Did not have enough markers after MAF cutoff",sep=""))}
if(is.null(Custom.out.skat)){ tmp.Custom.out.skat <- as.data.frame(cbind(tmp.Custom.out.skat$Test.Type <- "NA",
                                                                         tmp.Custom.out.skat$Q <-"NA",
                                                                         tmp.Custom.out.skat$p.value <- "NA",
                                                                         tmp.Custom.out.skat$param$liu_pval <- "NA",
                                                                         tmp.Custom.out.skat$param$Is_Converged <- "NA",
                                                                         tmp.Custom.out.skat$param$n.marker <- "NA",
                                                                         tmp.Custom.out.skat$param$n.marker.test <- "NA"
))
} else {tmp.Custom.out.skat <- cbind(tmp.Custom.out.skat$Test.Type <- as.character(Custom.out.skat$Test.Type),
                                     tmp.Custom.out.skat$Q <- as.numeric(Custom.out.skat$Q), #Q value for default skat
                                     tmp.Custom.out.skat$p.value <- as.numeric(Custom.out.skat$p.value), #defaul skat pvalue
                                     tmp.Custom.out.skat$liu_pval <- as.numeric(Custom.out.skat$param$liu_pval), #liu pvalue for defaul skat
                                     tmp.Custom.out.skat$Is_Converged <- as.numeric(Custom.out.skat$param$Is_Converged), #was it converged
                                     tmp.Custom.out.skat$n.marker <- as.numeric(Custom.out.skat$param$n.marker), #number of markers in total in the interval
                                     tmp.Custom.out.skat$n.marker.test <- as.numeric(Custom.out.skat$param$n.marker.test))
}
if(length(tmp.Custom.out.skat)<7) {print("length less than 7! all SNPs were filtered out because of MAF filtering")}
if(length(tmp.Custom.out.skat)<7) {tmp.Custom.out.skat <- NULL}
if(length(tmp.Custom.out.skat)<7){tmp.Custom.out.skat <- as.data.frame(cbind(tmp.Custom.out.skat$Test.Type <- "NA",
                                                                             tmp.Custom.out.skat$Q <-"NA",
                                                                             tmp.Custom.out.skat$p.value <- "NA",
                                                                             tmp.Custom.out.skat$param$liu_pval <- "NA",
                                                                             tmp.Custom.out.skat$param$Is_Converged <- "NA",
                                                                             tmp.Custom.out.skat$param$n.marker <- "NA",
                                                                             tmp.Custom.out.skat$param$n.marker.test <- "NA"
))}	



if(length(Custom.out.skato$param$rho_est)>1){print(paste("NOTE! Default weight SKATO analysis of ",GeneName," did not converge!",sep=""))}
if(length(Custom.out.skato$param$rho_est)>1){Custom.out.skato="NULL"}
if(is.null(Custom.out.skato)){ print(paste("ERROR! Default weight SKATO analysis of ",GeneName," failed. Did not have enough markers after MAF cutoff",sep=""))}
if(is.null(Custom.out.skato)){ tmp.Custom.out.skato <- as.data.frame(cbind(tmp.Custom.out.skato$p.value <-"NA",
                                                                           tmp.Custom.out.skato$param$minp <- "NA",
                                                                           tmp.Custom.out.skato$param$rho_est <- "NA"
))
} else { tmp.Custom.out.skato <- as.data.frame(cbind(tmp.Custom.out.skato$p.value <- as.numeric(Custom.out.skato$p.value),
                                                     tmp.Custom.out.skato$param$minp <- as.numeric(Custom.out.skato$param$minp),
                                                     tmp.Custom.out.skato$param$rho_est <- as.numeric(Custom.out.skato$param$rho_est)))
}
if(length(tmp.Custom.out.skato)<3) {print("length less than 3! all SNPs were filtered out because of MAF filtering")}
if(length(tmp.Custom.out.skato)<3) {tmp.Custom.out.skato <- NULL}
if(length(tmp.Custom.out.skato)<3) {tmp.Custom.out.skato <- as.data.frame(cbind(tmp.Custom.out.skato$p.value <-"NA",
                                                                                tmp.Custom.out.skato$param$minp <- "NA",
                                                                                tmp.Custom.out.skato$param$rho_est <- "NA"
))}

colnames(tmp.out.skat) <- c("Default_Test_Type","Default_Q","DefaultW_P.value","Default_liu_pval","Default_Is_Converged","Default_n.marker","Default_n.marker.test")
colnames(tmp.out.skato) <- c("DefaultW_skato_P.value","Default_Skato_minp","Default_skato_rho_est")
colnames(tmp.Custom.out.skat) <- c("Custom_Test.Type","Custom_Q","CustomW_P.value","Custom_liu_pval","Custom_Is_Converged","Custom_n.marker","Custom_nmarker.test")
colnames(tmp.Custom.out.skato) <- c("CustomW_skato_P.value","Custom_skato_minp","Custom_skato_rho_est")

#put everything into a single output file 
outT <- cbind(GeneName,
              phenoName,
              tmp.out.skat,
              tmp.out.skato,
              tmp.Custom.out.skat,
              tmp.Custom.out.skato)

IndSave <- paste("Chr",Chr,"_",GeneName,"_",PhenoName,"_",WeightName,"_model3log","_SKAT_SKATO_",mafCut,"_",rundate,".skat",sep="")
outT <- as.data.frame(outT)
write.table(outT,IndSave,row.names=F,col.names=T,quote=F,sep="\t")
t2 <- Sys.time() #get t2 
t2-t1
Close_SSD() #close the SSD file	

