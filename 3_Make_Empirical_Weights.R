### STEP 3: make empirical weights for each of the positions in the WES files based on the partitioned heritability analysis
# this will assign a weight to each variant in the genome for SKAT/SKAT-O prioritization.

###Load the libraries
library(R.utils)
library(readr)
library(data.table)

###Define directories
phenoName<-"pheno_name"
weightName<-"weight_name"
basedir<-"path_here"
workdir<-paste(basedir,phenoName,"/",sep="")
#setwd(workdir)
weightdir<-paste(workdir,"weights/",sep="")
ldscdir<-paste(workdir,"ldsc/",sep="")

###
mapPre <- "EmpWeights_UKBexome_chr" #this will be the prefix for the name of the output 

###Read the partitioned heritability file from XXX.R
h2File<-paste(ldscdir,weightName,"_partitionedh2Processed.results",sep="")
h2<-fread(file=h2File) 

h2Names<-names(h2) ### Store the column names into an object
head(h2) ###Custom annotations are in rows 1 to 4: PEC_enhancers, PFC_H3K27ac, TC_H3K27ac, and CBC_H3K27ac
dim(h2) ###101 rows and 10 columns. remember that we have additional rows here (such as MAF bins).

h2W <- h2 ### we will manipulate this file and keep h2 clean
h2W<-as.data.frame(h2,stringsAsFactors=F)  ###as dataframe
h2W$Category <- gsub("flanking","extend",h2W$Category) ### in the h2 file, extended regions names are called "flanking". this causes merging issues because
# these "flanking" regions are called "extend" in annotations. therefore, we fix the names here by gsubing "flank" with "extend". this should take care of it.

rownames(h2W) <- h2W$Category #add rownames 
h2Master<-rownames(h2W) #this is the list of all the partitioned heritability annotations. 
#101 annotations (including some that are NOT included in empirical weighting)


### We have some extreme outliers and need to remove annotation with outlier enrichment values
summary(h2W)

h2W<-h2W[-grep("MAF_Adj_Predicted_Allele_AgeL2_4",h2W$Category),] # enrichment -79141.15
h2W<-h2W[-grep("MAF_Adj_LLD_AFRL2_4",h2W$Category),] # enrichment -106.9946
h2W<-h2W[-grep("MAF_Adj_ASMCL2_4",h2W$Category),] #enrichment 1.856911e+13
dim(h2W) #98 rows amd 10 columns
h2Master<-rownames(h2W)# Remake the name of the columns object after filtering the outliers

## Make key for partitioned results output
partKey <- colnames(h2W) #partitioned heritability results columns
partKey <- partKey[2:10] #get the columns used for partitioning (excluding the "category" column)
###Make key of annotations
annKey <- rownames(h2W) #annotation names
annKey <- gsub("L2_4","",annKey,fixed = T) #fix the names for each partitioned category. they have "L2_4" suffix which we are removing
###############################################################################

## we will load the bim files we created in the previous step that define the annotation presence/absence (1/0) for each variant in the WES files 
bimPre<-"UKBexomeOQFE_" ### this is the prefix for the bim files - same as the bim files in the UKB folder
annotDir<-paste(weightdir,"annot/",sep="") #path to the directory for annotation maps for presence absence (1/0) of annotation classes

###Load annotated maps
mapExample<- "example.rds" # set the name for them. for example chrXXX where XXX will be the chromosome number 
chrKey <- 1:22
nChr <- 1:22

# open a loop to prepare the files generated in the previous step for each chromosome for to construct empirical weights for them
for (zz in nChr){ #open the for loop for each chromosome 
  chrT<-chrKey[zz] #starting with chromosome ChrT which is zz out of 1:22
  mapT<-gsub("XXX", chrT, mapExample) #gsub XXX with the chromosome number (1:22 no sex chr)
  inFileRDS<-paste(annotDir,mapT,sep="") #create an object for that  WES chromosome file for example
  load(inFileRDS) ### load it. object should be bimT after it is loaded into R
  #head(bimT)
  #dim(bimT)
  bimKey<-names(bimT)[-c(1:10)] #name of the annotations. we remove the first 10 columns because they represent information about the variants (SNP,chr,bp etc) rather than annotations
  bimNames<-names(bimT)  #names put in an object
  names(bimT) [bimNames %in% "extend"]<-"flanking" #make sure that you fix the "extend" "flanking" discrepancy described above to have uniform named for all annotations
  # for exmaple, WeakEnhancer_Hoffman.flanking.500 in annotation file and WeakEnhancer_Hoffman.extend.500 in bim file
  
  table(annKey %in% bimKey) #check the N for these
  table(bimKey %in% annKey) #check the N for these 
  
  bimKeep<-c(names(bimT)[c(1:10)],bimKey[bimKey %in% annKey]) #86 # get the annotations that we also have partition heritability for them. this excludes the ones that are only included in partitioned heritability analysis like MAF bins 
  bimT <- as.data.frame(bimT) #make it a data frame for further manipulation
  bimW<-bimT[,bimKeep] ### subset the  bim file to include only these
  annKey[!annKey %in% bimKey] #base, GERP.RSsup4, GERP.NS, MAFbin 1 to 10, Nucleotide_Diversity_10kb, Recomb_Rate_10kb, synonymous, non_synonymous, CpG_Content_50kb
  bimKey[!bimKey %in% annKey] #Vahedi_Tcell_SE_500bp, Vahedi_Tcell_SE, Vahedi_Tcell_TE_500bp, Vahedi_Tcell_TE
  
  # another reminder:  there is more in annKey because that includes extra annotations used solely for ldsc analysis such as MAF bins etc baselines
  
  ### Make h2W with only matching annotations
  keepAnn<-annKey[annKey %in% bimKey] #keeping those annotations which are also in bimT 
  h2MasterW<-h2Master[grep(paste(keepAnn,collapse="|"),h2Master)] # re-create a file with the list of these 
  tmp <- h2W[h2W$Category%in%h2MasterW,] #make sure the names are right
  h2MasterW <- tmp
  
  ###Make variable keys for pulling specific fields of interest.
  
  h2_enrich <- h2W[,c(1,5)] #get the enrichment values
  h2_enrich_p <- h2W[,c(1,7)] #get the enrichment p-values
  h2_coeff <- h2W[,c(1,8)] #get the coefficients
  
  length(h2_enrich)==length(h2_enrich_p) #  TRUE
  length(h2_enrich)==length(h2_coeff) #  TRUE
  
  
  t1_enrich <- h2_enrich #rename 
  t1_enrich <- as.data.frame(t1_enrich) #make it a dataframe
  t1_enrich$Category <- gsub("L2_4","",t1_enrich$Category) #clean the name
  names <- t1_enrich$Category #names 
  t1_enrich <- as.data.frame(t1_enrich$Enrichment) #enrichment values
  row.names(t1_enrich) <- names #rownames
  t1_enrich <- as.matrix(t1_enrich) #make it a matrix
  
  
  ###same as above but for enrichment p-values instead of enrichment values.
  t1_enrichP <- as.data.table(h2_enrich_p)
  t1_enrichP <- as.data.frame(t1_enrichP)
  t1_enrichP <- as.data.frame(t1_enrichP)
  t1_enrichP$Category <- gsub("L2_4","",t1_enrichP$Category)
  names <- t1_enrichP$Category
  t1_enrichP <- as.data.frame(t1_enrichP$Enrichment_p)
  row.names(t1_enrichP) <- names
  t1_enrichP <- as.matrix(t1_enrichP)
  
  
  dim(t1_enrichP)==dim(t1_enrich) #TRUE TRUE
  length(t1_enrichP)==length(t1_enrich) #TRUE 
  
  ###Make map 
  mapT<-t(bimW[,-c(1:10)]) # create a new map file by transposing the bim annotations
  mapT<-mapT[order(row.names(mapT)),] #order the file
  
  ###Verify names match after manipulations
  table(row.names(mapT)==row.names(t1_enrich)) #they don't match fixing it below
  ### working to fix the ordering
  t1_enrich <- as.data.frame(t1_enrich[match(rownames(mapT), rownames(t1_enrich)), ])
  t1_enrichP <- as.data.frame(t1_enrichP[match(rownames(mapT), rownames(t1_enrichP)), ])
  t1_enrich <- as.matrix(t1_enrich)
  t1_enrichP <- as.matrix(t1_enrichP)
  
  row.names(t1_enrich)==row.names(mapT) #fixed
  row.names(t1_enrichP)==row.names(mapT) #fixed
  table(row.names(mapT)==row.names(t1_enrich)) #all match now
  
  # we start creating the empirical weights now. we can have different models.
  
  ## model0: add up annotations 0/1 from the basic sets of annotations (no custom) if present
  # model0 includes all the baseline annotation categories WITHOUT the custom annotation categories. 
  # custom annotaions are: PEC_enhancers,CBC_H3K27ac,PFC_H3K27ac,TC_H3K27ac
  tmp <- mapT[-c(8,54,55,64),] #remove these 4
  mapT0 <- tmp # name the file
  model0 <- colSums(mapT0) #get the variant weights
  #table(model0)
  #summary(model0)
  
  ## Model1: Add up annotations 0/1 from the basic sets of annotations AND also add custom annotations if present
  model1<-colSums(mapT) #get the variant weights
  #table(model1)
  #summary(model1)
  #plot(density(as.matrix(mapT)))
  
  ## Model2: Add up annotations (0/1) only if enrichment is significant for that category 
  # we go with p-value of 0.16 (AIC) 
  t1_enrich0 <- t1_enrich
  t1_enrichP0 <- t1_enrichP
  pCut<-0.16/(length(t1_enrichP0)) #is that AIC value
  #table(t1_enrichP0<0.05)
  #keepVecT<-(as.numeric(t1_enrichP0<0.05))  ###This was used on weights generated 3/10/2022
  keepVecT<-(as.numeric(t1_enrichP0<pCut)) #keep the ones that are significant
  mapT_m1<-apply(mapT,2,function(x) {keepVecT*(x)}) #now create mapT_m1  where you only keep the significant ones
  model2<-apply(mapT_m1,2,sum) #get the variant weights for significant ones.
  #table(model2)
  #summary(model2)
  
  
  ###Output for weights before summing
  mapT_transposed <- (t(mapT)) #prep the file for saving
  mapOut2 <- cbind(bimT[,1:10],mapT_transposed) #static columns added 
  ###Output weight for model after summing the models
  mapOut<-cbind(bimT[,1:10],model0,model1,model2) #put them all into a single dataframe
  
  outFileRDS<-paste(outDir,mapPre,chrT,"_","GRCh38","_",chrT,".rds",sep="") #prep the output name of the weight file called mapOut for saving - this will be a separate file for each chromosome.
  save(mapOut, file=outFileRDS) #save it as an R object
  
  
} #close the loop. we shuold now have 1 file for each WES chromosoem file with empirical weights for them using these 3 models. the models can be improved upon. but forn ow we have these basic weights.


#######################
#######################
#######################
