### STEP 2: annotate the positions in the WES files for presence/absence (1/0) of an annotation categories.
# this will determine whether a particular variant will get a weight or not.

#load the libraries
library(R.utils)
library(readr)
library(data.table)
library(dplyr)
library(GenomicRanges)

# note that V2.2 annotations used for annotating the WES positions DOES NOT include MAF bins because they are only used 
# for partitioning the heritability. therefore, the number of annotation sources is less than the ones that went 
# into partitioning the heritability.

setwd("path_here") #set working directory
workDir<-"path_here" #working dir path to the bed files used for paritioned heritability
enrichDir<-"path_here" #where we will place the annotations

###Where are source annotation files in bed format
###Source https://alkesgroup.broadinstitute.org/LDSCORE/
bedDir<-"path_here" #base annotations + custom annotation bed files (excluding girdhar et al SCZ speicfic annots). they were lifted over to bulid 38 from build 37 using liftOver.
bedVersion<-"baseline_v2" #version of the annotation bed files
bedRaw<-list.files(bedDir,full.names=F, recursive = FALSE) #list all the bed files in the folder
length(bedRaw) #87 annotation
bedRaw
bedSuf<-".bed" #suffix 
bedFiles<-bedRaw[grep(bedSuf,bedRaw,fixed=T)] # assign bedRaw to bedFiles
nBed<-length(bedFiles) #87

genoDir2<-"path_here" #path to UKB variant positions in WES files
bimPre<-"UKBexomeOQFE_chr" ###Example UKBexomeOQFE_chr22.bim etc
bimSuf<-".bim"

outDir <- "path_here" #this is the output directory
chrKey<- 1:22 #list of the chromosomes
chrKey #1 to 22, sex chromosomes excluded
nChr<- 1:22 #24. 22 automosome and 2 sex but sex chromosomes excluded
class(nChr) #integer
nChr

# open a loop to annotate the bed files with ldsc values 
for (i in nChr){ #open the for loop for each chromosome 
  chrT<-chrKey[i] #starting with chromosome ChrT which is zz out of 1:22
  
  genoFile2<-paste(genoDir2,bimPre,chrT,bimSuf,sep="") # exome file to read the bim file for zz out of 1:22
  bimT<-fread(genoFile2) #fread the bim file
  names(bimT)<-c("chr_num","ID","cM","bp","REF", "ALT") #properly name each of the columns
  bimT$chr<-paste("chr",bimT$chr_num,sep="")# add a column where you write the choromosome name such as chr1 chr2 and so on.
  bimT$size<-(nchar(bimT$REF)-nchar(bimT$ALT)) #create a column called size, where you define the size of the base pairs for this position. for exmaple if it is SNP it will be 0, if insertion size is minus if deletion size is more than 1
  bimT$type<-NA #create a column called type which will be populated by the type of variation ie SNP, insertion, or deletion
  bimT$type[bimT$size==0]<-"SNP" #if size is zero, it is SNP polymorphism
  bimT$type[bimT$size<0]<-"INS" #if size is less than 0 or negative, it is insertion INS
  bimT$type[bimT$size>0]<-"DEL" # if the size is more than 0 or positive, it is deletion DEL
  
  bimT$bp_end<-bimT$bp #create a column for the base pair end position. remember this is important for indels. for SNPs it will be the same as bp start 
  bimT$bp_end[bimT$size>0]<-bimT$bp[bimT$size>0]+bimT$size[bimT$size>0] #for indels, if less than 0, or insertion in that location, bp end is the same.  if more than 0, and deletion in that location, bp end is whatever the bp position is + size column for that position
  
  bimW<-GRanges(seqnames=bimT$chr, #create a gRanges object. seqnames=chr number from bim file
                ranges = IRanges(start = bimT$bp, end=bimT$bp_end), # ranges is the strand and end position. remember for SNPs it is the same but it is different for indels 
                strand = "*") #asetrick for strand
  
  ###now that we have prepared the bim file for the variants, we can start to read bed files for annotation through a loop
  for (yy in 1:nBed){ #here we will read yy in 1:nBed with nBed being the total number of bed files 
    
    bedFileT<-bedFiles[yy] #select bedfile yy out of 1:nBed
    bedFullT<-paste(bedDir,bedFileT,sep="") #get the full path and name for that bed file
    
    annNameT<- gsub(".bed","",bedFileT) # name of that annotation file. this will essentially be the name of bed file minus the .bed suffix. as we will see later, it might differ from the original ld scores
    
    t0<-Sys.time() #start time 
    print(paste("Started reading annotation bed file ",yy," of ",nBed," at ",t0,sep="")) #print out that you have started reading the annotation file yy (1:nBeds)
    print(paste("Annotation names is ",annNameT,"",sep="")) #print out the name of the annotation file as created in line 106
    bedT<-fread(bedFullT) #fread the annotation file 
    t1<-Sys.time() # end time after reading
    print(t1-t0) #print the time
    
    if (ncol(bedT)==3) {names(bedT)<-c("chr","bp_start","bp_end")} # if you only have 3 columns, they represent chr , bp start and bp end for that particular annotation
    if (ncol(bedT)==4) {names(bedT)<-c("chr","bp_start","bp_end","value")} # same as line 86 + value
    if (ncol(bedT)==5) {names(bedT)<-c("chr","bp_start","bp_end","value","value2")} #same as line 87 + value2
    
    bedT$binSize<-bedT$bp_end-bedT$bp_start #create a column called binSize. here you define the size of the annotation which is defined as bp-end - bp-start. this varies a lot across different annotations depending on the bin size for that particular annotation source
    bedW<-GRanges(seqnames=bedT$chr, #create a gRanges file for the annotation bed file. seqname is chr
                  ranges = IRanges(start = bedT$bp_start, end=bedT$bp_end), #ranges is bp start and bp end
                  strand = "*") #strand is asetrick
    
    bim_v_bed<-subsetByOverlaps(bimW, bedW) # find the overlap between the bim file which induces the WES variants and bed file for the annotation bed file
    
    annT<-matrix(0,nrow(bimT),1) #create a matrix object with the number of rows based on the WES variants from the bim file
    annT[bimT$bp %in%  as.data.frame(bim_v_bed)$start,]<-1 #if the WES bim file bp position which is the variant, is present in bim_v_bed created in the line above, then give it the value of 1 indicating  presence, otherwise it stays zero indicating absence.
    bimT<-cbind(bimT,annT) # bind annT which defines the presence or absence of annotation for a particular annotation class to the bim file 
    names(bimT)[ncol(bimT)]<-annNameT #name the column for annotation properly based on the object annNameT
    
  }  ###END bed loop. we should have included all the annotation classes and cbind them to the bim file with proper column name for each chromosome
  # which means that we have a file called bimT which includes the annotation of WES variants with their presence or absence (1/0).
  # the file will have 95 columns (the first 10 are static variant information) and N rows based on the number of variants in that bim file
  
  outFileRDS<-paste(outDir,bimPre,chrT,"_",bedVersion,"_GRCh38","_",chrT,".rds",sep="") #prepare the annotation file for saving - this will be a separate file for each chromosome.
  save(bimT, file=outFileRDS) #save it as an R object
 
} ### END Chr loop. by now, we have  saved the annotated WES bim files for each chromosome. 

Sys.time()
