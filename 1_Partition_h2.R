### STEP 1: parition the heritability of your sumstats into functional categories using LDSC + including custom psychENCODE annotation categories

### The 1KG files used in the frqfile-chr and w-ld-chr commands were downloaded from:
### https://alkesgroup.broadinstitute.org/LDSCORE/

######### RUN IN TERMINAL ##############

# <activate LDSC>

<ldsc.py> \ #call LDSC
--h2 <GSCAN_DPW.sumstats.gz> \ #call the munnged sumstats for GSCAN-DPW
--ref-ld-chr <PEC_enhancers/PEC_enhancers.>,\ #path to psychENCODE enhancers annotation ld scores
<PFC/PFC_H3K27ac.>,\ #path to psychENCODE pre frontal cortex h3k27ac annotation ld scores
<TC/TC_H3K27ac.>,\ #path to  psychENCODE temporal cortex h3k27ac annotation ld scores
<CBC/CBC_H3K27ac.>,\ #path to psychENCODE cerebral cortex h3k27ac annotation ld scores
<1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.> \ # path to baseline ldscores
--out <GSCAN_DPW_partitionedh2> \ #output
--overlap-annot  \ #allow overlap
--frqfile-chr <1000G_Phase3_frq/1000G.EUR.QC.> \ #frequency
--w-ld-chr <1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.> \ #weights
--print-coefficients #get the coefficients
# note that V2.2 of LDSC annotations used for partitioning the heritability here, also includes some annotation sources 
# such as MAF bins that are not included in variant empirical weights.
# these categories are recommended to be included in paritioning the heritability using LDSC. 
