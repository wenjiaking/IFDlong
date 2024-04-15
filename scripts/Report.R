args = commandArgs(trailingOnly=TRUE)
if(!require("rlist")) {install.packages("rlist")}
if(!require("stringr")) {install.packages("stringr")}
if(!require("parallel")) {install.packages("parallel")}
library(rlist)
library(stringr)
library(parallel)

#source("/zfs2/sliu/wenjia/SimulationRO1/scripts/reportFuns.R")
#source("/zfs2/sliu/wenjia/SimulationRO1/scripts/reportFuns_nonparallel.R")

PATH=args[1]
sampleName=args[2]
Aligner=args[3]
buffer=args[4]
anchorLen=args[5]
refEXON=args[6] #"/zfs2/sliu/wenjia/SimulationRO1/refData/allexon_NO.bed"
refAA=args[7] #"/zfs2/sliu/wenjia/SimulationRO1/refData/isoformAA_Ref.txt"
refPseudo=args[8] #"/zfs2/sliu/wenjia/HCC/refData/pseudogenes.rds"
refRoot=args[9] #"/zfs2/sliu/wenjia/HCC/refData/rootNames.txt"
refHMmatch=args[10] #"/zfs2/sliu/wenjia/SimulationRO1/refData/Hm_Mm_match.rds"
species=args[11] #hg38 or mm10
sourceFuns=args[12]
ncores=args[13]

source(sourceFuns)
buffer=as.integer(buffer) #9 by default
anchorLen=as.integer(anchorLen) #10 by default
ncores=as.integer(ncores)

allCDS=read.table(refEXON)
colnames(allCDS)=c("chr","strat","end","name","score","strand")
allCDS$isoform=sapply(strsplit(allCDS$name,"__"),"[[",1)
isoformAA_Ref=read.table(refAA,header = T)


interSbedfile=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS.bed")
intergenebedfile=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_mapped_woSecond_geneTol500intersectS.bed")
coverSReOut=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_Rep.csv")


########isoform annotation##############
match_info=uncoverFilter(interSbedfile,intergenebedfile)
CDScover_source=coverRep.Gen(match_info, CDSref=allCDS,outfile=coverSReOut,tol=buffer,mc=ncores) #mc.cores=30 in parallel.

#####AA annotation###################
#chrN=paste0("chr",c(1:22,"X","Y"))
isoformName=sapply(strsplit(isoformAA_Ref$isoformID,"__"),"[[",1)
isoformAA=isoformAA_Ref
isoformAA$isoformID=isoformName
#isoformAAchrN=isoformAA_Ref[isoformChr %in% chrN,] #39865 isoforms in normal chr
#isoformAAchrN[isoformAAchrN$AAseq=="undivisible" | isoformAAchrN$AAseq=="unknown",c(1,3) ]

###############AAannotation generating functions#########################

coverSReOut=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_Rep.csv")
fullRepOut=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_fullRep.csv")
Rep=read.csv(coverSReOut)
nr=nrow(Rep)
#AArep.list=mclapply(1:nr, function(r) AAannote(Rep[r,],isoformAA),mc.cores=30)
AArep.list=lapply(1:nr, function(r) AAannote(Rep[r,],isoformAA))
AArep=do.call(rbind,AArep.list)
Rep$AAseq=AArep[,1]
Rep$AAnote=AArep[,2]
write.csv(Rep,fullRepOut,row.names = F)



########Filtering#######

###########Pseudo gene reference###############

fullRepOut=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_fullRep.csv")
fusionfiltReOut=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_fusionRep_anchor",anchorLen,"bp.filt.csv")

pseudogenes=ifelse(file.exists(refPseudo), readRDS(refPseudo), NULL)
rootNames=ifelse(file.exists(refRoot), read.table(refRoot,sep = "\t",header = T), NULL)
Hm_Mm_match=ifelse(file.exists(refHMmatch), readRDS(refHMmatch), NULL)

filtRep=filteredRep(reportPath =fullRepOut,fusionfiltPath = fusionfiltReOut,min.len=anchorLen,
                    pseudogenes=pseudogenes,rootNames=rootNames,species = species,Hm_Mm_match=Hm_Mm_match)


