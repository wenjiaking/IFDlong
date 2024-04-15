args = commandArgs(trailingOnly=TRUE)

if(!require("rlist")) {install.packages("rlist")}
if(!require("stringr")) {install.packages("stringr")}
if(!require("parallel")) {install.packages("parallel")}
library(rlist)
library(stringr)
library(parallel)

PATH=args[1]
sampleName=args[2]
Aligner=args[3]
buffer=args[4]
anchorLen=args[5]
refGTF=args[6] #"/zfs2/sliu/wenjia/short_seq_prac_wenjia/Collaboration_Wenjia_Silvia/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
sourceFuns=args[7]
ncores=args[8]

buffer=as.integer(buffer) #9 by default
anchorLen=as.integer(anchorLen) #10 by default
ncores=as.integer(ncores)

reportfile=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_Rep.csv")
Isof.quantfile=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_Isof_quant.csv")
Isof.quantRData=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_buffer",buffer,"bp_temp_Isof_quant.RData")
fusionfiltfile=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_fusionRep_anchor",anchorLen,"bp.filt.csv")
fusion.quantfile=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_Fusion_quant_anchor",anchorLen,"bp.csv")
fusion.quantRData=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_buffer",buffer,"bp_temp_Fusion_quant_anchor",anchorLen,"bp.RData")

############
Report=read.csv(reportfile)
fusionRep=read.csv(fusionfiltfile)
gtf.dat=read.table(refGTF,fill=TRUE)

source(sourceFuns)
#source("/zfs2/sliu/wenjia/SimulationRO1/scripts/reportFuns_nonparallel.R")
#source("/zfs2/sliu/wenjia/SimulationRO1/scripts/reportFuns.R")

quantList=quantList.gen(Report,fusionRep,Isof.quantfile,Isof.quantRData,fusion.quantfile,fusion.quantRData,gtf.dat,tol=1e-5,max.iter=200,mc=ncores)





