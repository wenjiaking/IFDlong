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
refEXON=args[4] #"/zfs2/sliu/wenjia/SimulationRO1/refData/allexon_NO.bed"
sourceFuns=args[5]

source(sourceFuns)
allCDS=read.table(refEXON)
colnames(allCDS)=c("chr","strat","end","name","score","strand")
allCDS$isoform=sapply(strsplit(allCDS$name,"__"),"[[",1)

interSbedfile=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS.bed")
bedfile=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_mapped_woSecond.bed")
uncoverSOut=paste0(PATH,"/",sampleName,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_EXONuncover.bed",sep = "")
uncoverReads=uncoverbed.gen(interSbedfile,bedfile,uncoverSOut)

