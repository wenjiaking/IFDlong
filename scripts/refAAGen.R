args = commandArgs(trailingOnly=TRUE)
if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if(!require("Biostrings")) {BiocManager::install("Biostrings")}
if(!require("dplyr")) {install.packages("dplyr")}
if(!require("seqRFLP")) {install.packages("seqRFLP")}
if(!require("rtracklayer")) {BiocManager::install("rtracklayer")}
if(!require("rlist")) {install.packages("rlist")}
if(!require("stringr")) {install.packages("stringr")}
if(!require("parallel")) {install.packages("parallel")}

library(dplyr)
library(seqRFLP)
library(parallel)
library(rlist)
library(rtracklayer)
library(Biostrings)

refCDSfa=args[1] #"$codeBase/refData/$ghc/CDS_all.fa"
outisoformfa=args[2] #"$codeBase/refData/$ghc/isoformAll_ref.fa"
outAA=args[3] #"$codeBase/refData/$ghc/isoformAA_Ref.txt"
ncores=args[4]

CDSfa=readDNAStringSet(refCDSfa)
CDSname=names(CDSfa)
CDSseq= paste(CDSfa)
fa_data=data.frame(CDSname, CDSseq)
nameSplit=strsplit(CDSname,"__") #list
transcriptid=sapply(nameSplit,"[[",1)
chr=sapply(nameSplit,"[[",2)
geneName=sapply(nameSplit,"[[",6)
strand=sapply(nameSplit,"[[",5)
Start=as.numeric(sapply(nameSplit,"[[",3))
fa_data$strand=strand
fa_data$transcriptid=transcriptid
fa_data$start=Start
fa_data$geneName=geneName
fa_data$chr=chr
pos_fa=fa_data[strand=="+",] 
neg_fa=fa_data[strand=="-",] 
pos_fa_order=pos_fa[order(pos_fa$start,decreasing=FALSE),]
neg_fa_order=neg_fa[order(neg_fa$start,decreasing=TRUE),]
pos_fa_order[,1]=as.character(pos_fa_order[,1])
pos_fa_order[,2]=as.character(pos_fa_order[,2])
neg_fa_order[,1]=as.character(neg_fa_order[,1])
neg_fa_order[,2]=as.character(neg_fa_order[,2])

isoform_seq=function(order_dat=pos_fa_order,id) {
  isoformSeq=c()
  ind=c()
  gene_name=c()
  start=c()
  chr=c()
  for (i in 1:dim(order_dat)[1]) {
    if (order_dat[i,4]==id) {
      isoformSeq=c(isoformSeq,order_dat[i,2])
      ind=c(ind,i)
      gene_name=c(gene_name,order_dat[i,6])
      start=c(start,order_dat[i,5])
      chr=c(chr,order_dat[i,7])
      print(paste0("chr: ",order_dat[i,7],"; start from ",order_dat[i,5],"; gene: ",order_dat[i,6]))
    }
  }
  return(list(ind,isoformSeq,gene_name,chr,start))
}

transcriptPos_tab=as.data.frame(table(pos_fa_order$transcriptid))
transcriptNeg_tab=as.data.frame(table(neg_fa_order$transcriptid))
library(parallel)
isoformSeq_pos=mclapply(transcriptPos_tab[,1],function(x) isoform_seq(order_dat=pos_fa_order,id=x),mc.cores=ncores)
isoformSeq_neg=mclapply(transcriptNeg_tab[,1],function(x) isoform_seq(order_dat=neg_fa_order,id=x),mc.cores=ncores)

isoform_pos=sapply(isoformSeq_pos,function(x) paste0(x[[2]],collapse = ""))
nCDS_pos=sapply(isoformSeq_pos,function(x) length(x[[2]]))
isoform_pos_gene=sapply(isoformSeq_pos,function(x) unique(x[[3]]))
isoform_pos_chr=sapply(isoformSeq_pos,function(x) unique(x[[4]]))

isoform_neg=sapply(isoformSeq_neg,function(x) paste0(x[[2]],collapse = ""))
nCDS_neg=sapply(isoformSeq_neg,function(x) length(x[[2]]))
isoform_neg_gene=sapply(isoformSeq_neg,function(x) unique(x[[3]]))
isoform_neg_chr=sapply(isoformSeq_neg,function(x) unique(x[[4]]))

isoformName_pos=paste0(transcriptPos_tab[,1],rep("__",length(transcriptPos_tab[,1])),isoform_pos_chr,rep("__",length(transcriptPos_tab[,1])),isoform_pos_gene,
                       rep("__",length(transcriptPos_tab[,1])),nCDS_pos)
isoformName_neg=paste0(transcriptNeg_tab[,1],rep("__",length(transcriptNeg_tab[,1])),isoform_neg_chr,rep("__",length(transcriptNeg_tab[,1])),isoform_neg_gene,
                       rep("__",length(transcriptNeg_tab[,1])),nCDS_neg)
isoformName=c(isoformName_pos,isoformName_neg)
isoformSeq=c(isoform_pos,isoform_neg)
isoformAll=data.frame(isoformName,isoformSeq)

library(seqRFLP)
dataframe2fas(isoformAll,outisoformfa)

############ Generate AA sequence reference ########################
library(Biostrings)
isoformAll=readDNAStringSet(outisoformfa)
isoformName=names(isoformAll)
isoformSeq= paste(isoformAll)
weirdisofID=unique(c(which(nchar(isoformSeq)%%3!=0),grep("NNN",isoformSeq))) #122
isoformAll.adj=isoformAll[-weirdisofID]
undivisibleIND=which(nchar(isoformSeq)%%3!=0 & !grepl("NNN",isoformSeq))
isoformAll.adj=c(isoformAll.adj,
                 substr(isoformAll[undivisibleIND],1,floor(nchar(isoformAll[which(nchar(isoformSeq)%%3!=0)])/3)*3))

AAall=translate(isoformAll.adj,no.init.codon=T,if.fuzzy.codon="X")
AAseq=paste(AAall)
AAisoformName=names(AAall)
# 29917 isoforms in total including 120 undividible isoforms but remove 2 "NNN" isoforms ("NM_028518__chr2__Col20a1(+)__35" "NM_010173__chr4__Faah(-)__15")
isoformAA_Ref=data.frame(isoformID=c(isoformName[-weirdisofID],
                                     isoformName[undivisibleIND]),
                         DNAseq=c(isoformSeq[-weirdisofID],
                                  isoformSeq[undivisibleIND]),
                         AAseq=AAseq,note=c(rep("divisible",(length(isoformName)-length(weirdisofID))),
                                            rep("undivisible",length(undivisibleIND))))
write.table(isoformAA_Ref, file=outAA,
            quote=F, sep="\t", row.names=F,col.names=T)

