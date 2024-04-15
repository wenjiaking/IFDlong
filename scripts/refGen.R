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

refGTF=args[1] #"$codeBase/refData/$ghc/genes.gtf"
outEXON=args[2] #"$codeBase/refData/$ghc/allexon_NO.bed"
outTol=args[3]
outCDS=args[4]


####convert genes.gtf reference file into txt file of blocks exon/CDS/start_codon/stop_codon#########
gtf <- rtracklayer::import(refGTF,format="gtf")
gtf_dat=as.data.frame(gtf)
colnames(gtf_dat)[1]="chr"
genes_annote=data.frame(gtf_dat[,c("chr","type","start","end","strand","score","gene_id","transcript_id")])
genes_annote$start=genes_annote$start-1

########Generate exon reference blocks####################
gtf_exon=genes_annote[genes_annote$type=="exon",]
exon_pos=gtf_exon[gtf_exon$strand=="+",]
exon_neg=gtf_exon[gtf_exon$strand=="-",]
exon_pos=exon_pos[order(exon_pos$chr,exon_pos$transcript_id, exon_pos$start, decreasing = FALSE),]
exon_neg=exon_neg[order(exon_neg$chr,exon_neg$transcript_id, exon_neg$start,decreasing = TRUE),]

allexon_feature=rbind(exon_pos, exon_neg)
allexon_NO=allexon_feature%>% group_by(transcript_id)%>% mutate(order=1:length(chr))
nexons=nrow(allexon_NO)
allexon_NO$name=paste0(allexon_NO$transcript_id, rep("__",nexons), allexon_NO$chr,rep("__",nexons), allexon_NO$start,rep("__",nexons),
                       allexon_NO$end,rep("__",nexons),allexon_NO$strand,rep("__",nexons), allexon_NO$gene_id,rep("__",nexons), allexon_NO$order)
write.table(allexon_NO[,c(1,3,4,10,6,5)],outEXON, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

######### generate gene range reference files #########
uniGenes=unique(genes_annote$gene_id)
nchrs=sapply(uniGenes, function(x) length(unique(genes_annote$chr[genes_annote$gene_id==x])))
mchrGene=uniGenes[nchrs==1]
genes_annote_F=genes_annote[which(genes_annote$gene_id %in% mchrGene),]

#######Generate the genes dictionary in bed format after filtering out those multiple-across chromosone genes########
gene_range=function(geneData) {
  geneData=geneData[order(geneData$start,decreasing = FALSE),]
  n_block=dim(geneData)[1]
  gene_info=c(geneData$chr[1],geneData$start[1],geneData$end[n_block],geneData$strand[1])
  return(gene_info)
}

gene_position_F=lapply(mchrGene,function(x) gene_range(genes_annote_F[genes_annote_F$gene_id==x,]))
genes_bed_F=data.frame(chr=sapply(gene_position_F,"[[",1),start=as.numeric(sapply(gene_position_F,"[[",2)),end=as.numeric(sapply(gene_position_F,"[[",3)),
                       strand=sapply(gene_position_F,"[[",4),gene_id=mchrGene)
ngenes=nrow(genes_bed_F)
genes_bed_F$name=paste0(genes_bed_F$gene_id,rep("__",ngenes),genes_bed_F$chr,rep("__",ngenes),genes_bed_F$start,rep("__",ngenes),genes_bed_F$end,rep("__",ngenes),genes_bed_F$strand)

write.table(data.frame(genes_bed_F$chr,sapply(genes_bed_F$start-500,function(x) max(0,x)),genes_bed_F$end+500,genes_bed_F$name,score=rep(NA,ngenes),genes_bed_F$strand),
            file=outTol,
            sep="\t",quote = FALSE,col.names = FALSE,row.names = FALSE) #tolerate +- 500 base error

############## Generate AA sequence for known isoforms ########
gtf_cds=gtf_dat[gtf_dat$type=="CDS",]
gtf_cds$start=gtf_cds$start-1
nCDS=nrow(gtf_cds)
Name=paste0(gtf_cds$transcript_id,rep("__",nCDS),gtf_cds$chr,rep("__",nCDS),gtf_cds$start,rep("__",nCDS),gtf_cds$end,rep("__",nCDS),
            gtf_cds$strand,rep("__",nCDS),gtf_cds$gene_name)

gtf_cds$Name=Name
gtf_bed=gtf_cds[,c(1:3,15,8,5)]
write.table(gtf_bed, file=outCDS, quote=F, sep="\t", row.names=F,col.names=F)
#run bash command: bedtools getfasta -fi /zfs2/sliu/reference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -bed /zfs2/sliu/wenjia/SimulationRO1/refData/MM10CDS_all.bed -split -s -name -fo /zfs2/sliu/wenjia/SimulationRO1/refData/MM10CDS_all.fa
