#################################################################################################
#####filter out reads with at least one block not 90% covered and save them in bed format#######
#################################################################################################
uncoverbed.gen=function(interbedfile,bedfile,outfile) {
  intersect_tab=read.table(interbedfile)
  colnames(intersect_tab)=c("chr","start","end","SampleID","score","strand","CDS_chr","CDS_start","CDS_end","CDS_name","CDS_score","CDS_strand","n_base")
  uncoverblocks=paste(intersect_tab$chr[intersect_tab$CDS_name=="."],intersect_tab$start[intersect_tab$CDS_name=="."],intersect_tab$end[intersect_tab$CDS_name=="."],
                      intersect_tab$SampleID[intersect_tab$CDS_name=="."],sep=":")
  bedPos_tab=read.table(bedfile)
  uncoverbed=bedPos_tab[match(uncoverblocks,paste(bedPos_tab$V1,bedPos_tab$V2,bedPos_tab$V3,bedPos_tab$V4,sep=":")),]
  #unmatched_reads=bedPos_tab[(bedPos_tab[,4]%in% unmatched_sampleID),]
  write.table(uncoverbed, outfile,sep="\t",col.names = FALSE, row.names = FALSE,quote = FALSE)
  print("save the uncovered blocks in bed file")
}

uncoverFilter=function(interbedfile,intergenebedfile) {
  intersect_tab=read.table(interbedfile)
  colnames(intersect_tab)=c("chr","start","end","SampleID","score","strand","CDS_chr","CDS_start","CDS_end","CDS_name","CDS_score","CDS_strand","n_base")
  match_intersect=intersect_tab[intersect_tab$CDS_name!=".",]
  CDS_info=strsplit(match_intersect$CDS_name, "__")
  match_info=match_intersect[,c(1:4,6:9,12,13)]
  match_info$gene=sapply(CDS_info, "[[",6)
  match_info$isoform=sapply(CDS_info,"[[",1)
  match_info$order=sapply(CDS_info,"[[",7)
  if (class(try(read.table(intergenebedfile),silent = T))=="try-error") {
    read_info=as.data.frame(match_info)
  }
  else {
    uncover_info=read.table(intergenebedfile)
    colnames(uncover_info)=c("chr","start","end","SampleID","score","strand","gene_chr","gene_start","gene_end","gene_name","gene_score","gene_strand","n_base")
    #Classify reads with all blocks 90% covered by EXON and save the report in txt file
    unmatch_info=uncover_info[,c(1:4,6:9,12,13)]
    unmatch_info$gene="undefined"
    if (sum(uncover_info$gene_name!=".")>0) {
      unmatch_info$gene[uncover_info$gene_name!="."]=sapply(strsplit(uncover_info$gene_name[uncover_info$gene_name!="."],"__"),"[[",1)
    }
    unmatch_info$isoform="undefined"
    unmatch_info$order="undefined"
    colnames(unmatch_info)=c("chr","start","end","SampleID","strand","CDS_chr","CDS_start","CDS_end","CDS_strand","n_base","gene","isoform","order")
    read_info=as.data.frame(rbind(match_info,unmatch_info))
  }
  read_info$CDS_strand[read_info$CDS_strand=="."]="undefined"
  read_info$CDS_chr[read_info$CDS_chr=="."]="undefined"
  print("Extract the CDS Covered Alignments Done!")
  return(read_info)
}
#####################Functions#######################
testContinue=function(vec) {
  numbers=as.numeric(vec)
  result <- diff(numbers)
  return(!any(abs(result)!=1))
}

isoformFilter=function(isoformData,tol=9) {
  nCDS=dim(isoformData)[1]
  strand=unique(isoformData$strand)
  if (strand=="-") isoformData=isoformData[order(isoformData$start,decreasing=T),]
  if (nCDS==1) {
    temp.qual=((isoformData$start+tol>=isoformData$CDS_start) & (isoformData$end-tol<=isoformData$CDS_end))
  }
  else {
    if (nCDS==2) {
      if (strand=="+") temp.qual=((isoformData$start[1]+tol>=isoformData$CDS_start[1]) & (abs(isoformData$end[1]-isoformData$CDS_end[1])<=tol) & 
                                    (abs(isoformData$start[2]-isoformData$CDS_start[2])<=tol) &(isoformData$end[2]-tol<=isoformData$CDS_end[2]))
      else temp.qual=((isoformData$end[1]-tol<=isoformData$CDS_end[1]) & (abs(isoformData$start[1]-isoformData$CDS_start[1])<=tol) &
                        (abs(isoformData$end[2]-isoformData$CDS_end[2])<=tol) &(isoformData$start[2]+tol>=isoformData$CDS_start[2]))
    }
    else {
      if (strand=="+") temp.qual=((isoformData$start[1]+tol>=isoformData$CDS_start[1]) & (abs(isoformData$end[1]-isoformData$CDS_end[1])<=tol) &
                                    (abs(isoformData$start[nCDS]-isoformData$CDS_start[nCDS])<=tol) &(isoformData$end[nCDS]-tol<=isoformData$CDS_end[nCDS]) &
                                    all((abs(isoformData$start[2:(nCDS-1)]-isoformData$CDS_start[2:(nCDS-1)])<=tol)) &
                                    all((abs(isoformData$end[2:(nCDS-1)]-isoformData$CDS_end[2:(nCDS-1)])<=tol)))
      else temp.qual=((isoformData$end[1]-tol<=isoformData$CDS_end[1]) & (abs(isoformData$start[1]-isoformData$CDS_start[1])<=tol) &
                        (abs(isoformData$end[nCDS]-isoformData$CDS_end[nCDS])<=tol) &(isoformData$start[nCDS]+tol>=isoformData$CDS_start[nCDS]) &
                        all((abs(isoformData$start[2:(nCDS-1)]-isoformData$CDS_start[2:(nCDS-1)])<=tol)) &
                        all((abs(isoformData$end[2:(nCDS-1)]-isoformData$CDS_end[2:(nCDS-1)])<=tol)))
      
    }
  }
  return(temp.qual)
}

perGene=function(geneData,allCDS,tol=9) {
  temp.tab=as.matrix(table(apply(geneData[,c("chr","start","end","strand")],1,function(x) paste0(x,collapse = ":")) ,geneData$isoform))
  nblock=dim(temp.tab)[1]
  temp.ind=which(apply(temp.tab,2,function(x) sum(x==0))==0) #what if there is no isoform cover all blocks? e.g. reference isoform missed CDS and became novel
  temp.isoform=colnames(temp.tab)[temp.ind]
  temp.strand=sapply(temp.isoform, function(x) length(unique(geneData$strand[geneData$isoform==x])))
  if (any(temp.strand==1)) {
    temp.ind=temp.ind[temp.strand==1]
  }
  else {
    temp.ind=NULL
  }
  if (length(temp.ind)==0) {
    print("There is no isoform which can cover all blocks of the read!")
    #temp.position=NA ###should be redefined later!
    temp.position=paste0(rownames(temp.tab),collapse = ";")
    temp.note="no full-covered isoform"
    qual.source="undefined"
    nCDS_ref=NA
    length_ref=NA
    NO.CDS=NA
    temp.type="novel with addition"
  }
  
  else {
    temp.isoform=colnames(temp.tab)[temp.ind]
    if (all(temp.isoform=="undefined")) {
      print("The read cannot be covered by EXONs!")
      temp.position=paste0(rownames(temp.tab),collapse = ";")
      temp.note="no full-covered isoform"
      qual.source="undefined"
      nCDS_ref=NA
      length_ref=NA
      NO.CDS=NA
      temp.type="novel with undefined"
    }
    else {
      #the possible problem of the apply is that there might only one block in geneData i.e. 
      #geneData[geneData$isoform==temp.isoform[1],c("chr","start","end","strand")] may not be a matrix
      temp.position=apply(geneData[geneData$isoform==temp.isoform[1],c("chr","start","end","strand")],1,function(x) paste0(x,collapse = ":"))
      temp.position=paste0(temp.position,collapse = ";")
      #print(temp.position)
      temp.continue=sapply(temp.isoform, function(x) testContinue(geneData$order[geneData$isoform==x]))
      temp.source=temp.isoform[temp.continue]
      #print(unique(geneData$SampleID))
      #print(temp.source)
      if (length(temp.source)==0) {
        temp.note="discontinuous CDS"
        temp.type="novel with deletion"
        #qual.source="undefined"
        #nCDS_ref=NA
        #length_ref=NA
        #NO.CDS=NA
        qual.source=temp.isoform
        nCDS_ref=sapply(qual.source,function(x) sum(allCDS$isoform==x))
        length_ref=sapply(qual.source,function(x) sum(allCDS$end[allCDS$isoform==x]-allCDS$strat[allCDS$isoform==x]))
        NO.CDS=sapply(qual.source,function(x) paste0(geneData$order[geneData$isoform==x],collapse = "-"))
      }
      else {
        temp.qual=sapply(temp.source,function(x) isoformFilter(geneData[geneData$isoform==x,],tol = tol))
        #print(temp.qual)
        if (any(temp.qual)) {
          qual.source=temp.source[temp.qual]
          temp.note="continuous CDS and edge-matching"
          temp.type="normal"
          nCDS_ref=sapply(qual.source,function(x) sum(allCDS$isoform==x))
          #print(nCDS_ref)
          length_ref=sapply(qual.source,function(x) sum(allCDS$end[allCDS$isoform==x]-allCDS$strat[allCDS$isoform==x]))
          NO.CDS=sapply(qual.source,function(x) paste0(geneData$order[geneData$isoform==x],collapse = "-"))
          #print(length_ref)
        }
        else {
          #qual.source="undefined"
          qual.source=temp.source
          temp.note="continuous CDS but edge-unmatching"
          temp.type="novel with deletion"
          #nCDS_ref=NA
          nCDS_ref=sapply(qual.source,function(x) sum(allCDS$isoform==x))
          #length_ref=NA
          #NO.CDS=NA
          length_ref=sapply(qual.source,function(x) sum(allCDS$end[allCDS$isoform==x]-allCDS$strat[allCDS$isoform==x]))
          NO.CDS=sapply(qual.source,function(x) paste0(geneData$order[geneData$isoform==x],collapse = "-"))
        }
        
      }
    }
    
  }
  
  #print(qual.source)
  return(c(paste0(qual.source,collapse = "||"),temp.position, nblock,paste0(nCDS_ref,collapse = "||"),paste0(length_ref,collapse = "||"),temp.note,paste0(NO.CDS,collapse = "||"),temp.type))
}

filtGene=function(temp.info) {
  blockname=paste(temp.info$chr,temp.info$start, temp.info$end,temp.info$strand,sep=":")
  gene_tab=table(blockname,temp.info$gene)
  gene_tab[gene_tab>1]=1 # In case one gene corresponds to multiple isoforms where the block will be aligned more than once
  
  totalGeneNum=dim(gene_tab)[2]
  totalGene=colnames(gene_tab)
  geneNum=0
  out=list()
  overlapping="N"
  while(geneNum < totalGeneNum & length(out)==0){
    geneNum=geneNum+1	
    indexMatrix = combn(totalGeneNum,geneNum)
    #print(indexMatrix)
    for(j in 1:ncol(indexMatrix)){
      
      if(all(rowSums(as.matrix(gene_tab[,indexMatrix[,j]]))==1)){
        #outgene=paste(totalGene[indexMatrix[,j]],collapse = "&")
        out=c(out,list(totalGene[indexMatrix[,j]]))  # replace
      }
    }
    #print(geneNum)
  }
  if (length(out)==0) {
    geneNum=0
    overlapping="Y"
    while(geneNum < totalGeneNum & length(out)==0){
      geneNum=geneNum+1	
      indexMatrix = combn(totalGeneNum,geneNum)
      #print(indexMatrix)
      for(j in 1:ncol(indexMatrix)){
        
        if(all(rowSums(as.matrix(gene_tab[,indexMatrix[,j]]))>=1)){
          #outgene=paste(totalGene[indexMatrix[,j]],collapse = "&")
          out=c(out,list(totalGene[indexMatrix[,j]]))  # replace
        }
      }
      #print(geneNum)
    }
  }
  return(list(outGene=out,geneNum=geneNum,nblock=dim(gene_tab)[1],position=paste0(rownames(gene_tab),collapse = ";"),overlapping=overlapping))
}



sample_source=function(temp.info,allCDS,tol=9) {
  #temp.info=match_info[match_info$SampleID==sampleID,]
  filtOut=filtGene(temp.info)
  geneNum=filtOut$geneNum
  outGene=filtOut$outGene
  overlapping=filtOut$overlapping
  if (length(outGene)==0) {
    print("There is no full-covered gene or fused genes!")
    #temp.position=NA ###should be redefined later!
    temp.position=filtOut$position
    temp.note="neither full-covered gene nor fused genes"
    qual.source="undefined"
    nCDS_ref=NA
    length_ref=NA
    NO.CDS=NA
    temp.type="novel"
    temp.label=c(paste0(qual.source,collapse = "||"),temp.position, filtOut$nblock,
                 paste0(nCDS_ref,collapse = "||"),paste0(length_ref,collapse = "||"),temp.note,paste0(NO.CDS,collapse = "||"),temp.type,NA,NA,NA)
  } else {
    if (overlapping=="Y") temp.fusion="fusion with overlapping"
    else {
      if (geneNum==1) temp.fusion="N"
      else temp.fusion="Y"
    }
    
    if (geneNum==1) {
      #temp.fusion="N"
      temp.strand=sapply(outGene,function(x) paste0(unique(temp.info$CDS_strand[temp.info$gene==x]),collapse = ""))
      temp.geneData=lapply(outGene,function(x) temp.info[temp.info$gene==x,])
      temp.source=lapply(temp.geneData, function(x) perGene(x,allCDS=allCDS,tol = tol))
      temp.label=c(paste0(sapply(temp.source,"[[",1),collapse = "#"), paste0(sapply(temp.source,"[[",2),collapse = "#"),
                   paste0(sapply(temp.source,"[[",3),collapse = "#"),paste0(sapply(temp.source,"[[",4),collapse = "#"),
                   paste0(sapply(temp.source,"[[",5),collapse = "#"),paste0(sapply(temp.source,"[[",6),collapse = "#"),
                   paste0(sapply(temp.source,"[[",7),collapse = "#"),paste0(sapply(temp.source,"[[",8),collapse = "#"),
                   paste0(unlist(outGene),collapse = "#"),temp.fusion,paste0(temp.strand,collapse = "#"))
    } else { 
      #when geneNum >1, the read should be a fusion. There may be some bugs since we have not simulated fusions
      #temp.fusion="Y"
      temp.strand=lapply(outGene,function(i) sapply(i,function(x) paste0(unique(temp.info$CDS_strand[temp.info$gene==x]),collapse = "")))
      temp.geneData=lapply(outGene, function(i) lapply(i,function(x) temp.info[temp.info$gene==x,]))
      temp.labels=lapply(1:length(temp.geneData),function(i) {
        temp.source=lapply(temp.geneData[[i]], function(x) perGene(x,allCDS=allCDS,tol=tol))
        return(c(paste0(sapply(temp.source,"[[",1),collapse = "&"), paste0(sapply(temp.source,"[[",2),collapse = "&"),
                 paste0(sapply(temp.source,"[[",3),collapse = "&"),paste0(sapply(temp.source,"[[",4),collapse = "&"),
                 paste0(sapply(temp.source,"[[",5),collapse = "&"),paste0(sapply(temp.source,"[[",6),collapse = "&"),
                 paste0(sapply(temp.source,"[[",7),collapse = "&"),paste0(sapply(temp.source,"[[",8),collapse = "&"),
                 paste0(outGene[[i]],collapse = "&"),temp.fusion,paste0(temp.strand[[i]],collapse = "&")))
      })
      temp.label=c(paste0(sapply(temp.labels,"[[",1),collapse = "#"),paste0(sapply(temp.labels,"[[",2),collapse = "#"),paste0(sapply(temp.labels,"[[",3),collapse = "#"),
                   paste0(sapply(temp.labels,"[[",4),collapse = "#"),paste0(sapply(temp.labels,"[[",5),collapse = "#"),paste0(sapply(temp.labels,"[[",6),collapse = "#"),
                   paste0(sapply(temp.labels,"[[",7),collapse = "#"),paste0(sapply(temp.labels,"[[",8),collapse = "#"),paste0(sapply(temp.labels,"[[",9),collapse = "#"),temp.fusion,
                   paste0(sapply(temp.labels,"[[",11),collapse = "#"))
    }
  }
  print(unique(temp.info$SampleID))
  return(temp.label)
}


##################################################################################
#################Generate reports for 90%covered samples########################
########################################################################
coverRep.Gen=function(match_info,CDSref,outfile,tol=9,mc=30) {
  sampleIDs=unique(match_info$SampleID) #3538
  sourceList=lapply(sampleIDs, function(x) sample_source(match_info[match_info$SampleID==x,],CDSref)) #for time comparison
  #names(sourceList)=sampleIDs
  #list.save(sourceList,"/zfs2/sliu/wenjia/SimulationRO1/S4/Output_EXON/Output_Len/NormFusionIsofLen200/MINIMAP2/sourceList.RData")
  print("All reads done!")
  # sampleReport=data.frame(SampleID=sampleIDs,gene=sapply(sourceList,"[[",9),gene_strand=sapply(sourceList,"[[",11),isoform=sapply(sourceList,"[[",1),position=sapply(sourceList,"[[",2),
  #                         nblock=sapply(sourceList,"[[",3),NO.CDS=sapply(sourceList,"[[",7),nCDS_isof=sapply(sourceList,"[[",4),length_isof=sapply(sourceList,"[[",5),
  #                         fusion=sapply(sourceList,"[[",10),note=sapply(sourceList,"[[",6),type=sapply(sourceList,"[[",8))
  sampleReport=data.frame(SampleID=sampleIDs,gene=sapply(sourceList,"[[",9),gene_strand=sapply(sourceList,"[[",11),isoform=sapply(sourceList,"[[",1),position=sapply(sourceList,"[[",2),
                          nblock=sapply(sourceList,"[[",3),NO.Exon=sapply(sourceList,"[[",7),nExon_isof=sapply(sourceList,"[[",4),length_isof=sapply(sourceList,"[[",5),
                          fusion=sapply(sourceList,"[[",10),note=sapply(sourceList,"[[",6),type=sapply(sourceList,"[[",8))
  
  write.csv(sampleReport,outfile,row.names = F)
  return(sampleReport)
}
#addRep=coverRep.Gen(match_info_add,allCDS,coverRep[1])


################################################
######Functions for AA annotation#############
##############################################
AAannote.isof=function(temp.isof,isoformAA) {
  temp.isofs=unlist(strsplit(temp.isof,"[||]"))
  temp.AAseqs=c()
  temp.AAnotes=c()
  for (i in temp.isofs) {
    if (i=="undefined") {
      temp.AAseq="undefined"
      temp.AAnote="undefined"
    }
    else {
      temp.AAind=match(i,isoformAA$isoformID)
      if (is.na(temp.AAind)) {
        temp.AAseq="noncoding"
        temp.AAnote="noncoding"
      }
      else {
        temp.AAseq=isoformAA$AAseq[temp.AAind]
        temp.AAnote=isoformAA$note[temp.AAind]
      }
    }
    temp.AAseqs=c(temp.AAseqs,temp.AAseq)
    temp.AAnotes=c(temp.AAnotes,temp.AAnote)
  }
  return(c(paste0(temp.AAseqs,collapse = "||"),paste0(temp.AAnotes,collapse = "||")))
}
AAannote=function(readRep,isoformAA) {
  temp.isof=as.character(readRep[4])
  temp.fusion=readRep[10]
  if (is.na(temp.fusion)) {
    temp.AA=c("undefined","undefined")
  }
  else {
    if (temp.fusion=="Y"|temp.fusion=="fusion with overlapping") {
      temp.combine=unlist(strsplit(temp.isof,"#"))
      comb.AAseqs=c()
      comb.AAnotes=c()
      for (j in temp.combine) {
        temp.parts=unlist(strsplit(j,"&"))
        part.AA=lapply(temp.parts,function(x) AAannote.isof(x,isoformAA))
        part.AAseqs=paste0(sapply(part.AA,"[[",1),collapse = "&")
        part.AAnotes=paste0(sapply(part.AA,"[[",2),collapse = "&")
        comb.AAseqs=c(comb.AAseqs,part.AAseqs)
        comb.AAnotes=c(comb.AAnotes,part.AAnotes)
      }
      
      temp.AA=c(paste0(comb.AAseqs,collapse = "#"),paste0(comb.AAnotes,collapse = "#"))
    }
    else {
      temp.combine=unlist(strsplit(temp.isof,"#"))
      comb.AA=lapply(temp.combine,function(x) AAannote.isof(x,isoformAA))
      temp.AA=c(paste0(sapply(comb.AA,"[[",1),collapse = "#"),paste0(sapply(comb.AA,"[[",2),collapse = "#"))
    }
  }
  
  return(temp.AA)
}

#####Fusion filtering ######
partlength=function(parts) {
  partAs=sapply(parts,"[[",1)
  partBs=sapply(parts,"[[",2)
  blocksAs=strsplit(partAs,";")
  blocksBs=strsplit(partBs,";")
  partAslen=sapply(blocksAs,function(p) {
    blockAs=strsplit(p,":")
    blockAstart=as.numeric(sapply(blockAs,"[[",2))
    blockAend=as.numeric(sapply(blockAs,"[[",3))
    blocklen=abs(blockAend-blockAstart)
    partlen=sum(blocklen)
    return(partlen)
  })
  
  partBslen=sapply(blocksBs,function(p) {
    blockBs=strsplit(p,":")
    blockBstart=as.numeric(sapply(blockBs,"[[",2))
    blockBend=as.numeric(sapply(blockBs,"[[",3))
    blocklen=abs(blockBend-blockBstart)
    partlen=sum(blocklen)
    return(partlen)
  })
  partlen=cbind(partAslen,partBslen)
  colnames(partlen)=c("partA_len","partB_len")
  return(partlen)
}

filteredRep=function(reportPath,fusionfiltPath,min.len=10,pseudogenes,rootNames,species="hg38",Hm_Mm_match) {
  report=read.csv(reportPath)
  allannot=strsplit(report$position, "#")
  nannot=sapply(allannot,length)
  #print(table(nannot))
  allparts=strsplit(sapply(allannot,"[[",1), "&")
  nparts=sapply(allparts,length)
  #print(which(nparts>2))
  #print(table(nparts))
  
  #criteria I
  if (sum(nparts > 1) > 0) {
    
    if (sum(nparts==2)>0) {
      fusionPart=partlength(allparts[nparts==2])
      Rep_2p=report[nparts==2,]
      report$fusionlen="pass"
      report$fusionlen[(report$SampleID %in% Rep_2p$SampleID[apply(fusionPart,1,min)<=min.len])]="failed"
    }
    if (sum(nparts > 2)>0) {
      partslen.list=lapply(allparts[nparts > 2],function(x) {
        np=length(x)
        partslen=c()
        for (i in 1:np) {
          blocks=strsplit(x[i],";")
          partlen=sapply(blocks,function(p) {
            blocklist=strsplit(p,":")
            blockstart=as.numeric(sapply(blocklist,"[[",2))
            blockend=as.numeric(sapply(blocklist,"[[",3))
            blocklen=abs(blockend-blockstart)
            partlen=sum(blocklen)
            return(partlen)
          })
          partslen[i]=partlen
        }
        return(partslen)
      })
      Rep_ge3p=report[nparts > 2,]
      if ("fusionlen" %in% colnames(report)) {
        report$fusionlen[(report$SampleID %in% Rep_ge3p$SampleID[sapply(partslen.list,min)<=min.len])]="failed"
        report$fusionlen[!(report$SampleID %in% Rep_2p$SampleID|report$SampleID %in% Rep_ge3p$SampleID)]=NA
      }
      else {
        report$fusionlen="pass"
        report$fusionlen[(report$SampleID %in% Rep_ge3p$SampleID[sapply(partslen.list,min)<=min.len])]="failed"
        report$fusionlen[!(report$SampleID %in% Rep_ge3p$SampleID)]=NA
      }
    }
  }
  else {
    report$fusionlen=NA
  }
  
  #criteria II
  if (!is.null(pseudogenes)) {
    pseudoID=sapply(1:length(report$gene),function(i) {
      x=report$gene[i]
      annolist=strsplit(x,"#")
      annopseudo=sapply(annolist,function(g) {
        genes=unlist(strsplit(g,"&"))
        #npseudo=sum(genes %in% unique(pseudoName))
        npseudo=sum(genes %in% unique(pseudogenes))
      })
      return(all(annopseudo!=0))
    })
    report$pseudogene="pass"
    report$pseudogene[pseudoID]="failed" 
  }
  else {
    report$pseudogene="pass"
  }
  
  #criteria III
  if (!is.null(rootNames)) {
    if (sum(report$fusion!="N")>0) {
      if (species=="mm10") {
        familyID=sapply(1:sum(report$fusion!="N"),function(i) {
          x=report$gene[report$fusion!="N"][i]
          annolist=strsplit(x,"#")
          annofam=sapply(annolist,function(g) {
            genes=unlist(strsplit(g,"&"))
            Hm_match=keep(Hm_Mm_match,function(x) x$mouse %in% genes)
            if (length(Hm_match)==0) return(0)
            FAMlist=sapply(rootNames$Common.root.gene.symbol,function(g) {
              return(sum(sapply(Hm_match,function(i) length(str_which(unlist(i),paste0("^",g,"[^a-zA-Z]+"))))))
            })
            FAMid=sapply(FAMlist,function(j) {
              return(j==length(genes))
            })
            nfam=sum(FAMid)
            return(nfam)
          })
          return(all(annofam!=0))
        })
      }
      else {
        familyID=sapply(1:sum(report$fusion!="N"),function(i) {
          x=report$gene[report$fusion!="N"][i]
          annolist=strsplit(x,"#")
          annofam=sapply(annolist,function(g) {
            genes=unlist(strsplit(g,"&"))
            FAMlist=sapply(rootNames$Common.root.gene.symbol,function(g) grep(paste0("^",g,"[^a-zA-Z]+"),genes))
            FAMid=sapply(FAMlist,function(j) {
              return(length(j)==length(genes))
            })
            nfam=sum(FAMid)
          })
          return(all(annofam!=0))
        })
      }
      report$FamGene="pass"
      report$FamGene[report$fusion=="N"]=NA
      report$FamGene[report$fusion!="N"][familyID]="failed" #1481
    }
    
    else {
      report$FamGene=NA
    }
  }
  else {
    report$FamGene="pass"
    report$FamGene[report$fusion=="N"]=NA
  }
  
  write.csv(report,reportPath,row.names = F)
  write.csv(report[report$fusionlen!="failed" & report$pseudogene!="failed" & report$FamGene!="failed" & report$fusion!="N" ,],fusionfiltPath,row.names = F)
  return(report)
}

###Fusion filtering >=3 genes ######
#####Fusion filtering ######

# filteredRep=function(reportPath,fusionfiltPath,min.len=10,pseudogenes,rootNames) {
#   report=read.csv(reportPath)
#   positions=report$position
#   allannot=strsplit(positions, "#")
#   nannot=sapply(allannot,length)
#   #print(table(nannot))
#   allparts=strsplit(sapply(allannot,"[[",1), "&")
#   nparts=sapply(allparts,length)
#   parts=allparts[nparts > 2]
#   partslen.list=lapply(parts,function(x) {
#     np=length(x)
#     partslen=c()
#     for (i in 1:np) {
#       blocks=strsplit(x[i],";")
#       partlen=sapply(blocks,function(p) {
#         blocklist=strsplit(p,":")
#         blockstart=as.numeric(sapply(blocklist,"[[",2))
#         blockend=as.numeric(sapply(blocklist,"[[",3))
#         blocklen=abs(blockend-blockstart)
#         partlen=sum(blocklen)
#         return(partlen)
#       })
#       partslen[i]=partlen
#     }
#     return(partslen)
#   })
#   #criteria I
#   Rep_ge3p=report[nparts > 2,]
#   Rep_ge3p$fusionlen="pass"
#   Rep_ge3p$fusionlen[sapply(partslen.list,min)<=min.len]="failed"
#   
#   #criteria II
#   pseudoID=sapply(1:length(Rep_ge3p$gene),function(i) {
#     x=Rep_ge3p$gene[i]
#     annolist=strsplit(x,"#")
#     annopseudo=sapply(annolist,function(g) {
#       genes=unlist(strsplit(g,"&"))
#       #npseudo=sum(genes %in% unique(pseudoName))
#       npseudo=sum(genes %in% unique(pseudogenes))
#     })
#     return(all(annopseudo!=0))
#   })
#   Rep_ge3p$pseudogene="pass"
#   Rep_ge3p$pseudogene[pseudoID]="failed" 
#   
#   #criteria III
#   familyID=sapply(Rep_ge3p$gene,function(x) {
#     annolist=strsplit(x,"#")
#     annofam=sapply(annolist,function(g) {
#       genes=unlist(strsplit(g,"&"))
#       FAMlist=sapply(rootNames$Common.root.gene.symbol,function(g) grep(paste0("^",g,"[^a-zA-Z]+"),genes))
#       FAMid=sapply(FAMlist,function(j) {
#         return(length(j)==length(genes))
#       })
#       nfam=sum(FAMid)
#     })
#     return(all(annofam!=0))
#   })
#   Rep_ge3p$FamGene="pass"
#   Rep_ge3p$FamGene[familyID]="failed" 
#   
#   write.csv(report,reportPath,row.names = F)
#   write.csv(report[report$fusionlen!="failed" & report$pseudogene!="failed" & report$FamGene!="failed" & report$fusion!="N" ,],fusionfiltPath,row.names = F)
#   return(report)
# }



###########################################
######isoform and fusion quantification#############
###########################################

#note that only those reads with multiple isoform annotation has missing z in Q function
# the reads with discontinuous/novel isoform annotation are considered as other isoform (which means known z)
EMperGene=function(read_annot,tol=1e-5,max.iter=200) {
  isoforms=strsplit(unlist(strsplit(read_annot$isoform,"#")),"\\|\\|")
  theta0=prop.table(table(unlist(isoforms)))
  if (!("undefined" %in% names(theta0))) {
    theta0["undefined"]=0
  }
  nr=nrow(read_annot)
  mat_update=matrix(,nrow=nr,ncol = length(theta0))
  colnames(mat_update)=names(theta0)
  theta=Inf
  niter=0
  
  while(sqrt(sum((as.vector(theta0)-as.vector(theta))^2))>tol & (niter < max.iter)) {
    if (niter>0) {theta0=theta}
    niter=niter+1
    for (i in 1:nr) {
      r=read_annot[i,]
      
      # len=sum(sapply(strsplit(strsplit(r[1,5],"#")[[1]],";"),function(x) {
      #   pos=unlist(strsplit(x,":"))
      #   lens=abs(as.numeric(pos[2])-as.numeric(pos[3]))
      #   return(lens)
      # }))
      # print(len)
      #f=summary(km_fit,len-1)$surv-summary(km_fit,len)$surv
      zknown=(!grepl("#",r[1,4])) & (!grepl("\\|\\|",r[1,4]))
      # & (r[1,11]=="continuous CDS and edge-matching")
      notes=unlist(strsplit(r[1,11],"#"))
      if (any(notes=="continuous CDS and edge-matching")) {
        if(zknown) {mat_update[i,r[1,4]]=1}
        else {
          isofs=unlist(strsplit(unlist(strsplit(r[1,4],"#")),"\\|\\|"))
          mat_update[i,isofs]=theta0[isofs]/sum(theta0[isofs])
        }
      }
      else {
        mat_update[i,"undefined"]=1
      }
      
    }
    
    theta=apply(mat_update,2,function(t) sum(t,na.rm = T)/length(t))
    #print(theta)
    
  }
  #print(niter)
  return(list(prop=theta,counts=theta*nr,niter=niter))
}

quantBygene=function(report.multigenes,report.uniquegene,tol=1e-5,max.iter=200) {
  readBYgeneList=list()
  n_multigenes=length(unique(report.multigenes$gene))
  quantBYgenesList=lapply(1:n_multigenes, function(i) {
      g=unique(report.multigenes$gene)[i]
      genes=unlist(strsplit(g,"#"))
      read.index=which(report.uniquegene$gene %in% genes)
      if (length(read.index)!=0) {
        reads=rbind(report.multigenes[report.multigenes$gene==g,],report.uniquegene[read.index,])
      }
      else {
        reads=rbind(report.multigenes[report.multigenes$gene==g,])
      }

      quant=EMperGene(read_annot = reads,tol=tol,max.iter=max.iter)
      return(quant)
    })
  
  names(quantBYgenesList)=unique(report.multigenes$gene)
  multigenes=unique(unlist(strsplit(names(quantBYgenesList),"#")))
  report.uniquegene=report.uniquegene[-which(report.uniquegene$gene %in% multigenes),]
  
  quantBYgeneList=lapply(1:length(unique(report.uniquegene$gene)),function(j) {
    g=unique(report.uniquegene$gene)[j]
    reads=report.uniquegene[report.uniquegene$gene==g,]
    quant=EMperGene(read_annot = reads,tol=tol,max.iter=max.iter)
    return(quant)
  })
  
  names(quantBYgeneList)=unique(report.uniquegene$gene)
  
  return(c(quantBYgenesList,quantBYgeneList))
}

EMperGenefusion=function(read_annot,tol=1e-5,max.iter=200) {
  genes=unlist(strsplit(read_annot$isoform,"#"))
  isoforms=unlist(as.vector(sapply(genes, function(g) {
    pair=unlist(strsplit(g,"&"))
    # head=unlist(strsplit(pair[1],"\\|\\|"))
    # tail=unlist(strsplit(pair[2],"\\|\\|")) #only for fusions of two genes/isoforms
    # return(as.vector(sapply(head,function(x) paste0(x,"&",tail))))
    pairlist=lapply(pair, function(x) unlist(strsplit(x,"\\|\\|")))
    head=pairlist[[1]]
    for (i in 2:length(pairlist)) {
      head=as.vector(sapply(head,function(x) paste0(x,"&",pairlist[[i]])))
    }
    return(head)
    
  })))
  isoforms[grepl("undefined",isoforms)]="undefined"
  isoforms[grepl("NA",isoforms)]="undefined"
  theta0=prop.table(table(unlist(isoforms)))
  if (!("undefined" %in% names(theta0))) {
    theta0["undefined"]=0
  }
  nr=nrow(read_annot)
  mat_update=matrix(,nrow=nr,ncol = length(theta0))
  colnames(mat_update)=names(theta0)
  theta=Inf
  niter=0
  
  while(sqrt(sum((as.vector(theta0)-as.vector(theta))^2))>tol & (niter < max.iter)) {
    if (niter>0) {theta0=theta}
    niter=niter+1
    for (i in 1:nr) {
      r=read_annot[i,]
      
      # len=sum(sapply(strsplit(strsplit(r[1,5],"#")[[1]],";"),function(x) {
      #   pos=unlist(strsplit(x,":"))
      #   lens=abs(as.numeric(pos[2])-as.numeric(pos[3]))
      #   return(lens)
      # }))
      # print(len)
      #f=summary(km_fit,len-1)$surv-summary(km_fit,len)$surv
      zknown=(!grepl("#",r[1,4])) & (!grepl("\\|\\|",r[1,4]))
      # & (r[1,11]=="continuous CDS and edge-matching")
      notes=unlist(strsplit(r[1,11],"#"))
      ng=str_count(notes[1],"&")
      if (any(notes==paste0(rep("continuous CDS and edge-matching",ng+1),collapse = "&") & !is.na(notes))) {
        if(zknown) {mat_update[i,r[1,4]]=1}
        else {
          isofmulti=unlist(strsplit(r[1,4],"#"))
          isofs=unlist(as.vector(sapply(isofmulti, function(g) {
            pair=unlist(strsplit(g,"&"))
            # head=unlist(strsplit(pair[1],"\\|\\|"))
            # tail=unlist(strsplit(pair[2],"\\|\\|"))
            # return(as.vector(sapply(head,function(x) paste0(x,"&",tail))))
            pairlist=lapply(pair, function(x) unlist(strsplit(x,"\\|\\|")))
            head=pairlist[[1]]
            for (i in 2:length(pairlist)) {
              head=as.vector(sapply(head,function(x) paste0(x,"&",pairlist[[i]])))
            }
            return(head)
          })))
          isofs[grepl("undefined",isofs)]="undefined"
          isofs[grepl("NA",isofs)]="undefined"
          isofs=unique(isofs)
          mat_update[i,isofs]=theta0[isofs]/sum(theta0[isofs])
        }
      }
      else {
        mat_update[i,"undefined"]=1
      }
      
    }
    
    theta=apply(mat_update,2,function(t) sum(t,na.rm = T)/length(t))
    #print(theta)
    
  }
  #print(niter)
  return(list(prop=theta,counts=theta*nr,niter=niter))
}

quantBygenefusion=function(report.multigenes,report.uniquegene,tol=1e-5,max.iter=200) {
  readBYgeneList=list()
  n_multigenes=length(unique(report.multigenes$gene))
  
  quantBYgenesList=lapply(1:n_multigenes, function(i) {
    g=unique(report.multigenes$gene)[i]
    genes=unlist(strsplit(g,"#"))
    read.index=which(report.uniquegene$gene %in% genes)
    if (length(read.index)!=0) {
      reads=rbind(report.multigenes[report.multigenes$gene==g,],report.uniquegene[read.index,])
    }
    else {
      reads=rbind(report.multigenes[report.multigenes$gene==g,])
    }
    
    quant=EMperGenefusion(read_annot = reads,tol=tol,max.iter=max.iter)
    return(quant)
  })
  
  names(quantBYgenesList)=unique(report.multigenes$gene)
  multigenes=unique(unlist(strsplit(names(quantBYgenesList),"#")))
  if (length(which(report.uniquegene$gene %in% multigenes))!=0) {
    report.uniquegene=report.uniquegene[-which(report.uniquegene$gene %in% multigenes),]
  }
  
  
  quantBYgeneList=lapply(1:length(unique(report.uniquegene$gene)),function(j) {
    g=unique(report.uniquegene$gene)[j]
    reads=report.uniquegene[report.uniquegene$gene==g,]
    quant=EMperGenefusion(read_annot = reads,tol=tol,max.iter=max.iter)
    return(quant)
  })
  
  names(quantBYgeneList)=unique(report.uniquegene$gene)
  
  return(c(quantBYgenesList,quantBYgeneList))
}

quantList.gen=function(Report,fusionRep,Isof.quantfile,Isof.quantRData,fusion.quantfile,fusion.quantRData,gtf.dat,tol=1e-5,max.iter=200,mc=30){
  isof.ind=which(Report$fusion=="N")
  # fusion.ind=which(Report$fusion=="Y")
  if (length(isof.ind)!=0) {
    ######### Isoform quantification #############
    Report.isof=Report[Report$fusion=="N",]
    multigenes.ind=grep("#",Report.isof$gene)
    Report.multigenes=Report.isof[multigenes.ind,]
    Report.uniquegene=Report.isof[-multigenes.ind,]
    quantList=quantBygene(report.multigenes=Report.multigenes,report.uniquegene=Report.uniquegene,tol=tol,max.iter=max.iter)
    list.save(quantList,Isof.quantRData)
    # quantList=list.load(quantRData)
    
    count=c()
    prop=c()
    group=c()
    isoform=c()
    for (i in 1:length(quantList)) {
      count=c(count,quantList[[i]]$counts)
      isoform=c(isoform,names(quantList[[i]]$prop))
      prop=c(prop,quantList[[i]]$prop)
      group=c(group,rep(names(quantList)[i],length(quantList[[i]]$counts)))
    }
    
    gene=lapply(isoform,function(x) unique(gtf.dat$V13[gtf.dat$V16==x | gtf.dat$V19==x]))
    gene[sapply(gene,length)==0]="undefined"
    isof_quant.dat=data.frame(isoform=isoform, gene=unlist(gene),group=group,prop=prop,count=count)
    write.csv(isof_quant.dat,file=Isof.quantfile)
  }
  else {
    isof_quant.dat=NULL
    message("No isoform detected in the data")
  }
  
  if (is.data.frame(fusionRep)) {
    ######## Fusion quantification ################
    multigenes.ind=grep("#",fusionRep$gene)
    fusionRep.multigenes=fusionRep[multigenes.ind,]
    fusionRep.uniquegene=fusionRep[-multigenes.ind,]
    
    quantList=quantBygenefusion(report.multigenes=fusionRep.multigenes,report.uniquegene=fusionRep.uniquegene,tol=tol,max.iter=max.iter)
    
    list.save(quantList,fusion.quantRData)
    # quantList=list.load(quantRData)
    
    count=c()
    prop=c()
    group=c()
    isoform=c()
    for (i in 1:length(quantList)) {
      count=c(count,quantList[[i]]$counts)
      isoform=c(isoform,names(quantList[[i]]$prop))
      prop=c(prop,quantList[[i]]$prop)
      group=c(group,rep(names(quantList)[i],length(quantList[[i]]$counts)))
    }
    
    gene=lapply(isoform,function(x) {
      isofs=unlist(strsplit(x,"&"))
      genes=c()
      for (i in 1:length(isofs)) {
        if (isofs[i]=="undefined") {
          genes=c(genes,"undefined")
        }
        else {
          genes=c(genes,unique(gtf.dat$V13[gtf.dat$V16==isofs[i] | gtf.dat$V19==isofs[i]]))
        }
      }
      return(paste0(genes,collapse ="&"))
    })
    gene=unlist(gene)
    fusion_quant.dat=data.frame(isoform=isoform,gene=gene, group=group,prop=prop,count=count)
    
    if (length(isof.ind)!=0) {
      isof_counts=lapply(fusion_quant.dat$isoform,function(x) {
        isofs=unlist(strsplit(x,"&"))
        counts=c()
        for (i in 1:length(isofs)) {
          if (isofs[i]=="undefined") {
            counts=c(counts,"undefined")
          }
          else {
            if (is.na(match(isofs[i],isof_quant.dat$isoform))) {
              counts=c(counts,0)
            } else {
              counts=c(counts,round(isof_quant.dat$count[match(isofs[i],isof_quant.dat$isoform)],2))
            }
          }
        }
        return(paste0(counts,collapse ="&"))
      })
      isof_counts=unlist(isof_counts)
      
      gene_counts=lapply(fusion_quant.dat$gene,function(x) {
        genes=unlist(strsplit(x,"&"))
        counts=c()
        for (i in 1:length(genes)) {
          if (genes[i]=="undefined") {
            counts=c(counts,"undefined")
          }
          else {
            if (is.na(match(genes[i],isof_quant.dat$gene))) {
              counts=c(counts,0)
            }else {
              counts=c(counts,round(sum(isof_quant.dat$count[match(genes[i],isof_quant.dat$gene)]),2))
            }
          }
        }
        return(paste0(counts,collapse ="&"))
      })
      gene_counts=unlist(gene_counts)
    } 
    else {
      isof_counts=rep("0&0",dim(fusion_quant.dat)[1])
      gene_counts=rep("0&0",dim(fusion_quant.dat)[1])
    }
    
    fusion_quant.dat$isoform_counts=isof_counts
    fusion_quant.dat$gene_counts=gene_counts
    write.csv(fusion_quant.dat,fusion.quantfile)
    
  }
  else {
    message("No fusion detected in the data")
    fusion_quant.dat=NULL
  }
  
  return(list(fusion_quant=fusion_quant.dat,isof_quant=isof_quant.dat))
}
