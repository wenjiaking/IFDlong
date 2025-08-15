
suppressPackageStartupMessages({
  install_and_load <- function(pkgs) {
    # Set CRAN mirror if not already set
    if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
      options(repos = c(CRAN = "https://cloud.r-project.org"))
    }
    
    for (p in pkgs) {
      if (!requireNamespace(p, quietly = TRUE)) {
        install.packages(p)
      }
      library(p, character.only = TRUE)
    }
  }
  
  install_and_load(c("Rcpp", "data.table", "parallel", "stringr", "rlist"))
})




#### parameter
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 12)

PATH <- args[1]
sampleName <- args[2]
Aligner <- args[3]
buffer <- as.integer(args[4])
anchorLen <- as.integer(args[5])
refEXON <- args[6] #"./refData/allexon_NO.bed"
refAA <- args[7] #"./refData/isoformAA_Ref.txt"
refPseudo <- args[8] #"./refData/pseudogenes.rds"
refRoot <- args[9] #"./refData/rootNames.txt"
refHMmatch <- args[10] #"./refData/Hm_Mm_match.rds"
species <- args[11] #hg38 or mm10
ncores <- as.integer(args[12])
speedup <- args[13]

# Compile C++ (point to where you saved the file)
Rcpp::sourceCpp(speedup)

### ref files
allCDS <- fread(refEXON, col.names = c("chr", "start", "end", "name", "score", "strand"))
allCDS[, isoform := tstrsplit(name, "__", keep=1)]

isoformAA_Ref <- fread(refAA, header=TRUE)
isoformAA_Ref[, isoformID := tstrsplit(isoformID, "__", keep=1)]

read_if_exists <- function(path, fun) if (file.exists(path)) fun(path) else NULL
pseudogenes <- read_if_exists(refPseudo, readRDS)
rootNames   <- read_if_exists(refRoot, fread)
Hm_Mm_match <- read_if_exists(refHMmatch, readRDS)

#### input files
interSbedfile=paste0(PATH,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS.bed")
intergenebedfile=paste0(PATH,"/",Aligner,"/",sampleName,"_mapped_woSecond_geneTol500intersectS.bed")
coverSReOut=paste0(PATH,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_Rep.csv")


uncoverFilter <- function(interbedfile, intergenebedfile) {
  intersect_tab <- read.table(interbedfile, stringsAsFactors = FALSE)
  colnames(intersect_tab) <- c("chr", "start", "end", "SampleID", "score", "strand", 
                               "CDS_chr", "CDS_start", "CDS_end", "CDS_name", 
                               "CDS_score", "CDS_strand", "n_base")
  
  match_intersect <- subset(intersect_tab, CDS_name != ".")
  CDS_info <- strsplit(match_intersect$CDS_name, "__")
  
  match_info <- match_intersect[, c(1:4, 6:9, 12, 13)]
  match_info$gene    <- sapply(CDS_info, `[[`, 6)
  match_info$isoform <- sapply(CDS_info, `[[`, 1)
  match_info$order   <- sapply(CDS_info, `[[`, 7)
  
  if (inherits(try(read.table(intergenebedfile), silent = TRUE), "try-error")) {
    read_info <- match_info
  } else {
    uncover_info <- read.table(intergenebedfile, stringsAsFactors = FALSE)
    colnames(uncover_info) <- c("chr", "start", "end", "SampleID", "score", "strand", 
                                "gene_chr", "gene_start", "gene_end", "gene_name", 
                                "gene_score", "gene_strand", "n_base")
    
    unmatch_info <- uncover_info[, c(1:4, 6:9, 12, 13)]
    unmatch_info$gene <- ifelse(uncover_info$gene_name != ".", 
                                sapply(strsplit(uncover_info$gene_name[uncover_info$gene_name != "."], "__"), `[[`, 1), 
                                "undefined")
    
    unmatch_info$isoform <- "undefined"
    unmatch_info$order   <- "undefined"
    colnames(unmatch_info) <- c("chr", "start", "end", "SampleID", "strand", 
                                "CDS_chr", "CDS_start", "CDS_end", "CDS_strand", 
                                "n_base", "gene", "isoform", "order")
    
    read_info <- rbind(match_info, unmatch_info)
  }
  
  read_info$CDS_strand[read_info$CDS_strand == "."] <- "undefined"
  read_info$CDS_chr[read_info$CDS_chr == "."]       <- "undefined"
  
  message(" Extract the CDS Covered Alignments Done!")
  return(as.data.frame(read_info))
}

testContinue <- function(vec) {
  numbers <- suppressWarnings(as.numeric(vec))  # suppress coercion warnings
  if (any(is.na(numbers))) {
    # If any NA after coercion, consider it not continuous or handle appropriately
    return(FALSE)
  }
  result <- diff(numbers)
  return(!any(abs(result) != 1))
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


# ---- Wire in the fast versions ----
testContinue <- testContinueCpp
isoformFilter <- isoformFilterCpp


# perGene: call the C++ version
perGene <- function(geneData, allCDS, tol = 9) {
  # Ensure required columns exist
  stopifnot(all(c("chr","start","end","strand","isoform","order",
                  "CDS_start","CDS_end") %in% names(geneData)))
  stopifnot(all(c("isoform","start","end") %in% names(allCDS)))
  perGeneCpp(geneData, allCDS, tol)
}

# filtGene: call the C++ version (same return shape as your R list)
filtGene <- function(temp.info) {
  res <- filtGeneCpp(temp.info)
  # R side returns a list with the same names; no change required
  res
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

coverRep.Gen=function(match_info,CDSref,outfile,tol=9,mc=30) {
  sampleIDs=unique(match_info$SampleID) #3538
  sourceList=mclapply(sampleIDs, function(x) sample_source(match_info[match_info$SampleID==x,],CDSref,tol = tol),mc.cores=mc)
  print("All reads done!")
  sampleReport=data.frame(SampleID=sampleIDs,gene=sapply(sourceList,"[[",9),gene_strand=sapply(sourceList,"[[",11),isoform=sapply(sourceList,"[[",1),position=sapply(sourceList,"[[",2),
                          nblock=sapply(sourceList,"[[",3),NO.Exon=sapply(sourceList,"[[",7),nExon_isof=sapply(sourceList,"[[",4),length_isof=sapply(sourceList,"[[",5),
                          fusion=sapply(sourceList,"[[",10),note=sapply(sourceList,"[[",6),type=sapply(sourceList,"[[",8))
  
  write.csv(sampleReport,outfile,row.names = F)
  return(sampleReport)
}

####### AA
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


###### filter 
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
            nfam=sum(unlist(FAMid))
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
            nfam=sum(unlist(FAMid))
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

########isoform annotation##############
match_info=uncoverFilter(interSbedfile,intergenebedfile)
#CDScover_source=coverRep.Gen(match_info, CDSref=allCDS,outfile=coverSReOut,tol=buffer,mc=ncores) #mc.cores=30 in parallel.


#print(ls())

sampleIDs <- unique(match_info$SampleID)
chunk_size <- ceiling(length(sampleIDs) / 60)
chunks <- split(sampleIDs, ceiling(seq_along(sampleIDs) / chunk_size))

# Create list of data frames for each chunk
match_split <- split(match_info, match_info$SampleID)

match_info_chunks <- lapply(chunks, function(ids) {
  do.call(rbind, match_split[ids])
})


options(error = recover)
results <- mclapply(seq_along(match_info_chunks), function(i) {
  chunk <- match_info_chunks[[i]]
  outfilename <- file.path(normalizePath(PATH), Aligner, paste0("coverSReOut_chunk_", i, ".txt"))
  tryCatch({
    coverRep.Gen(chunk, CDSref=allCDS, outfile=outfilename, tol=buffer, mc=ncores)
    # Optionally, comment out or remove the file.exists check if it's causing errors
    # stopifnot(file.exists("expected_input_file.txt"))
    TRUE  # Indicate success
  }, error = function(e) {
    message(sprintf("Error in chunk %d: %s", i, e$message))
    #save(chunk, file = paste0("chunk_", i, "_error.RData")) # save for debugging
    FALSE # Indicate failure but continue
  })
}, mc.cores = ncores)


chunk_files <- paste0(PATH,"/",Aligner,"/", "coverSReOut_chunk_", seq_along(match_info_chunks), ".txt")
merged_data <- do.call(rbind, lapply(chunk_files, read.csv, stringsAsFactors = FALSE))

write.csv(merged_data, coverSReOut, row.names = FALSE)

file.remove(chunk_files)


#####AA annotation###################
isoformName=sapply(strsplit(isoformAA_Ref$isoformID,"__"),"[[",1)
isoformAA=isoformAA_Ref
isoformAA$isoformID=isoformName


###############AAannotation generating functions#########################
coverSReOut=paste0(PATH,"/", Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_Rep.csv")
fullRepOut=paste0(PATH,"/", Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_fullRep.csv")

Rep=read.csv(coverSReOut)
nr=nrow(Rep)

AArep.list=lapply(1:nr, function(r) AAannote(Rep[r,],isoformAA))
AArep=do.call(rbind,AArep.list)
Rep$AAseq=AArep[,1]
Rep$AAnote=AArep[,2]
write.csv(Rep,fullRepOut,row.names = F)


########  Filtering   Pseudo gene reference  ###############
fullRepOut=paste0(PATH,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_fullRep.csv")
fusionfiltReOut=paste0(PATH,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_fusionRep_anchor",anchorLen,"bp.filt.csv")

pseudogenes=ifelse(file.exists(refPseudo), readRDS(refPseudo), NULL)
rootNames=ifelse(file.exists(refRoot), read.table(refRoot,sep = "\t",header = T), NULL)
Hm_Mm_match=ifelse(file.exists(refHMmatch), readRDS(refHMmatch), NULL)

filtRep=filteredRep(reportPath =fullRepOut,fusionfiltPath = fusionfiltReOut,min.len=anchorLen,
                    pseudogenes=pseudogenes,rootNames=rootNames,species = species,Hm_Mm_match=Hm_Mm_match)












