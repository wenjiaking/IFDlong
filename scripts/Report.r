
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
  
  install_and_load(c("data.table", "parallel", "stringr", "rlist", "dplyr", "purrr"))
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

###### 
process_isoform_group <- function(oiso, oord, CDS_chr, CDS_start, CDS_end, CDS_strand) {
  if (any(oiso == "undefined")) {
    return(list(
      isoform = "undefined",
      nblock = length(oiso),
      NO.Exon = NA_character_,
      note = "no full-covered isoform",
      type = "novel with addition",
      position = paste(paste(CDS_chr, CDS_start, CDS_end, CDS_strand, sep = ":"), collapse = ";")
    ))
  }
  
  oord <- as.numeric(oord)
  iso_split <- split(oord, oiso)
  iso_counts <- vapply(iso_split, length, integer(1))
  iso_continuous <- vapply(iso_split, function(x) all(diff(sort(x)) == 1), logical(1))
  
  if (any(iso_continuous)) {
    continuous_iso <- names(iso_counts[iso_continuous])
    chosen_iso <- continuous_iso[which.max(iso_counts[continuous_iso])]
  } else {
    chosen_iso <- names(iso_counts)[which.max(iso_counts)]
  }
  
  iso_ordered <- names(sort(iso_counts, decreasing = TRUE))
  
  note <- if (any(iso_continuous)) "continuous CDS and edge-matching" else "discontinuous CDS"
  type <- if (any(iso_continuous)) "normal" else "novel with deletion"
  
  list(
    isoform = paste(iso_ordered, collapse = "||"),
    nblock = length(iso_split[[chosen_iso]]),
    NO.Exon = paste(sort(as.numeric(iso_split[[chosen_iso]])), collapse = "-"),
    note = note,
    type = type,
    position = paste(paste(CDS_chr, CDS_start, CDS_end, CDS_strand, sep = ":"), collapse = ";")
  )
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
Sys.time()
print(dim(match_info))

block_match <- function(chr_val, start_val, end_val, isoform_val, ref_df, buffer=9) {
ref_sub <- ref_df[chr == chr_val & isoform == isoform_val]
if (nrow(ref_sub) == 0) return(FALSE)
any((start_val >= (ref_sub$start - buffer) & start_val <= (ref_sub$end + buffer)) &
      (end_val >= (ref_sub$start - buffer) & end_val <= (ref_sub$end + buffer)))
}

gene_counts <- match_info%>%
 filter(gene != "undefined")%>%
 group_by(SampleID) %>%
 summarise(n_gene = n_distinct(gene), .groups = "drop")

# Join back to the original data
match_info_with_counts <- match_info%>%
 filter(gene != "undefined")%>%
 left_join(gene_counts, by = "SampleID")

# Separate into two data frames
df1 <- match_info_with_counts %>% filter(n_gene == 1)%>%select(-n_gene)  # Only one gene
df2 <- match_info_with_counts %>% filter(n_gene > 1)%>%select(-n_gene)   # fusion

print(dim(df1))
print(dim(df2))

Sys.time()
df1_summary <- df1%>%
 group_by(SampleID)%>%
 summarise(
   gene = first(gene),
   gene_strand = first(strand),
   
   isoform_list = list(isoform),
   order_list = list(order),
   
   # --- isoform selection logic with continuity priority ---
   isoform = {
     oiso <- unlist(isoform_list)
     oord <- as.numeric(unlist(order_list))
     if (any(oiso == "undefined")) {
       "undefined"
     } else {
       # Split exon orders by isoform
       iso_split <- split(oord, oiso)
       
       # Determine continuity for each isoform
       iso_continuous <- sapply(iso_split, function(x) all(diff(sort(as.numeric(x))) == 1))
       iso_counts <- sapply(iso_split, length)
       
       # Select best isoform: continuous first, then longest, then first alphabetically
       if (any(iso_continuous)) {
         continuous_iso <- names(iso_counts[iso_continuous])
         chosen_iso <- continuous_iso[which.max(iso_counts[continuous_iso])]
       } else {
         chosen_iso <- names(iso_counts)[which.max(iso_counts)]
       }
       
       # Combine all isoforms (sorted by rule) for output
       iso_ordered <- names(sort(iso_counts, decreasing = TRUE))
       paste(iso_ordered, collapse = "||")
     }
   },
   
   position = paste(paste(CDS_chr, CDS_start, CDS_end, CDS_strand, sep = ":"), collapse = ";"),
   
   # --- nblock: number of exons for the chosen isoform ---
   nblock = {
     oiso <- unlist(isoform_list)
     oord <- as.numeric(unlist(order_list))
     if (any(oiso == "undefined")) {
       n()
     } else {
       iso_split <- split(oord, oiso)
       iso_continuous <- sapply(iso_split, function(x) all(diff(sort(as.numeric(x))) == 1))
       iso_counts <- sapply(iso_split, length)
       if (any(iso_continuous)) {
         continuous_iso <- names(iso_counts[iso_continuous])
         chosen_iso <- continuous_iso[which.max(iso_counts[continuous_iso])]
       } else {
         chosen_iso <- names(iso_counts)[which.max(iso_counts)]
       }
       length(iso_split[[chosen_iso]])
     }
   },
   
   # --- NO.Exon: exon orders for the chosen isoform ---
   NO.Exon = {
     oiso <- unlist(isoform_list)
     oord <- as.numeric(unlist(order_list))
     if (any(oiso == "undefined")) {
       NA
     } else {
       iso_split <- split(oord, oiso)
       iso_continuous <- sapply(iso_split, function(x) all(diff(sort(as.numeric(x))) == 1))
       iso_counts <- sapply(iso_split, length)
       if (any(iso_continuous)) {
         continuous_iso <- names(iso_counts[iso_continuous])
         chosen_iso <- continuous_iso[which.max(iso_counts[continuous_iso])]
       } else {
         chosen_iso <- names(iso_counts)[which.max(iso_counts)]
       }
       paste(sort(as.numeric(iso_split[[chosen_iso]])), collapse = "-")
     }
   },
   
   nExon_isof = NA,
   length_isof = NA,
   fusion = "N",
   
   note = {
     oiso <- unlist(isoform_list)
     oord <- as.numeric(unlist(order_list))
     if (any(oiso == "undefined")) {
       "no full-covered isoform"
     } else {
       iso_split <- split(oord, oiso)
       iso_continuous <- sapply(iso_split, function(x) all(diff(sort(as.numeric(x))) == 1))
       if (any(iso_continuous)) {
         "continuous CDS and edge-matching"
       } else {
         "discontinuous CDS"
       }
     }
   },
   
   type = {
     oiso <- unlist(isoform_list)
     oord <- as.numeric(unlist(order_list))
     if (any(oiso == "undefined")) {
       "novel with addition"
     } else {
       iso_split <- split(oord, oiso)
       iso_continuous <- sapply(iso_split, function(x) all(diff(sort(as.numeric(x))) == 1))
       if (any(iso_continuous)) {
         "normal"
       } else {
         "novel with deletion"
       }
     }
   },
   
   .groups = "drop"
 ) %>%
 select(-isoform_list, -order_list)
Sys.time()



# v1 correct slower
# #### revise the continuous CDS match or unmatch
df1_summary_cont <- df1_summary%>%filter(note == "continuous CDS and edge-matching")
dim(df1_summary_cont)
#head(df1_summary_cont)


Sys.time()
df1_summary_cont <- df1_summary_cont %>%
rowwise() %>%
mutate(
  # Split position string into individual blocks as a list
  blocks = list(str_split(position, ";")[[1]]),
  # Check each block against reference
  matched = all(sapply(blocks[[1]], function(b) {
    parts <- str_split(b, ":")[[1]]
    chr <- parts[1]
    start <- as.integer(parts[2])
    end <- as.integer(parts[3])
    block_match(chr, start, end, isoform, allCDS)
  })),
  # Update note if any block does not match
  note = ifelse(matched, note, "continuous CDS and edge-unmatching")
) %>%
ungroup() %>%
select(-blocks, -matched)
Sys.time()

df1_summary[df1_summary$note == "continuous CDS and edge-matching", ] <- df1_summary_cont

print(table(df1_summary$note))


#### add ref 
isoform_summary <- allCDS[, .(
 nExon_isof = .N,  
 length_isof = sum(end - start)), by = isoform]

#### add nExon_isof and length_isof cols
df1_summary <- df1_summary %>%
 rowwise() %>%
 mutate(
   nExon_isof = if_else(
     isoform == "undefined", 
     NA_character_,
     paste(
       str_split(isoform, "\\|\\|")[[1]] %>%
         sapply(function(x) {
           val <- isoform_summary$nExon_isof[isoform_summary$isoform == x]
           if(length(val) == 0) NA else val
         }),
       collapse = "||"
     )
   ),
   length_isof = if_else(
     isoform == "undefined",
     NA_character_,
     paste(
       str_split(isoform, "\\|\\|")[[1]] %>%
         sapply(function(x) {
           val <- isoform_summary$length_isof[isoform_summary$isoform == x]
           if(length(val) == 0) NA else val
         }),
       collapse = "||"
     )
   )
 ) %>%
 ungroup()

#head(df1_summary)
#table(df1_summary$note)
#table(df1_summary$fusion)

#########################
####### the df2 fusion part
#head(df2)

Sys.time()
df2_summary <- df2%>%
 group_by(SampleID, gene)%>%
 summarise(
   gene = first(gene),
   gene_strand = first(strand),
   
   isoform_list = list(isoform),
   order_list = list(order),
   
   # --- isoform selection logic with continuity priority ---
   isoform = {
     oiso <- unlist(isoform_list)
     oord <- as.numeric(unlist(order_list))
     if (any(oiso == "undefined")) {
       "undefined"
     } else {
       # Split exon orders by isoform
       iso_split <- split(oord, oiso)
       
       # Determine continuity for each isoform
       iso_continuous <- sapply(iso_split, function(x) all(diff(sort(as.numeric(x))) == 1))
       iso_counts <- sapply(iso_split, length)
       
       # Select best isoform: continuous first, then longest, then first alphabetically
       if (any(iso_continuous)) {
         continuous_iso <- names(iso_counts[iso_continuous])
         chosen_iso <- continuous_iso[which.max(iso_counts[continuous_iso])]
       } else {
         chosen_iso <- names(iso_counts)[which.max(iso_counts)]
       }
       
       # Combine all isoforms (sorted by rule) for output
       iso_ordered <- names(sort(iso_counts, decreasing = TRUE))
       paste(iso_ordered, collapse = "||")
     }
   },
   
   position = paste(paste(CDS_chr, CDS_start, CDS_end, CDS_strand, sep = ":"), collapse = ";"),
   
   # --- nblock: number of exons for the chosen isoform ---
   nblock = {
     oiso <- unlist(isoform_list)
     oord <- as.numeric(unlist(order_list))
     if (any(oiso == "undefined")) {
       n()
     } else {
       iso_split <- split(oord, oiso)
       iso_continuous <- sapply(iso_split, function(x) all(diff(sort(as.numeric(x))) == 1))
       iso_counts <- sapply(iso_split, length)
       if (any(iso_continuous)) {
         continuous_iso <- names(iso_counts[iso_continuous])
         chosen_iso <- continuous_iso[which.max(iso_counts[continuous_iso])]
       } else {
         chosen_iso <- names(iso_counts)[which.max(iso_counts)]
       }
       length(iso_split[[chosen_iso]])
     }
   },
   
   # --- NO.Exon: exon orders for the chosen isoform ---
   NO.Exon = {
     oiso <- unlist(isoform_list)
     oord <- as.numeric(unlist(order_list))
     if (any(oiso == "undefined")) {
       NA
     } else {
       iso_split <- split(oord, oiso)
       iso_continuous <- sapply(iso_split, function(x) all(diff(sort(as.numeric(x))) == 1))
       iso_counts <- sapply(iso_split, length)
       if (any(iso_continuous)) {
         continuous_iso <- names(iso_counts[iso_continuous])
         chosen_iso <- continuous_iso[which.max(iso_counts[continuous_iso])]
       } else {
         chosen_iso <- names(iso_counts)[which.max(iso_counts)]
       }
       paste(sort(as.numeric(iso_split[[chosen_iso]])), collapse = "-")
     }
   },
   
   nExon_isof = NA,
   length_isof = NA,
   fusion = "N",
   
   note = {
     oiso <- unlist(isoform_list)
     oord <- as.numeric(unlist(order_list))
     if (any(oiso == "undefined")) {
       "no full-covered isoform"
     } else {
       iso_split <- split(oord, oiso)
       iso_continuous <- sapply(iso_split, function(x) all(diff(sort(as.numeric(x))) == 1))
       if (any(iso_continuous)) {
         "continuous CDS and edge-matching"
       } else {
         "discontinuous CDS"
       }
     }
   },
   
   type = {
     oiso <- unlist(isoform_list)
     oord <- as.numeric(unlist(order_list))
     if (any(oiso == "undefined")) {
       "novel with addition"
     } else {
       iso_split <- split(oord, oiso)
       iso_continuous <- sapply(iso_split, function(x) all(diff(sort(as.numeric(x))) == 1))
       if (any(iso_continuous)) {
         "normal"
       } else {
         "novel with deletion"
       }
     }
   },
   
   .groups = "drop"
 ) %>%
 select(-isoform_list, -order_list)
Sys.time()


df2_summary_cont <- df2_summary%>%filter(note == "continuous CDS and edge-matching")
dim(df2_summary_cont)
#head(df2_summary_cont)


Sys.time()
df2_summary_cont <- df2_summary_cont%>%
 rowwise() %>%
 mutate(
   # Split position string into individual blocks as a list
   blocks = list(str_split(position, ";")[[1]]),
   # Check each block against reference
   matched = all(sapply(blocks[[1]], function(b) {
     parts <- str_split(b, ":")[[1]]
     chr <- parts[1]
     start <- as.integer(parts[2])
     end <- as.integer(parts[3])
     block_match(chr, start, end, isoform, allCDS)
   })),
   # Update note if any block does not match
   note = ifelse(matched, note, "continuous CDS and edge-unmatching")
 ) %>%
 ungroup() %>%
 select(-blocks, -matched)
Sys.time()

df2_summary[df2_summary$note == "continuous CDS and edge-matching", ] <- df2_summary_cont

print(table(df2_summary$note))



#### add nExon_isof and length_isof cols
df2_summary <- df2_summary %>%
 rowwise() %>%
 mutate(
   nExon_isof = if_else(
     isoform == "undefined", 
     NA_character_,
     paste(
       str_split(isoform, "\\|\\|")[[1]] %>%
         sapply(function(x) {
           val <- isoform_summary$nExon_isof[isoform_summary$isoform == x]
           if(length(val) == 0) NA else val
         }),
       collapse = "||"
     )
   ),
   length_isof = if_else(
     isoform == "undefined",
     NA_character_,
     paste(
       str_split(isoform, "\\|\\|")[[1]] %>%
         sapply(function(x) {
           val <- isoform_summary$length_isof[isoform_summary$isoform == x]
           if(length(val) == 0) NA else val
         }),
       collapse = "||"
     )
   )
 ) %>%
 ungroup()

head(df2_summary)



########### v2
head(df2_summary)

df_counts <- df2_summary%>%
  group_by(SampleID)%>%
  mutate(n_rows = n())

dim(df_counts)

# df1: groups where row count > 2
df_counts21 <- df_counts%>%
  filter(n_rows > 2) %>%
  select(-n_rows)


combine_all_cols <- function(df) {
  # Ensure all columns have names
  if (any(names(df) == "")) {
    names(df)[names(df) == ""] <- paste0("V", seq_len(sum(names(df) == "")))
  }
  
  # Group by position
  pos_groups <- df %>%
    group_by(position) %>%
    summarise(across(everything(), ~list(.x)), .groups = "drop")
  
  # If only 1 unique position, combine all rows with &
  if (nrow(pos_groups) == 1) {
    combined <- pos_groups %>%
      mutate(across(everything(), ~paste(.x[[1]], collapse = "&")))
    return(combined)
  }
  
  # All pairwise position combinations
  idx <- combn(nrow(pos_groups), 2)
  
  # Initialize a named list to hold final column values
  final_combined <- setNames(vector("list", length = ncol(pos_groups)), names(pos_groups))
  final_combined <- map(final_combined, ~character(0))
  
  for (k in seq_len(ncol(idx))) {
    i <- idx[, k]
    row1 <- pos_groups[i[1], ]
    row2 <- pos_groups[i[2], ]
    
    for (colname in names(pos_groups)) {
      cross <- expand.grid(row1[[colname]][[1]], row2[[colname]][[1]], stringsAsFactors = FALSE)
      combined <- paste0(cross$Var1, "&", cross$Var2)
      final_combined[[colname]] <- c(final_combined[[colname]], combined)
    }
  }
  
  # Collapse each column by #
  final_combined <- map(final_combined, ~paste(.x, collapse = "#"))
  
  # Return as flat tibble
  tibble::as_tibble(final_combined)
}

# Apply per SampleID
df2_summary_grouped1 <- df_counts21 %>%
  group_by(SampleID) %>%
  filter(n() >= 3) %>%
  group_modify(~combine_all_cols(.x))


###################### 
# df2: groups where row count == 2
df_counts22 <- df_counts%>%
  filter(n_rows == 2)%>%
  select(-n_rows)

######### fusion 1x1 
Sys.time()
df2_summary_grouped <- df_counts22%>%
 group_by(SampleID) %>%
 summarise(
   # Check if positions are all identical
   same_position = n_distinct(position) == 1,
   
   # Merge each column accordingly
   gene = if(same_position) paste(gene, collapse = "#") else paste(gene, collapse = "&"),
   gene_strand = if(same_position) paste(gene_strand, collapse = "#") else paste(gene_strand, collapse = "&"),
   isoform = if(same_position) paste(isoform, collapse = "#") else paste(isoform, collapse = "&"),
   position = if(same_position) paste(position, collapse = "#") else paste(position, collapse = "&"),
   nblock = if(same_position) paste(nblock, collapse = "#") else paste(nblock, collapse = "&"),
   NO.Exon = if(same_position) paste(NO.Exon, collapse = "#") else paste(NO.Exon, collapse = "&"),
   nExon_isof = if(same_position) paste(nExon_isof, collapse = "#") else paste(nExon_isof, collapse = "&"),
   length_isof = if(same_position) paste(length_isof, collapse = "#") else paste(length_isof, collapse = "&"),
   note = if(same_position) paste(note, collapse = "#") else paste(note, collapse = "&"),
   type = if(same_position) paste(type, collapse = "#") else paste(type, collapse = "&"),
   
   # Fusion logic
   fusion = if(same_position) "N" else "Y",
   .groups = "drop"
 )

#head(df2_summary_grouped)
print(table(df2_summary_grouped$fusion))

collist <- c("SampleID","gene","gene_strand","isoform","position","nblock","NO.Exon","nExon_isof","length_isof","fusion","note","type")
df1_summary <- df1_summary[, collist]
df2_summary_grouped <- df2_summary_grouped[, collist]
df2_summary_grouped1 <- df2_summary_grouped1[, collist]

com <- rbind(df1_summary, df2_summary_grouped, df2_summary_grouped1)
head(com)

write.table(com, coverSReOut, row.names = FALSE, sep = ',')
Sys.time()
                        
                 
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


print("finish!!!")









