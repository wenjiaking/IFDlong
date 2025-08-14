

suppressPackageStartupMessages({
  install_and_load <- function(pkgs) {
    for (p in pkgs) {
      if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
      library(p, character.only = TRUE)
    }
  }
  install_and_load(c("rlist", "stringr", "parallel", "data.table", "dplyr"))
})


#### parameter
args = commandArgs(trailingOnly=TRUE)

PATH <- args[1]
sampleName <- args[2]
Aligner <- args[3]
buffer <- as.integer(args[4])
anchorLen <- as.integer(args[5])
refGTF <- args[6] #refData/genes.gtf"
ncores=args[7]


reportfile=paste0(PATH,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_Rep.csv")
Isof.quantfile=paste0(PATH,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_Isof_quant.csv")
Isof.quantRData=paste0(PATH,"/",Aligner,"/",sampleName,"_buffer",buffer,"bp_temp_Isof_quant.RData")
fusionfiltfile=paste0(PATH,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_fusionRep_anchor",anchorLen,"bp.filt.csv")
fusion.quantfile=paste0(PATH,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_Fusion_quant_anchor",anchorLen,"bp.csv")
fusion.quantRData=paste0(PATH,"/",Aligner,"/",sampleName,"_buffer",buffer,"bp_temp_Fusion_quant_anchor",anchorLen,"bp.RData")

############ load 
Report <- read.csv(reportfile)
fusionRep <- read.csv(fusionfiltfile)
gtf.dat <- read.table(refGTF, fill = TRUE)



EMperGene <- function(read_annot, tol=1e-5, max.iter=200) {
  isoforms <- tryCatch({
    strsplit(unlist(strsplit(read_annot$isoform, "#")), "\\|\\|")
  }, error = function(e) {
    if (grepl("non-character argument", e$message)) {
      message("Warning: 'isoform' column contains non-character data. Skipping this chunk.")
      return(NULL)
    } else {
      stop(e)  # re-throw other errors
    }
  })
  
  if (is.null(isoforms)) {
    return(NULL)  # Skip processing if error happened
  }
  
  theta0 <- prop.table(table(unlist(isoforms)))
  
  if (!("undefined" %in% names(theta0))) {
    theta0["undefined"] <- 0
  }
  
  nr <- nrow(read_annot)
  mat_update <- matrix(0, nrow=nr, ncol=length(theta0))
  colnames(mat_update) <- names(theta0)
  
  theta <- rep(0, length(theta0))
  names(theta) <- names(theta0)
  
  niter <- 0
  
  while (sqrt(sum((theta0 - theta)^2)) > tol && niter < max.iter) {
    if (niter > 0) theta0 <- theta
    niter <- niter + 1
    
    for (i in 1:nr) {
      r <- read_annot[i, ]
      
      zknown <- (!grepl("#", r[1, 4])) & (!grepl("\\|\\|", r[1, 4]))
      notes = tryCatch({
        unlist(strsplit(as.character(r[1, 11]), "#"))
      }, error = function(e) {
        message("Warning: non-character data in r[1,11], skipping this row")
        return(NA_character_)
      })
      
      if (any(is.na(notes))) {
        next
      }
      
      if (any(notes == "continuous CDS and edge-matching")) {
        if (zknown) {
          mat_update[i, r[1, 4]] <- 1
        } else {
          isofs <- unlist(strsplit(unlist(strsplit(r[1, 4], "#")), "\\|\\|"))
          # Normalize theta0 over isofs to avoid NA if sum==0
          s <- sum(theta0[isofs])
          if (s > 0) {
            mat_update[i, isofs] <- theta0[isofs] / s
          } else {
            mat_update[i, isofs] <- 1 / length(isofs)
          }
        }
      } else {
        mat_update[i, "undefined"] <- 1
      }
    }
    
    theta <- colSums(mat_update, na.rm = TRUE) / nr
  }
  
  return(list(prop = theta, counts = theta * nr, niter = niter))
}

quantBygene <- function(report.multigenes, report.uniquegene, tol=1e-5, max.iter=200, mc=30) {
  require(parallel)
  
  readBYgeneList <- list()
  n_multigenes <- length(unique(report.multigenes$gene))
  
  quantBYgenesList <- mclapply(1:n_multigenes, function(i) {
    g <- unique(report.multigenes$gene)[i]
    genes <- unlist(strsplit(g, "#"))
    
    read.index <- which(report.uniquegene$gene %in% genes)
    if (length(read.index) != 0) {
      reads <- rbind(report.multigenes[report.multigenes$gene == g, ], report.uniquegene[read.index, ])
    } else {
      reads <- report.multigenes[report.multigenes$gene == g, ]
    }
    
    quant <- EMperGene(read_annot = reads, tol = tol, max.iter = max.iter)
    return(quant)
  }, mc.cores = mc)
  
  names(quantBYgenesList) <- unique(report.multigenes$gene)
  
  multigenes <- unique(unlist(strsplit(names(quantBYgenesList), "#")))
  report.uniquegene <- report.uniquegene[!report.uniquegene$gene %in% multigenes, ]
  
  quantBYgeneList <- mclapply(1:length(unique(report.uniquegene$gene)), function(j) {
    g <- unique(report.uniquegene$gene)[j]
    reads <- report.uniquegene[report.uniquegene$gene == g, ]
    quant <- EMperGene(read_annot = reads, tol = tol, max.iter = max.iter)
    return(quant)
  }, mc.cores = mc)
  
  names(quantBYgeneList) <- unique(report.uniquegene$gene)
  
  return(c(quantBYgenesList, quantBYgeneList))
}


EMperGenefusion = function(read_annot, tol=1e-5, max.iter=200) {
  genes = tryCatch({
    unlist(strsplit(read_annot$isoform, "#"))
  }, error = function(e) {
    if (grepl("non-character argument", e$message)) {
      message("Warning: 'isoform' column contains non-character data. Skipping this chunk.")
      return(NULL)  # or return NA or empty vector if preferred
    } else {
      stop(e)  # re-throw other errors
    }
  })
  
  if (is.null(genes)) {
    return(NULL)  # Skip processing if error happened
  }
  
  isoforms = unlist(as.vector(sapply(genes, function(g) {
    pair = unlist(strsplit(g, "&"))
    pairlist = lapply(pair, function(x) unlist(strsplit(x, "\\|\\|")))
    head = pairlist[[1]]
    for (j in 2:length(pairlist)) {
      head = as.vector(sapply(head, function(x) paste0(x, "&", pairlist[[j]])))
    }
    return(head)
  })))
  isoforms[grepl("undefined", isoforms)] = "undefined"
  isoforms[grepl("NA", isoforms)] = "undefined"
  theta0 = prop.table(table(unlist(isoforms)))
  if (!("undefined" %in% names(theta0))) {
    theta0["undefined"] = 0
  }
  nr = nrow(read_annot)
  mat_update = matrix(0, nrow=nr, ncol=length(theta0))  # initialize with zeros
  colnames(mat_update) = names(theta0)
  theta = Inf
  niter = 0
  
  while(sqrt(sum((as.vector(theta0) - as.vector(theta))^2)) > tol & (niter < max.iter)) {
    if (niter > 0) { theta0 = theta }
    niter = niter + 1
    for (i in 1:nr) {
      r = read_annot[i,]
      zknown = (!grepl("#", r[1,4])) & (!grepl("\\|\\|", r[1,4]))
      notes = unlist(strsplit(r[1,11], "#"))
      ng = str_count(notes[1], "&")
      if (any(notes == paste0(rep("continuous CDS and edge-matching", ng + 1), collapse = "&") & !is.na(notes))) {
        if (zknown) {
          mat_update[i, r[1,4]] = 1
        } else {
          isofmulti = unlist(strsplit(r[1,4], "#"))
          isofs = unlist(as.vector(sapply(isofmulti, function(g) {
            pair = unlist(strsplit(g, "&"))
            pairlist = lapply(pair, function(x) unlist(strsplit(x, "\\|\\|")))
            head = pairlist[[1]]
            for (j in 2:length(pairlist)) {
              head = as.vector(sapply(head, function(x) paste0(x, "&", pairlist[[j]])))
            }
            return(head)
          })))
          isofs[grepl("undefined", isofs)] = "undefined"
          isofs[grepl("NA", isofs)] = "undefined"
          isofs = unique(isofs)
          mat_update[i, isofs] = theta0[isofs] / sum(theta0[isofs])
        }
      } else {
        mat_update[i, "undefined"] = 1
      }
    }
    theta = apply(mat_update, 2, function(t) sum(t, na.rm = TRUE) / length(t))
  }
  return(list(prop = theta, counts = theta * nr, niter = niter))
}

quantBygenefusion=function(report.multigenes,report.uniquegene,tol=1e-5,max.iter=200,mc=30) {
  readBYgeneList=list()
  n_multigenes=length(unique(report.multigenes$gene))
  
  quantBYgenesList=mclapply(1:n_multigenes, function(i) {
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
  },mc.cores = mc)
  
  names(quantBYgenesList)=unique(report.multigenes$gene)
  multigenes=unique(unlist(strsplit(names(quantBYgenesList),"#")))
  if (length(which(report.uniquegene$gene %in% multigenes))!=0) {
    report.uniquegene=report.uniquegene[-which(report.uniquegene$gene %in% multigenes),]
  }
  
  
  quantBYgeneList=mclapply(1:length(unique(report.uniquegene$gene)),function(j) {
    g=unique(report.uniquegene$gene)[j]
    reads=report.uniquegene[report.uniquegene$gene==g,]
    quant=EMperGenefusion(read_annot = reads,tol=tol,max.iter=max.iter)
    return(quant)
  },mc.cores = mc)
  
  names(quantBYgeneList)=unique(report.uniquegene$gene)
  
  return(c(quantBYgenesList,quantBYgeneList))
}


quantList.gen <- function(Report, fusionRep, Isof.quantfile, Isof.quantRData, fusion.quantfile, fusion.quantRData, gtf.dat, tol=1e-5, max.iter=200, mc=30) {
  isof.ind <- which(Report$fusion == "N")
  
  if (length(isof.ind) != 0) {
    Report.isof <- Report[Report$fusion == "N", ]
    multigenes.ind <- grep("#", Report.isof$gene)
    Report.multigenes <- Report.isof[multigenes.ind, ]
    Report.uniquegene <- Report.isof[-multigenes.ind, ]
    
    quantList <- quantBygene(report.multigenes = Report.multigenes, report.uniquegene = Report.uniquegene, tol = tol, max.iter = max.iter, mc = mc)
    
    # list.save(quantList, Isof.quantRData) # Uncomment if using rlist package or equivalent
    
    count <- c()
    prop <- c()
    group <- c()
    isoform <- c()
    
    for (i in seq_along(quantList)) {
      count <- c(count, quantList[[i]]$counts)
      isoform <- c(isoform, names(quantList[[i]]$prop))
      prop <- c(prop, quantList[[i]]$prop)
      group <- c(group, rep(names(quantList)[i], length(quantList[[i]]$counts)))
    }
    
    gene <- mclapply(isoform, function(x) unique(gtf.dat$V13[gtf.dat$V16 == x | gtf.dat$V19 == x]), mc.cores = mc)
    gene[sapply(gene, length) == 0] <- "undefined"
    
    isof_quant.dat <- data.frame(isoform = isoform, gene = unlist(gene), group = group, prop = prop, count = count)
    write.csv(isof_quant.dat, file = Isof.quantfile, row.names = FALSE)
    
  } else {
    isof_quant.dat <- NULL
    message("No isoform detected in the data")
  }
  
  if (is.data.frame(fusionRep)) {
    ######## Fusion quantification ################
    multigenes.ind=grep("#",fusionRep$gene)
    fusionRep.multigenes=fusionRep[multigenes.ind,]
    fusionRep.uniquegene=fusionRep[-multigenes.ind,]
    
    quantList=quantBygenefusion(report.multigenes=fusionRep.multigenes,report.uniquegene=fusionRep.uniquegene,tol=tol,max.iter=max.iter,mc=mc)
    
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
    
    gene=mclapply(isoform,function(x) {
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
    },mc.cores = mc)
    gene=unlist(gene)
    fusion_quant.dat=data.frame(isoform=isoform,gene=gene, group=group,prop=prop,count=count)
    
    if (length(isof.ind)!=0) {
      isof_counts=mclapply(fusion_quant.dat$isoform,function(x) {
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
      },mc.cores = mc)
      isof_counts=unlist(isof_counts)
      gene_counts=mclapply(fusion_quant.dat$gene,function(x) {
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
      },mc.cores = mc)
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
  
 
  #return(list(isof_quant = isof_quant.dat))
  return(list(fusion_quant=fusion_quant.dat,isof_quant=isof_quant.dat))
}


#### test
quantList=quantList.gen(Report,fusionRep,Isof.quantfile,Isof.quantRData, fusion.quantfile,fusion.quantRData,gtf.dat,tol=1e-5,max.iter=200,mc=ncores)

# 
# #### Split Report into 20 Parts
# split_report_parts <- Report %>%
#   group_by(gene)%>%
#   group_split() %>%
#   { split(., cut(seq_along(.), 20, labels = FALSE)) } %>%
#   lapply(dplyr::bind_rows)
# 
# 
# #### Process Each Report Chunk then merged
# merged_isof_quant <- list()
# merged_fusion_quant <- list()
# 
# for (m in seq_along(split_report_parts)) {
#   message("Processing chunk: ", m)
#   Report_chunk <- split_report_parts[[m]]
#   
#   Isof.quantfile_chunk <- paste0(PATH, "/", Aligner, "/", sampleName, "_", m, "_mapped_woSecond_intersectS_buffer", buffer, "bp_Isof_quant.csv")
#   Isof.quantRData_chunk <- paste0(PATH, "/", Aligner, "/", sampleName, "_", m, "_buffer", buffer, "bp_temp_Isof_quant.RData")
#   
#   result <- quantList.gen(Report_chunk, 
#                           fusionRep, Isof.quantfile, Isof.quantRData, 
#                           fusion.quantfile, fusion.quantRData, gtf.dat, tol=1e-5,
#                           max.iter=200, mc=ncores)
# 
#   
#   if (!is.null(result$isof_quant)) {
#     merged_isof_quant[[m]] <- result$isof_quant
#   }
#   if (!is.null(result$fusion_quant)) {
#     merged_fusion_quant[[m]] <- result$fusion_quant
#   }
# }
# 
# # Combine all chunks into one data frame each
# final_isof_quant <- do.call(rbind, merged_isof_quant)
# final_fusion_quant <- do.call(rbind, merged_fusion_quant)
# 
# # Optionally, write to CSV
# write.csv(final_isof_quant, paste0(PATH, "/", Aligner, "/", sampleName, "_merged_Isof_quant.csv"), row.names = FALSE)
# write.csv(final_fusion_quant, paste0(PATH, "/", Aligner, "/", sampleName, "_merged_Fusion_quant.csv"), row.names = FALSE)
# 
# 
# 
# 