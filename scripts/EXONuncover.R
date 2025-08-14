

if(!require("data.table")) {install.packages("data.table")}
suppressMessages(library(data.table))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript script.R <mainPath> <sampleName> <Aligner> <refEXON>")
}

mainPath <- args[1]
sampleName <- args[2]
Aligner <- args[3]
refEXON <- args[4]

# Output file paths
intersectBedFile <- file.path(mainPath, Aligner, paste0(sampleName, "_mapped_woSecond_intersectS.bed"))
bedFile <- file.path(mainPath, Aligner, paste0(sampleName, "_mapped_woSecond.bed"))
outputFile <- file.path(mainPath, Aligner, paste0(sampleName, "_mapped_woSecond_intersectS_EXONuncover.bed"))

# Load reference CDS
cat("Loading reference exon file:", refEXON, "\n")
cds <- fread(refEXON, col.names = c("chr", "start", "end", "name", "score", "strand"))
cds[, isoform := tstrsplit(name, "__")[[1]]]

# Function to find uncovered BED blocks
extract_uncovered_bed <- function(interFile, bedFile, outFile) {
  cat("Reading intersect BED file:", interFile, "\n")
  inter <- fread(interFile, header = FALSE)
  setnames(inter, c("chr", "start", "end", "sampleID", "score", "strand",
                    "CDS_chr", "CDS_start", "CDS_end", "CDS_name", "CDS_score", "CDS_strand", "n_base"))

  # Identify blocks not overlapping exons
  uncovered_keys <- inter[CDS_name == ".", 
                          paste(chr, start, end, sampleID, sep = ":")]

  cat("Reading original BED file:", bedFile, "\n")
  bed <- fread(bedFile, header = FALSE)

  # Create match keys
  bed_keys <- paste(bed[[1]], bed[[2]], bed[[3]], bed[[4]], sep = ":")

  # Match and extract uncovered reads
  uncovered_bed <- bed[bed_keys %in% uncovered_keys]

  # Write uncovered reads to output
  fwrite(uncovered_bed, file = outFile, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  cat("Uncovered blocks saved to:", outFile, "\n")
}

# Run the extraction
extract_uncovered_bed(intersectBedFile, bedFile, outputFile)

