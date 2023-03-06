# Adaptation of the DADA2 pipeline to the NYC Fish experiment.
#
# The starting point is a set of Illumina-sequenced paired-end fastq files
# that have been split (or demultiplexed) by sample and which have barcodes/adapters already removed.
# The end product will be a sequence table, analogous to the ubiquitous "OTU table",
# which records the number of times sample sequences were observed in each sample. The key
# difference between the output of DADA2 and standard OTU analyses is that DADA2 infers sample
# sequences exactly rather than clustering sequencines into fuzzy OTUs which hide and complicate biological variation.
#
# This relies on Bioconductor packages(http://bioconductor.org/) in addition to base R.
# The FastQ files use in this experiment are available in the SRA under
# BioProject PRJNA358446.
#
#
#BiocManager::install(c("dada2","ShortRead", "phyloseq", "Biostrings"))
    library(dada2)
    library(ShortRead)
    library(phyloseq)
    library(Biostrings)
    library(dplyr)

# path the where the where the fastq data lives
#path <- "data-raw/fastqs/"
path <- "/home/lodopore/stockle/"
fns <- list.files(path)

# Filtering and Trimming
# First we read in the file names for all the fastq files and do a little string manipulation
# to get lists of the forward and reverse fastq files in matched order:
fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs) # Sort should keep them paired in order
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)

# Perform filtering and trimming
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
              #truncLen=c(100,100),
              trimLeft=c(18,18),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
save(out, file="out.file")

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])



seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)


sum(seqtab.nochim)/sum(seqtab)

save(seqtab, file="seqtab.file")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)




# Make "otu_table"
seqs <- colnames(seqtab)
otab <- otu_table(seqtab, taxa_are_rows=FALSE)
colnames(otab) <- paste0("Seq_", seq(ncol(otab)))
#Write fastas to test file
writeFasta <- function(seqs, output) {

  seqsout <- mapply( function(idx, sequence) paste0(">Seq_",idx,"\n",sequence,"\n"),
                     seq(length(seqs)),
                     seqs)
  write(paste0(seqsout), file = output, sep = "")

}

seqs_for_blast <- DNAStringSet(seqs)
names(seqs_for_blast) <-  sapply(seq(length(seqs)),function(x) {paste0("Seq_",x)})

# Write the fasta sequences and the OTU table
write.table(otab,  file="otutable.txt")
writeFasta(seqs, "tax_sequences.fasta")





# getN <- function(x) sum(getUniques(x))

# track <- cbind(
#   out, 
# sapply(dadaFs, getN), 
# sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
# colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
# rownames(track) <- sample.names
# head(track)



# # for(i in seq_along(fnFs)) {
# #   fastqPairedFilter(
# #     paste0(path, c(fnFs[i], fnRs[i])), c(filtFs[i], filtRs[i]),
# #     maxN=0, 
# #     #maxEE=2, 
# #     #truncQ=2, 
# #     #trimLeft=c(10, 10), 
# #     #truncLen=c(100,100), 
# #     compress=TRUE, 
# #     verbose=TRUE)
# # }

# # Dereplication
# #derepFs <- lapply(filtFs, derepFastq, verbose=TRUE)
# #derepRs <- lapply(filtRs, derepFastq, verbose=TRUE)

# # Name the derep-class objects by the sample names
# sam_names <- sapply(strsplit(fnFs, "/"), tail, n=1)
# sam_names <- sapply(strsplit(sam_names, "_"), `[`, 1)
# names(derepFs) <- sam_names
# names(derepRs) <- sam_names

# # Sample Inference
# dadaFs <- dada(derepFs, err=inflateErr(tperr1,3), errorEstimationFunction=loessErrfun, selfConsist = TRUE)
# dadaRs <- dada(derepRs, err=inflateErr(tperr1,3), errorEstimationFunction=loessErrfun, selfConsist = TRUE)

# # Identify chimeras
# bimFs <- sapply(dadaFs, isBimeraDenovo, verbose=TRUE)
# bimRs <- sapply(dadaRs, isBimeraDenovo, verbose=TRUE)
# print(unname(sapply(bimFs, mean)), digits=2)
# print(unname(sapply(bimRs, mean)), digits=2)

# # Merge paired reads
# mergers <- mapply(mergePairs, dadaFs, derepFs, dadaRs, derepRs, SIMPLIFY=FALSE)

# # Remove chimeras
# mergers.nochim <- mapply(function(mm, bF, bR) mm[!bF[mm$forward] & !bR[mm$reverse],], mergers, bimFs, bimRs, SIMPLIFY=FALSE)

# # Constructing the sequence table
# seqtab <- makeSequenceTable(mergers.nochim)
# table(nchar(colnames(seqtab)))
# length(unique(substr(colnames(seqtab), 1, 100)))
# dim(seqtab)

# # Make "otu_table"
# seqs <- colnames(seqtab)
# otab <- otu_table(seqtab, taxa_are_rows=FALSE)
# colnames(otab) <- paste0("Seq_", seq(ncol(otab)))
# #Write fastas to test file
# writeFasta <- function(seqs, output) {

#   seqsout <- mapply( function(idx, sequence) paste0(">Seq_",idx,"\n",sequence,"\n"),
#                      seq(length(seqs)),
#                      seqs)
#   write(paste0(seqsout), file = output, sep = "")

# }

# seqs_for_blast <- DNAStringSet(seqs)
# names(seqs_for_blast) <-  sapply(seq(length(seqs)),function(x) {paste0("Seq_",x)})

# # Write the fasta sequences and the OTU table
# writeFasta(seqs,tax_sequences)
# write.table(otab,  file="data/otutable.txt")

