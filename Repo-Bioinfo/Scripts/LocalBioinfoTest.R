# gzip -d *.gz to unzip compressed archive files on a server
# Add Cutadap work flow

#install the latest R version
install.packages("installr")

library(installr)

updateR()

#install phyloseq

R.Version()
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

library("phyloseq"); packageVersion("phyloseq")

#install and load other libraries needed

library (tidyverse)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")
library (dada2);packageVersion("dada2")

library (digest)

BiocManager::install("seqinr")
library (seqinr)
library (knitr)

BiocManager::install("kableExtra")
library (kableExtra); packageVersion("kableExtra")

#Zach for Mark:
BiocManager::install(c("ShortRead", "Biostrings"))
library(ShortRead)
library(Biostrings)
library(dplyr)

getwd()
path <- "C:/Users/Yuan.Liu/Desktop/Runs-combined"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#We start by visualizing the quality profiles of the forward reads:
plotQualityProfile(fnFs[1:2])

#Now we visualize the quality profile of the reverse reads:
plotQualityProfile(fnRs[1:2])

#Assign the filenames for the filtered fastq.gz files.
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#We'll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the maximum number of "expected errors" allowed in a read, which is a better filter than simply averaging quality scores.

#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(140,140),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

#Note: truncLen=(150,150) since the green line is stable since cycle 75 [Mark used 100,100]. Try maxEE=c(1,1). Add trimLeft=c(18, 18) to trim primers


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(140,140), trimLeft=c(18,18), 
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
head(out)

# Learning the error rates
errF <- learnErrors(filtFs, multithread=TRUE)  #?Why 6 samples were used for learning the error rates (Tutorial used 20)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)



#We are now ready to apply the core sample inference algorithm to the filtered and trimmed sequence data.

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)  #?? By default, the dada function processes each sample independently.Can try dada(..., pool=TRUE) or dada(...,pool="pseudo") later!
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Inspecting the returned dada-class object:

dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#???remove non-target-length sequences from your sequence table: keep 135-200bp for now
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 122:150]
table(nchar(getSequences(seqtab2)))
dim(seqtab2)
#removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
table(nchar(getSequences(seqtab.nochim)))
sum(seqtab.nochim)/sum(seqtab2)  #Here chimeras make up about xx% of the merged sequence variants, but when we account for the abundances of those variants we see they account for only about 4% of the merged sequence reads.
#Zach for Mark:
save(seqtab, file="seqtab.file")


#Track reads through the pipeline
#As a final check of our progress, we'll look at the number of reads that made it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)  #!Most of reads were removed as chimeric, upstream processing may need to be revisited. In almost all cases this is caused by primer sequences with ambiguous nucleotides that were not removed prior to beginning the DADA2 pipeline

#To calculate total filtered read #:
sum(track[,6])

write.csv(track,"track.csv") #This is to further take out low-yield NextSeq libraries in Excel.

#Zach for Mark:
# Make "otu_table"
seqs <- colnames(seqtab.nochim)
otab <- otu_table(seqtab.nochim, taxa_are_rows=FALSE)
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



#Assign Taxonomy
#On the file server, zip and keep the fasta file: [yliu@smbfish RefLib]$ gzip -k 2018mar6ECO_refs_for_DADA2updated2020_04_03apr.fas
#.fa is practically the same as .fas, so I simply changed Mark's file name to *.fa.
# Tax information in Mark's ref library is not formmatted to be compatible with assignTaxonomy though. 

taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/Yuan.Liu/Desktop/Runs-combined/tax/Ref.fa.gz", multithread=TRUE)

