# Grubisic et al. < add ms title >

# R script for DADA2 analysis of DNA metabarcoding data 

# Script modified from the DADA2 Pipeline Tutorial (1.16) at https://benjjneb.github.io/dada2/tutorial.html

# load dada2 package

library(dada2); packageVersion("dada2")

# set path

path <- "~/data/Grubisic_et_al.data" 
list.files(path)

# read in file names and obtain matched lists

fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

# check read quality

plotQualityProfile(fnFs[1:12])

plotQualityProfile(fnRs[1:12])

# filter and trim
# check that negative controls were good (reads.out 0) and re-run this step after removing

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(20,20), truncLen=c(260,200), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
head(out)

# error rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# dereplication

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]

# merge paired reads

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

# sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

# chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# check and track reads

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# extract asv table with sequence and abundance per sample

write.table(t(seqtab.nochim), "seqtab_nochim.tsv")

# assign taxonomy

saveRDS(seqtab.nochim, "seqtab.nochim.RDS")

taxa <- assignTaxonomy(seqtab.nochim, "~/path/pr2_version_4.12.0_18S_dada2.fasta.gz")

saveRDS(taxa, "taxa.RDS")

# write to csv

write.csv(cbind(t(seqtab.nochim), taxa), "~/path/seqtab_nonchim_id.csv", quote=FALSE)

