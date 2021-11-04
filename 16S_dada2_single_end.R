library(dada2); packageVersion("dada2")
path <- "path/to/file" 
list.files(path)
fnFs <- sort(list.files(path, pattern="_S.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
rdataFile <- "D:/bioinformatic/Spongy/Contigs_scaffolds/18S_sponges/2018_F_trimm/Spongy_2018.RData"
save.image(rdataFile)
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
save.image(rdataFile)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaFs[[1]]
save.image(rdataFile)
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
table(nchar(getSequences(seqtab)))
save.image(rdataFile)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
save.image(rdataFile)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
save.image(rdataFile)

#save information:frequency ASV from seqtab.nochim
seqtab.nochim.export <- t(seqtab.nochim)
rownames(seqtab.nochim.export) <- paste("ASV", seq(nrow(seqtab.nochim.export)))
write.table(seqtab.nochim.export, sep="\t", file = "D:/bioinformatic/Spongy/Contigs_scaffolds/18S_sponges/seqtabnochim_n.csv")

#save information: make fasta file with ASV sequences
library("ape")
fastaMx <- as.matrix(colnames(seqtab.nochim))
rownames(fastaMx) <- rownames(seqtab.nochim.export)
write.dna(fastaMx, file = "D:/bioinformatic/Spongy/Contigs_scaffolds/18S_sponges/seqtabnochim.fasta" , format="fasta", nbcol = -1, colw = 700)

##mmseq2 clustering and forming new OTU table

#return data to DADA2
seqtab.nochimOTU <- read.table(sep = "\t", file = "path/to/file", header = TRUE, check.names=FALSE, row.names=1)
seqtab.nochimOTU <- t(seqtab.nochimOTU)

#taxonomy assignment
taxa <- assignTaxonomy(seqtab.nochimNEW, "D:/bioinformatic/silva_nr_v1??_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "D:/bioinformatic/silva_species_assignment_v1??.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
dim(taxa)

#if required save taxonomy with and without nt
write.table(taxa, sep="\t", file = "path/to/file")
write.table(taxa.print, sep="\t", file = "path/to/file")


#save information:frequency OTU from seqtab.nochi
seqtab.nochim.exportOTU <- t(seqtab.nochimOTU)
rownames(seqtab.nochim.exportOTU) <- paste("OTU", seq(nrow(seqtab.nochim.exportOTU)))
write.table(seqtab.nochim.exportOTU, sep="\t", file = "path/to/file")










