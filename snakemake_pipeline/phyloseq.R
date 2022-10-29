### good phyloseq script
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library("tibble") 
library(RColorBrewer)

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  path <- args[1]
}

directory = main()

otu_mat = read.csv2(file = paste(directory, '/results/all_OTU_frequency.tsv', sep=''),
                    header=TRUE, check.names = F, sep='\t')
tax_mat = read.csv2(file = paste(directory, '/results/all_phylogeny.tsv', sep=''),
                    header=TRUE, check.names = F, sep='\t')

tax_mat[1] = paste("OTU_", seq(nrow(tax_mat)), sep = '')

sample_names <- colnames(otu_mat)

print(sample_names)

colnames(otu_mat) = c('OTU', tail(sample_names, -1))


OTU_name = tax_mat[1]
colnames(OTU_name) = 'OTU'

rownames(tax_mat) = OTU_name$OTU
tax_mat[1] = NULL
rownames(otu_mat) = OTU_name$OTU
otu_mat[1] = NULL

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)


### get phyloseq object
carbom <- phyloseq(OTU, TAX)
carbom

sample_names(carbom)
rank_names(carbom)


### Normalize number of reads in each sample using median sequencing depth
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)


### plots figures
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.05) > 0, TRUE)
carbom_abund

total = median(sample_sums(carbom_abund))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom_abund = transform_sample_counts(carbom_abund, standf)


fig_1 <- plot_bar(carbom_abund, fill = "Class") + 
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") + 
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20),
        axis.text = element_text(size = 8), axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 65, hjust = 1), axis.text.y = element_text(size = 12)) + 
  scale_y_continuous(expand = c(0.03,0.03), breaks = c(0, 13000, 26000, 39000, 52000), 
                     labels = c('0', '25', '50', '75', '100'))
ggsave(paste(directory, "/results/phyloseq.png", sep=''))

fig_2 <- plot_heatmap(carbom_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Class", taxa.order = "Class", 
             trans=NULL, low="beige", high="red", na.value="beige") + 
  theme(legend.title = element_text(size=20), legend.text = element_text(size=20),
        axis.text.x = element_text(size = 6, angle = 65, vjust = 1, hjust = 1), 
        axis.title = element_text(size = 20))
ggsave(paste(directory, "/results/heat_map.png", sep=''))

