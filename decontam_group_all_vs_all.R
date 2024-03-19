library(phyloseq); packageVersion("phyloseq")
#library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library(jsonlite)


### Prepare the data

# metadata
samples_all = read.csv2('work/results_21_10_22/Metadata_DA.tsv',
                        header=TRUE, check.names = F, sep='\t', row.names = 1)

# otu tables
otu_mat = read.csv2('work/results_21_10_22/OTU_new/all_OTU_frequency.tsv',
                    header=TRUE, check.names = F, sep='\t')


#rename otu_table
sample_id = colnames(otu_mat)
id_names = c()
for (i in sample_id){
  id_names = c(id_names, substr(i, 2, 6))
}
colnames(otu_mat) = id_names

### taxonomy table
tax_mat = read.csv2('work/results_21_10_22/OTU_new/all_phylogeny.tsv',
                    header=TRUE, check.names = F, sep='\t')
tax_mat[1] = paste("OTU_", seq(nrow(tax_mat)), sep = '')

# filter metadata with IDs 
sample_id <- colnames(otu_mat)
samples_df = data.frame()
for (i in sample_id){
  print(i)
  row = samples_all[samples_all$Sample_ID == i, ]
  samples_df = rbind(samples_df, row)
}

# rename rownames in otu and tax tables
OTU_name = tax_mat[1]
colnames(OTU_name) = 'OTU'
rownames(tax_mat) = OTU_name$OTU
tax_mat[1] = NULL
rownames(otu_mat) = OTU_name$OTU
otu_mat[1] = NULL


otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE) # otu table for phyloseq object
TAX = tax_table(tax_mat) # tax table for phyloseq object
samples = sample_data(samples_df) # metadata for phyloseq object



### some comments for me
#contam_otu = data.frame(matrix(0, nrow=nrow(as.data.frame(otu_table(carbom))), 
                               #ncol=ncol(as.data.frame(otu_table(carbom)))))
#colnames(contam_otu) = colnames(as.data.frame(otu_table(carbom)))
#rownames(contam_otu) = rownames(as.data.frame(otu_table(carbom)))


animal_soil = subset(samples, samples$Expedition.station.number == "76" | 
                       samples$Expedition.station.number == "M2")
#methods = c("either", "prevalence")
#threshold_array = c(0.1, 0.2, 0.3, 0.4, 0.5)

methods = c("either")
threshold_array = c(0.5)

experiment_list = list()
experiment_list[[1]] = animal_soil$Sample_ID
experiment_list[[2]] = water_all$Sample_ID

for (method in methods){
  print(paste("Start method:", method))
  for (threshold in threshold_array){
    print(paste("Use threshold:", threshold))
    
    output_dir = "work/results_21_10_22/" # set your directory
    count = 0
    start = 1
    end = 0
    
    #print(length(s_list))
    #print(s_list)

    sample_ids = c()
    
    contamin_list = list()
    
    for (experiment in experiment_list){
      #print(control)
      print(experiment)
      count = count + 1
      end = end + length(experiment) - 8
      
      
      current_set = subset(samples, samples$Sample_ID %in% experiment)
      print(paste("Current set:", count))
      
      sample_set = subset(current_set, current_set$Sample_type != "control")$Sample_ID
      print(sample_set)
      
      
      current_samples = current_set
      carbom <- phyloseq(OTU, TAX, current_samples)
      
      sample_data(carbom)$is.neg <- sample_data(carbom)$Sample_type == "control";
      #sample_data(carbom)
      
      ### adjust arguments in decontam
      contamdf.either <- isContaminant(carbom, method=method, neg="is.neg",
                                       conc="DNA.concentration..ng.mkl",
                                       threshold=threshold, normalize = T)
      
      print(table(contamdf.either$contaminant))
      #print(which(contamdf.either$contaminant))
      
      sample_ids = c(sample_ids, sample_set)
      
      for (j in rep(start:end)){
        contamin_list[[j]] = c(which(contamdf.either$contaminant))
        print(j)
      }
      
      start = start + end
    }
    
    
    #print(contamin_list)
    #print(sample_ids)
    
    ### save the contaminant OTUs in json file
    names(contamin_list) = sample_ids
    json_string <- toJSON(contamin_list, auto_unbox = TRUE, pretty = TRUE)
    write(json_string, paste0(output_dir,"/contamin_list_",method,"_",as.character(threshold),"_grouped_all_vs_all_2.json"))
    print("Save JSON...")
    
    
    ### filter the contaminant OTUs
    carbom <- phyloseq(OTU, TAX, all_samples)
    carbom
    
    otu_df = as.data.frame(otu_table(carbom))
    
    
    for (column in names(contamin_list)){
      #print(column)
      for (contamin in contamin_list[[column]]){
        otu_df[contamin, column] = 0 #set zero frequency for a contaminant OTU
      }
    }
    
    
    ### save otu table
    write.table(otu_df, 
                paste0(output_dir,"/OTUs_decontam_",method,"_",as.character(threshold),"_grouped_all_vs_all_2.tsv"), 
                sep = "\t", row.names = TRUE, col.names = TRUE)
    print("Save filtered OTU table...")
    
    otu_df = as.data.frame(otu_table(carbom))
    
    for (column in names(contamin_list)){
      #print(column)
      for (contamin in contamin_list[[column]]){
        contam_otu[contamin, column] = otu_df[contamin, column] #set zero frequency for a contaminant OTU
      }
    }
    
    
    ### save otu table
    write.table(contam_otu, 
                paste0(output_dir,"/OTUs_contamin_",method,"_",as.character(threshold),"_grouped_all_vs_all.tsv"), 
                sep = "\t", row.names = TRUE, col.names = TRUE)
    print("Save contaminated OTU table...")
    
  }
}



#########################
#   THE END OF SCRIPT   #
#########################


