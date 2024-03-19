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

### subset all controls and samples
all_samples = subset(samples, (samples$Sample_type == "sample"))
all_control = subset(samples, (samples$Sample_type == "control"))

###########################
### Script for decontam ###
###########################


methods = c("either", "prevalence") # choose the methods for decontam
threshold_array = c(0.1, 0.2, 0.3, 0.4, 0.5) # set a list with threshold

for (method in methods){
  print(paste("Start method:", method))
  for (threshold in threshold_array){
    print(paste("Use threshold:", threshold))
    
    output_dir = "work/results_21_10_22/OTU_new/" #set your output directory
    
    ### create arrays with indexes for controls and samples
    control_list = c(all_control$Sample_ID)
    control_list = c(control_list, control_list[1:8])
    #print(control_list)
    
    s_list = c(all_samples$Sample_ID)
    #print(length(s_list))
    
    sample_ids = c()
    
    contamin_list = list()
    
    count = rep(0:(length(control_list)-1))
    #print(count)
    
    ### receive contaminant OTUs using decontam
    for (i in count) {
      print(paste("Filter a group of samples: ", i+1))
      control = as.character(control_list[i+1])
      sample_a = as.character(s_list[3*i + 1])
      sample_b = as.character(s_list[3*i + 2])
      sample_c = as.character(s_list[3*i + 3])
      #print(paste(control, sample_a, sample_b, sample_c))
      
      name_list = c(sample_a, sample_b, sample_c)
      current_set = subset(samples, samples$Sample_ID == control | samples$Sample_ID == sample_a |
                             samples$Sample_ID == sample_b | samples$Sample_ID == sample_c)
      
      current_samples = current_set
      carbom <- phyloseq(OTU, TAX, current_samples)
      #carbom
      
      sample_data(carbom)$is.neg <- sample_data(carbom)$Sample_type == "control";
      #sample_data(carbom)
      
      ### adjust arguments in decontam
      contamdf.either <- isContaminant(carbom, method=method, neg="is.neg",
                                       conc="DNA.concentration..ng.mkl",
                                       threshold=threshold, normalize = T)
      
      table(contamdf.either$contaminant)
      #print(which(contamdf.either$contaminant))
      
      sample_ids = c(sample_ids, name_list)
      
      contamin_list[[3*i+1]] = c(which(contamdf.either$contaminant))
      contamin_list[[3*i+2]] = c(which(contamdf.either$contaminant))
      contamin_list[[3*i+3]] = c(which(contamdf.either$contaminant))
      
    }
    
    #print(contamin_list)
    #print(sample_ids)
    
    ### save the contaminant OTUs in json file
    names(contamin_list) = s_list
    json_string <- toJSON(contamin_list, auto_unbox = TRUE, pretty = TRUE)
    write(json_string, paste0(output_dir,"/contamin_list_",method,"_",as.character(threshold),".json"))
    print("Save JSON...")
    
    
    ### filter the contaminant OTUs
    carbom <- phyloseq(OTU, TAX, all_samples)
    carbom
    
    otu_df = as.data.frame(otu_table(carbom))
    
    
    for (column in rep(1:length(contamin_list))){
      #print(column)
      for (contamin in contamin_list[column]){
        for (row in contamin){
          otu_df[contamin, column] = 0 #set zero frequency for a contaminant OTU
        }
      }
    }
    
    
    
    
    ### save otu table
    write.table(otu_df, 
                paste0(output_dir,"/OTU_frequency_decontam_",method,"_",as.character(threshold),"_by_groups.tsv"), 
                sep = "\t", row.names = TRUE, col.names = TRUE)
    print("Save filtered OTU table...")
    
  }
}


#########################
#   THE END OF SCRIPT   #
#########################


