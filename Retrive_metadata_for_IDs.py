###############################################
##Dmitry Sutormin, 2023##
##Retrive metadata for a list of sample IDs##

#Takes table with metadata for Evrogen as a list of IDs.
#Takes a full table with the whole metadata (Atlas).
#Retrive metadata for the list of IDs.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#################
### Variables to be defined.
#################

#Input Evrogen table.
Evrogen_table_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Sequencing\Evrogen_sequencing\Docs_2023\Sequencing_May_2023\\Evrogen_covering_June_1_2023.xlsx"

#Input Atlas database.
Atlas_db_inpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Atlas_database\Atlas_database_table_dump\\Atlas_database_dump_19_06_2023.xlsx"

#Path to the output file.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Sequencing\Evrogen_sequencing\Docs_2023\Sequencing_May_2023\\Evrogen_June_1_2023_with_metadata.xlsx"


#################
### Functions and analysis.
#################


########
# Read Evrogen metadata, retrive IDs and type of sequencing.
########

def read_Evrogen_metadata(evrogen_inpath):
    
    Evrogen_data=pd.read_excel(evrogen_inpath, sheet_name="Лист1", header=0)
    ID_seq_type_tuples=Evrogen_data[['Название образца*', 'Тип секвенирования']]
    Amplic_16S_type='Ампликон V3-V4'
    Shotgun_type='Shot-gun секвенирование метагеномной ДНК'
    ID_seq_type_tuples_16S=ID_seq_type_tuples[ID_seq_type_tuples['Тип секвенирования']==Amplic_16S_type]
    print(ID_seq_type_tuples_16S)
    
    ID_seq_type_tuples_Shotgun=ID_seq_type_tuples[ID_seq_type_tuples['Тип секвенирования']==Shotgun_type]
    print(ID_seq_type_tuples_Shotgun)
    
    Amplic_16S_IDs_list=ID_seq_type_tuples_16S['Название образца*'].tolist()
    Shotgun_IDs_list=ID_seq_type_tuples_Shotgun['Название образца*'].tolist()
    
    print(len(Amplic_16S_IDs_list), len(Shotgun_IDs_list))
    
    return Amplic_16S_IDs_list, Shotgun_IDs_list


########
# Read the whole Atlas metadata, return metadata for the list of IDs.
########

def read_Atlas_db_return_metadata(atlas_db_inpath, Amplic_16S_IDs_list, Shotgun_IDs_list, output_path):
    
    Atlas_database=pd.read_excel(atlas_db_inpath, sheet_name="Ответы на форму (1)", header=0)
    
    Atlas_subset_16S=Atlas_database[Atlas_database['Sample_ID'].isin(Amplic_16S_IDs_list)]
    Atlas_subset_Shotgun=Atlas_database[Atlas_database['Sample_ID'].isin(Shotgun_IDs_list)]
    
    with pd.ExcelWriter(output_path) as writer:  
        Atlas_subset_16S.to_excel(writer, sheet_name='16S_samples')
        Atlas_subset_Shotgun.to_excel(writer, sheet_name='Shotgun_samples')    
    
    return


########
# Wrapper function.
########

def Wrapper(evrogen_inpath, atlas_db_inpath, output_path):
    
    #Read Evrogen metadata, return IDs.
    Amplic_16S_IDs_list, Shotgun_IDs_list=read_Evrogen_metadata(evrogen_inpath)
    
    #Read Atlas database, return metadata for the list of IDs.
    read_Atlas_db_return_metadata(atlas_db_inpath, Amplic_16S_IDs_list, Shotgun_IDs_list, output_path)
    
    return

Wrapper(Evrogen_table_path, Atlas_db_inpath, Output_path)