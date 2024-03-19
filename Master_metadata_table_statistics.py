###############################################
##Dmitry Sutormin, 2023##
##Calculates statistics for master metadata table##

#Takes master table with the whole metadata (Atlas).
#Calculates statistics.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openpyxl import load_workbook

#################
### Variables to be defined.
#################

#Input Atlas database.
Atlas_db_inpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Atlas_database\Atlas_database_table_dump\\Samples_for_shot_gun_sequencing_09_03_2024.xlsx"

#Path to the output file.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Atlas_database\Atlas_database_table_dump\Atlas_statistics\\"

if not os.path.isdir(Output_path):
    os.mkdir(Output_path)
    
    
def write_sheet_to_excel(excel_path, df_to_write, df_sheet_name):
    
    if not os.path.isfile(excel_path):
        
        writer=pd.ExcelWriter(excel_path, engine='xlsxwriter')
        df_to_write.to_excel(writer, sheet_name=df_sheet_name, index=False)
        writer.close()         
        
    else:
    
        book=load_workbook(excel_path)
        writer=pd.ExcelWriter(excel_path, engine = 'openpyxl')
        writer.book=book
        df_to_write.to_excel(writer, sheet_name=df_sheet_name, index=False)
        writer.close()    
    
    return


def expeditions_sample_types_barplot(Atlas_subset_main_sample, Color_dict, excel_path):
    
    plot_path=f'{excel_path[:-5]}_expeditions_stacked_barplot'
    
    List_of_expeditions=list(set(Atlas_subset_main_sample['Expedition_ID'].tolist()))
    
    Num_samples_ar=[]
    
    for expedition_name in List_of_expeditions:
        number_of_samples=Atlas_subset_main_sample[Atlas_subset_main_sample['Expedition_ID']==expedition_name].shape[0]
        Num_samples_ar.append(number_of_samples)
        
    Num_samples_ar_sorted, List_of_expeditions_sorted = (list(t) for t in zip(*sorted(zip(Num_samples_ar, List_of_expeditions), reverse=True))) # Sorting method taken from here: https://stackoverflow.com/questions/9764298/given-parallel-lists-how-can-i-sort-one-while-permuting-rearranging-the-other
    
    
    fig, ax = plt.subplots(figsize=(6.8, len(List_of_expeditions_sorted)*0.2))
    
    for expedition_name in List_of_expeditions_sorted:
        
        Expedition_data=Atlas_subset_main_sample[Atlas_subset_main_sample['Expedition_ID']==expedition_name]
        
        samples_num_cumul=0
        
        if Expedition_data.shape[0]>0:
            
            sample_types_ar_sorted=sorted(list(set(Expedition_data['Sample_type_extended'].tolist())), reverse=True)
        
            for sample_type in sample_types_ar_sorted:
                
                sample_number=Expedition_data[Expedition_data['Sample_type_extended']==sample_type].shape[0]
                
                if sample_number>1:
                    ax.barh(expedition_name, sample_number, 0.65, color=Color_dict[sample_type], left=samples_num_cumul)
                else:
                    ax.barh(expedition_name, sample_number, 0.65, color=Color_dict[sample_type], left=samples_num_cumul)            
                samples_num_cumul+=sample_number
                
        else:
            
            ax.bar(expedition_name, 0, 0.65)
    
    color_code_cumul=0
    
    for sample_type in Color_dict.keys():
        ax.barh('Color code', 1, 0.65, color=Color_dict[sample_type], left=color_code_cumul, label=sample_type)
        color_code_cumul+=1
    
    ax.set_ylim([-0.5, len(List_of_expeditions_sorted)-0.5])
    ax.tick_params(axis='x', which='major', labelsize=10)
    plt.xticks(rotation=90)
    ax.tick_params(axis='y', which='major', labelsize=8)
    ax.set_xlabel('Number of samples', size=12)
    #ax.set_ylabel('Expeditions', size=12)
    ax.set_title('Sequenced samples statistics', size=14)
    ax.spines['top'].set_color(None)
    ax.spines['right'].set_color(None)
    ax.invert_yaxis()
    
    #Place legend outside of a graph. Taken from: https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
    #box=ax.get_position()
    #ax.set_position([box.x0, box.y0, box.width * 0.4, box.height])
    ax.legend(fontsize=10, ncol=1, loc='lower right', frameon=False, markerscale=1, handlelength=0.7, handletextpad=0.3, columnspacing=0.7, labelcolor='k') #bbox_to_anchor=(1, 0.95), 
    plt.tight_layout() #rect=[0,0,0.7,1]
    plt.savefig(f'{plot_path}.png', dpi=300)
    plt.savefig(f'{plot_path}.svg', dpi=300)
    plt.close()
    
    return


def prep_points_for_Qgis(Atlas_subset_main_sample, excel_path):
    
    points_path=f'{excel_path[:-5]}_Qgis_points.csv'
    fileout=open(points_path, 'w')
    
    List_of_expeditions=list(set(Atlas_subset_main_sample['Expedition_ID'].tolist()))
    
    Num_samples_ar=[]
    
    fileout.write('Lat\tLon\tName\n')
    
    for expedition_name in List_of_expeditions:  
        
        Expedition_data=Atlas_subset_main_sample[Atlas_subset_main_sample['Expedition_ID']==expedition_name]
        
        coordinates_list=list(set(Expedition_data['Coordinates, decimal'].tolist()))
        
        start=0
        
        for coord_pair in coordinates_list:
            print(coord_pair)
            coord_pair=coord_pair.split(' N, ')
            
            if start==0:
                fileout.write(f'{coord_pair[0]}\t{coord_pair[1][:-2]}\t{expedition_name} ({Expedition_data.shape[0]})\n')
                
            else:
                fileout.write(f'{coord_pair[0]}\t{coord_pair[1][:-2]}\n')
            
            start+=1
    
    fileout.close()
    
    return

    
def calc_write_statistics(Atlas_subset_main_sample, action, excel_path):
    
    Color_dict={'soil, soil' : '#755a30', 'soil, rhizosphere' : '#241c0f', 
                'sea, water' : '#3ba6e0', 'sea, soil' : '#9898a1', 'sea, seaweed' : '#e85568', 'sea, plant' : '#2be322', 'sea, organism' : '#ffbe31', 
                'plant, plant' : '#69992e', 'laboratory control, laboratory water' : '#73e3c3', 'sea, bacterial mat': '#477038', 
                'pond, fresh water' : '#56b58d'}
    
    # Statistics for expeditions.
    Atlas_subset_main_sample['Sample_type_extended']=Atlas_subset_main_sample['Biome'] + ', ' + Atlas_subset_main_sample['Sample type']
    List_of_sample_types=sorted(list(set(Atlas_subset_main_sample['Sample_type_extended'].tolist())), reverse=True)
    print(f'Ordered list of sample types collected: {List_of_sample_types}')
    
    List_of_expeditions=list(set(Atlas_subset_main_sample['Expedition_ID'].tolist()))
    if np.nan in List_of_expeditions:
        List_of_expeditions.remove(np.nan)
    List_of_expeditions=sorted(List_of_expeditions)
    
    Expedition_dict={'Expedition_name' : [],
                     'Number_of_samples' : []}
    print('\nExpedition_name Number_of_samples')  
    for expedition_name in List_of_expeditions:
        number_of_samples=Atlas_subset_main_sample[Atlas_subset_main_sample['Expedition_ID']==expedition_name].shape[0]
        print(f'{expedition_name} {number_of_samples}')
        Expedition_dict['Expedition_name'].append(expedition_name)
        Expedition_dict['Number_of_samples'].append(number_of_samples)
        
    print(f'Number of expeditions organized: {len(List_of_expeditions)}')
    print(f'Number of samples {action}: {Atlas_subset_main_sample.shape[0]}\n')
    Expedition_dict['Expedition_name'].append('Number of expeditions organized:')
    Expedition_dict['Number_of_samples'].append(len(List_of_expeditions))  
    Expedition_dict['Expedition_name'].append(f'Number of samples {action}:')
    Expedition_dict['Number_of_samples'].append(Atlas_subset_main_sample.shape[0])  
    Expedition_df=pd.DataFrame.from_dict(Expedition_dict)
    write_sheet_to_excel(excel_path, Expedition_df, f'Expeditions_{action}_stat')
    
    if len(List_of_expeditions)>0:
        # Plot numbers of samples of different types with stacked barplot.
        expeditions_sample_types_barplot(Atlas_subset_main_sample, Color_dict, excel_path)
        
        # Prepare points for Qgis.
        prep_points_for_Qgis(Atlas_subset_main_sample, excel_path)
    
    # Statistics for locations.
    Locations_dict={'Locations' : [],
                    'Number_of_samples' : []}    
    List_of_coordinates=list(set(Atlas_subset_main_sample['Coordinates, decimal'].tolist()))
    if np.nan in List_of_coordinates:
        List_of_coordinates.remove(np.nan)
    List_of_coordinates=sorted(List_of_coordinates)

    print(f'Number of expeditions stations {action}: {len(List_of_coordinates)}\n')  
    Locations_dict['Locations'].append(f'Number of expeditions stations {action}:')
    Locations_dict['Number_of_samples'].append(len(List_of_coordinates)) 
    Location_df=pd.DataFrame.from_dict(Locations_dict)
    write_sheet_to_excel(excel_path, Location_df, f'Locations_{action}_stat')    
    
    # Statistics for sample types.     
    List_of_sample_types=list(set(Atlas_subset_main_sample['Sample type'].tolist()))
    if np.nan in List_of_sample_types:
        List_of_sample_types.remove(np.nan)
    List_of_sample_types=sorted(List_of_sample_types)
    
    Sample_type_dict={'Sample_type' : [],
                    'Number_of_samples' : []}     
    print('Sample_type Number_of_samples')     
    for sample_type in List_of_sample_types:
        number_of_samples=Atlas_subset_main_sample[Atlas_subset_main_sample['Sample type']==sample_type].shape[0]
        print(f'{sample_type} {number_of_samples}')
        Sample_type_dict['Sample_type'].append(sample_type)
        Sample_type_dict['Number_of_samples'].append(number_of_samples)        

    print(f'Number of sample types {action}: {len(List_of_sample_types)}\n')  
    Sample_type_dict['Sample_type'].append(f'Number of sample types {action}:')
    Sample_type_dict['Number_of_samples'].append(len(List_of_sample_types))  
    Sample_type_df=pd.DataFrame.from_dict(Sample_type_dict)
    write_sheet_to_excel(excel_path, Sample_type_df, f'Sample_types_{action}_stat')        
    
    # Statistics for seaweeds.
    Atlas_subset_main_sample_seaweed=Atlas_subset_main_sample[Atlas_subset_main_sample['Sample type']=='seaweed']
    List_of_seaweeds_main_sample=list(set(Atlas_subset_main_sample_seaweed['Host organism'].tolist()))
    if np.nan in List_of_seaweeds_main_sample:
        List_of_seaweeds_main_sample.remove(np.nan)
    List_of_seaweeds_main_sample=sorted(List_of_seaweeds_main_sample)
    
    Seaweed_dict={'Seaweed_name' : [],
                  'Number_of_samples' : []}       
    print('Seaweed_name Number_of_samples')  
    for Seaweed_name_main_sample in List_of_seaweeds_main_sample:
        number_of_samples=Atlas_subset_main_sample_seaweed[Atlas_subset_main_sample_seaweed['Host organism']==Seaweed_name_main_sample].shape[0]
        print(f'{Seaweed_name_main_sample} {number_of_samples}')
        Seaweed_dict['Seaweed_name'].append(Seaweed_name_main_sample)
        Seaweed_dict['Number_of_samples'].append(number_of_samples)            

    print(f'Number of seaweed species {action}: {len(List_of_seaweeds_main_sample)}') 
    print(f'Number of seaweed samples {action}: {Atlas_subset_main_sample_seaweed.shape[0]}\n')
    Seaweed_dict['Seaweed_name'].append(f'Number of seaweed species {action}:')
    Seaweed_dict['Number_of_samples'].append(len(List_of_seaweeds_main_sample)) 
    Seaweed_dict['Seaweed_name'].append(f'Number of seaweed samples {action}:')
    Seaweed_dict['Number_of_samples'].append(Atlas_subset_main_sample_seaweed.shape[0])      
    Seaweed_df=pd.DataFrame.from_dict(Seaweed_dict)
    write_sheet_to_excel(excel_path, Seaweed_df, f'Seaweed_{action}_stat')     
    
    # Statistics for animals, Main samples.
    Atlas_subset_main_sample_organism=Atlas_subset_main_sample[Atlas_subset_main_sample['Sample type']=='organism']
    List_of_organism_main_sample=list(set(Atlas_subset_main_sample_organism['Host organism'].tolist()))
    if np.nan in List_of_organism_main_sample:
        List_of_organism_main_sample.remove(np.nan)
    List_of_organism_main_sample=sorted(List_of_organism_main_sample)
    
    Animal_dict={'Animal_name' : [],
                 'Number_of_samples' : []}      
    print('Animal_name Number_of_samples')  
    for Organism_name_main_sample in List_of_organism_main_sample:
        number_of_samples=Atlas_subset_main_sample_organism[Atlas_subset_main_sample_organism['Host organism']==Organism_name_main_sample].shape[0]
        print(f'{Organism_name_main_sample} {number_of_samples}')
        Animal_dict['Animal_name'].append(Organism_name_main_sample)
        Animal_dict['Number_of_samples'].append(number_of_samples)          

    print(f'Number of animal species {action}: {len(List_of_organism_main_sample)}') 
    print(f'Number of animal samples {action}: {Atlas_subset_main_sample_organism.shape[0]}\n')  
    Animal_dict['Animal_name'].append(f'Number of animal species {action}:')
    Animal_dict['Number_of_samples'].append(len(List_of_organism_main_sample)) 
    Animal_dict['Animal_name'].append(f'Number of animal samples {action}:')
    Animal_dict['Number_of_samples'].append(Atlas_subset_main_sample_organism.shape[0])      
    Animal_df=pd.DataFrame.from_dict(Animal_dict)
    write_sheet_to_excel(excel_path, Animal_df, f'Animal_{action}_stat')       
    
    # Statistics for plants, Main samples.
    Atlas_subset_main_sample_plant=Atlas_subset_main_sample[Atlas_subset_main_sample['Sample type']=='plant']
    List_of_plant_main_sample=list(set(Atlas_subset_main_sample_plant['Host organism'].tolist()))
    if np.nan in List_of_plant_main_sample:
        List_of_plant_main_sample.remove(np.nan)
    List_of_plant_main_sample=sorted(List_of_plant_main_sample)
    
    Plant_dict={'Plant_name' : [],
                'Number_of_samples' : []}      
    print('Plant_name Number_of_samples')  
    for Plant_name_main_sample in List_of_plant_main_sample:
        number_of_samples=Atlas_subset_main_sample_plant[Atlas_subset_main_sample_plant['Host organism']==Plant_name_main_sample].shape[0]
        print(f'{Plant_name_main_sample} {number_of_samples}')
        Plant_dict['Plant_name'].append(Plant_name_main_sample)
        Plant_dict['Number_of_samples'].append(number_of_samples)  

    print(f'Number of plant species {action}: {len(List_of_plant_main_sample)}') 
    print(f'Number of plant samples {action}: {Atlas_subset_main_sample_plant.shape[0]}\n')  
    Plant_dict['Plant_name'].append(f'Number of plant species {action}:')
    Plant_dict['Number_of_samples'].append(len(List_of_plant_main_sample)) 
    Plant_dict['Plant_name'].append(f'Number of plant samples {action}:')
    Plant_dict['Number_of_samples'].append(Atlas_subset_main_sample_plant.shape[0])      
    Plant_df=pd.DataFrame.from_dict(Plant_dict)
    write_sheet_to_excel(excel_path, Plant_df, f'Plant_{action}_stat')      
    
    # Statistics for terrestrial soils, Main samples.
    Atlas_subset_main_sample_ter_soil=Atlas_subset_main_sample[(Atlas_subset_main_sample['Sample type']=='soil') & 
                                                               (Atlas_subset_main_sample['Expedition type']=='Terrestrial') & 
                                                               (Atlas_subset_main_sample['Location']!='Skoltech')]
    List_of_soil_locations_main_sample=list(set(Atlas_subset_main_sample_ter_soil['Location'].tolist()))
    if np.nan in List_of_soil_locations_main_sample:
        List_of_soil_locations_main_sample.remove(np.nan)
    List_of_soil_locations_main_sample=sorted(List_of_soil_locations_main_sample)
    
    Ter_soil_dict={'Reserve_name' : [],
                   'Number_of_samples' : []}      
    print('Reserve_name Number_of_samples')  
    for Soil_loc_name_main_sample in List_of_soil_locations_main_sample:
        number_of_samples=Atlas_subset_main_sample_ter_soil[Atlas_subset_main_sample_ter_soil['Location']==Soil_loc_name_main_sample].shape[0]
        print(f'{Soil_loc_name_main_sample} {number_of_samples}')
        Ter_soil_dict['Reserve_name'].append(Soil_loc_name_main_sample)
        Ter_soil_dict['Number_of_samples'].append(number_of_samples)         

    print(f'Number of natural reserves sampled: {len(List_of_soil_locations_main_sample)}') 
    print(f'Number of soil samples {action} in natural reserves: {Atlas_subset_main_sample_ter_soil.shape[0]}\n')  
    Ter_soil_dict['Reserve_name'].append(f'Number of natural reserves sampled:')
    Ter_soil_dict['Number_of_samples'].append(len(List_of_soil_locations_main_sample)) 
    Ter_soil_dict['Reserve_name'].append(f'Number of soil samples {action} in natural reserves:')
    Ter_soil_dict['Number_of_samples'].append(Atlas_subset_main_sample_ter_soil.shape[0])      
    Ter_soil_df=pd.DataFrame.from_dict(Ter_soil_dict)
    write_sheet_to_excel(excel_path, Ter_soil_df, f'Soil_reserve_{action}_stat')        
    
    
    return


########
# Read the whole Atlas metadata, calculate statistics.
########

def read_Atlas_db_return_statistics(atlas_db_inpath, output_path):
    
    Atlas_database=pd.read_excel(atlas_db_inpath, sheet_name="General_table", header=0)  #Ответы на форму (1)
    Atlas_master_table_name=atlas_db_inpath.split("\\")[-1][:-5]
    Excel_path=os.path.join(output_path, f'{Atlas_master_table_name}_stat.xlsx')   
    
    # Statistics for Main samples.
    Atlas_subset_main_sample=Atlas_database[Atlas_database['Processing step']=='Main sample']
    Action='collected'
    calc_write_statistics(Atlas_subset_main_sample, Action, Excel_path)
    
    # Statistics for sequenced samples with 16S data avaliable.
    Atlas_subset_16S=Atlas_database[Atlas_database['Processing step']=='DNA'] #[Atlas_database['Raw sequencing data'].str.len()>0]
    Action='16S_sequenced'
    calc_write_statistics(Atlas_subset_16S, Action, Excel_path)    

    
    
    
    
    #with pd.ExcelWriter(output_path) as writer:  
    #    Atlas_subset_16S.to_excel(writer, sheet_name='16S_samples')
    #    Atlas_subset_Shotgun.to_excel(writer, sheet_name='Shotgun_samples')    
    
    return

read_Atlas_db_return_statistics(Atlas_db_inpath, Output_path)