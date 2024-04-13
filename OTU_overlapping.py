###############################################
##Dmitry Sutormin, 2023##
##DNA extraction kits comparison##

#Takes table with OTU or ASV data after decontamination and calculates overlapping of sets between different samples groups.
###############################################

#######
#Packages to be imported.
#######

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib_venn import venn2, venn3, venn3_circles, venn2_circles
import seaborn as sns
import matplotlib.pyplot as plt


#Path to OTU/ASV table.
#Either paired 0.5.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTU_frequency_decontam_either_05.tsv"
#Either paired 0.1.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTU_frequency_decontam_either_01_grouped.tsv"
#Prevalence paired 0.5.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTU_frequency_decontam_prev_05.tsv"
#Either all vs all 0.5.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTU_frequency_decontam_either_05_all_vs_all.tsv"
#Prevalence all vs all 0.5.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTU_frequency_decontam_prev_05_all_vs_all.tsv"
#Either all vs all 0.1.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTU_frequency_decontam_either_01_all_vs_all.tsv"
#Prevalence all vs all 0.1.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTU_frequency_decontam_prev_01_all_vs_all.tsv"
#No filtration.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\\all_OTU_frequency_no_controls.tsv"
#Prevalence groups 0.1.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTU_frequency_decontam_prevalence_0.1_by_groups.tsv"
#Prevalence groups 0.5.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTU_frequency_decontam_prevalence_0.5_by_groups.tsv"
#Either groups 0.1.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTU_frequency_decontam_either_0.1_by_groups.tsv"
#Either groups 0.5.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTU_frequency_decontam_either_0.2_by_groups.tsv"

#Either grouped by sample type all vs all 0.5.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTUs_decontam_either_0.5_grouped_all_vs_all.tsv"
#Either grouped by sample type all vs all 0.3.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTUs_decontam_either_0.3_grouped_all_vs_all.tsv"
#Prevalence grouped by sample type all vs all 0.5.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTUs_decontam_prevalence_0.5_grouped_all_vs_all.tsv"

#Either grouped by control all vs all 0.5.
Dataset_name='OTUs_decontam_either_0.5_grouped_by_control_all_vs_all'
Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTUs_decontam_either_0.5_grouped_by_control_all_vs_all.tsv"

#Either grouped by control all vs all 0.3.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTUs_decontam_either_0.3_grouped_by_control_all_vs_all.tsv"
#Either grouped by control all vs all 0.1.
#Table_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\OTUs_decontam_either_0.1_grouped_by_control_all_vs_all.tsv"



#Path to metadata.
Table_metadata_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\Metadata_conc.tsv"

#Output path.
Pathout="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\OTU_replicating\\"

#######
## Prepare OTU/ASV data and metadata.
#######

#Read input metadata table.
Input_metadata=pd.read_csv(Table_metadata_path, sep='\t', header=0, index_col=0)
Sample_type_ar=list(Input_metadata.loc[:, "Sample_type"])
Kit_type_list=list(Input_metadata.loc[:, "Kits_type"])
Samples_ID_ar=list(Input_metadata.index)

Sample_kit_types=[]
Sample_kit_types_dict={}
for i in range(len(Sample_type_ar)):
    Sample_kit_type=Sample_type_ar[i]+'_'+Kit_type_list[i]+'_'+str((i%3)+1)
    Sample_kit_types.append(Sample_kit_type)
    Sample_kit_types_dict[Samples_ID_ar[i]]=Sample_kit_type
print(Sample_kit_types_dict)

#Read input data table.
Input_data=pd.read_csv(Table_data_path, sep='\t', header=0, index_col=0)
Input_data_col_names_old=list(Input_data.columns)

#Rename samples.
Input_data_col_names_new=[]
for col_name_old in Input_data_col_names_old:
    col_name_new=Sample_kit_types_dict[int(col_name_old)]
    Input_data_col_names_new.append(col_name_new)

print(Input_data_col_names_new)

Input_data_nh=Input_data.set_axis(Input_data_col_names_new, axis=1, inplace=False)
print(Input_data_nh)

#######
## Average tech replicates.
#######

## Method 1. OTU/ASV should be in all 3 tech replicates to be included - Absolute replication.
## Method 2. OTU/ASV should be in either of 3 tech replicates to be included - Relaxed replication.

#Prepare kit names lists.
Kit_types_list_short=[]
Sample_types_list_short=[]
Sample_kit_combinations=[]
for sample_name in Input_data_col_names_new:
    sample_name=sample_name[:-2]
    print(sample_name)
    kit_name=sample_name.split('_')[1]
    if (kit_name not in Kit_types_list_short) and (kit_name not in ['PowerFilter']):
        Kit_types_list_short.append(kit_name)
    sample_type=sample_name.split('_')[0]
    if sample_type not in Sample_types_list_short:
        Sample_types_list_short.append(sample_type)
    if sample_name in Sample_kit_combinations:
        continue
    else:
        Sample_kit_combinations.append(sample_name)


Number_of_plots=8  
Number_of_series=3
fig, plots=plt.subplots(Number_of_series, Number_of_plots, figsize=(Number_of_plots*2.5, Number_of_series*2.5), dpi=100)
        
Replicated_OTUs_dict_abs={}
Replicated_OTUs_dict_rel={}

Abs_replication_df=pd.DataFrame(columns=Sample_types_list_short, index=Kit_types_list_short)
Abs_replication_df=Abs_replication_df.astype(float)
print(Abs_replication_df)

j=0
for sk_comb in Sample_kit_combinations:
    sample_type=sk_comb.split('_')[0]
    kit_type=sk_comb.split('_')[1]
    if kit_type=='PowerFilter':
        kit_type='PowerSoil'
    Tech_rep_ar=[]
    for i in range(3):
        Tech_rep_ar.append(sk_comb+'_'+str(i+1))
    Tech_set_df=Input_data_nh.loc[:, Tech_rep_ar]
    Tech_set_df_filt_abs=Tech_set_df[(Tech_set_df[Tech_rep_ar[0]]>0) & (Tech_set_df[Tech_rep_ar[1]]>0) & (Tech_set_df[Tech_rep_ar[2]]>0)] #Absolute replication.
    Tech_set_df_filt_rel=Tech_set_df[(Tech_set_df[Tech_rep_ar[0]]>0) | (Tech_set_df[Tech_rep_ar[1]]>0) | (Tech_set_df[Tech_rep_ar[2]]>0)] #Relaxed replication.
    Abs_rep_list=list(Tech_set_df_filt_abs.index)
    Rel_rep_list=list(Tech_set_df_filt_rel.index)
    Replicated_OTUs_dict_abs[sk_comb]=Abs_rep_list
    Replicated_OTUs_dict_rel[sk_comb]=Rel_rep_list
    
    print(f'List of triplicated OTUs for sample {sk_comb}: {list(Tech_set_df_filt_abs.index)}')
    
    Rep_1=set(list(Tech_set_df[Tech_set_df[Tech_rep_ar[0]]>0].index))
    Rep_2=set(list(Tech_set_df[Tech_set_df[Tech_rep_ar[1]]>0].index))
    Rep_3=set(list(Tech_set_df[Tech_set_df[Tech_rep_ar[2]]>0].index))
    
    x_coord=j%8
    y_coord=j//8
    
    print(x_coord, y_coord)
    if y_coord==0:
        plots[y_coord, x_coord].set_title(sk_comb.split('_')[1], loc='center', y=1.1, fontsize=16)
        
    if x_coord==0:
        plots[y_coord, x_coord].text(-1.2, 0.3, sk_comb.split('_')[0], verticalalignment='center', rotation=90, size=16)
        
    venn3([Rep_1, Rep_2, Rep_3], set_labels = ('Rep 1', 'Rep 2', 'Rep 3'), ax=plots[y_coord, x_coord])
    venn3_circles([Rep_1, Rep_2, Rep_3], ax=plots[y_coord, x_coord])
    Abs_replication_level=int(100*float(len(Abs_rep_list))/len(Rel_rep_list))
    plots[y_coord, x_coord].text(0.4, -0.8, f'{len(Abs_rep_list)}/{len(Rel_rep_list)} OTUs\ntriplicated ({Abs_replication_level}%)', size=8)
    
    #Store replication level.
    Abs_replication_df.loc[kit_type, sample_type]=Abs_replication_level
    j+=1
    
#plt.tight_layout()
plt.savefig(f'{Pathout}Tech_replicates_venn.png', dpi=320)
plt.savefig(f'{Pathout}Tech_replicates_venn.svg', dpi=320)
plt.close()

print(Abs_replication_df)
Min=Abs_replication_df.min().min()
Max=Abs_replication_df.max().max()

fig, ax=plt.subplots(figsize = (3, 4))

#sns.heatmap(Control_df, annot=True)
ax.set_title('Reproducibility rate')
sns.heatmap(Abs_replication_df, annot=True, square=True, vmin=Min, vmax=Max, cmap='Blues')

plt.tight_layout()
plt.savefig(f'{Pathout}Absolute_reproducibility_heatmap.png', dpi=300)
plt.savefig(f'{Pathout}Absolute_reproducibility_heatmap.svg', dpi=300)


    
## Investigate reproducibility between kits using Absolute replication.

OTU_replication_numbers_by_st={}
sample_type_list=['soil', 'water', 'organism']
for sample_type in sample_type_list:
    OTU_replication_numbers={}
    for sk_comb, OTU_list in Replicated_OTUs_dict_abs.items():
        if sample_type in sk_comb:
            for OTU in OTU_list:
                if OTU in OTU_replication_numbers:
                    OTU_replication_numbers[OTU]+=1
                else:
                    OTU_replication_numbers[OTU]=1
    OTU_replication_numbers_by_st[sample_type]=OTU_replication_numbers



Rep_num_by_sample_type={}
Universal_OTU_num_dict={}
for sample_type, st_data in OTU_replication_numbers_by_st.items():
    Universal_OTU_num_dict[sample_type]=0
    Rep_num_ar=[]
    for OTU_name, OTU_rep_num in st_data.items():
        Rep_num_ar.append(OTU_rep_num)
        if OTU_rep_num==len(Kit_types_list_short):
            Universal_OTU_num_dict[sample_type]+=1
    Rep_num_by_sample_type[sample_type]=Rep_num_ar
    
    fig=plt.figure(figsize=(3, 3), dpi=100)
    plt.hist(Rep_num_ar, bins=range(1, 8 + 2, 1), facecolor='blue', edgecolor='black', alpha=0.5)
    plt.xticks(np.array(range(1, 8 + 1, 1))+0.5, range(1, 8 + 1, 1))
    plt.xlim([0,10])
    plt.xlabel('Number of kits')
    plt.ylabel('Number of OTUs')
    plt.title(f'{sample_type} samples')
    plt.tight_layout()
    plt.savefig(f'{Pathout}Absolute_replication_for_{sample_type}.png', dpi=300, figsize=(4, 4))
    plt.savefig(f'{Pathout}Absolute_replication_for_{sample_type}.svg', dpi=300, figsize=(4, 4))
    plt.close()           
    
#print(Rep_num_by_sample_type)
#print(Kit_types_list_short)
print(Universal_OTU_num_dict)


#######
#Make barplot of the number of OTUs of different types: universal, other, unique.
#######

def barplot_OTUs_types(otu_num_dict, Types_name_ar, Color_ar, kit_types_list, Sample_type, Dataset_name, path_out):
    
    fig, ax = plt.subplots()
    
    j=0
    x_coords=[]
    for kit_name, OTU_num_ar in otu_num_dict.items():
    
        Num_cumul=0
        for i in range(len(OTU_num_ar)):
            if j==0:
                ax.bar(j, OTU_num_ar[i], 0.65, color=Color_ar[i], bottom=Num_cumul, tick_label=kit_name, label=f'{Types_name_ar[i]}')
            else:
                ax.bar(j, OTU_num_ar[i], 0.65, color=Color_ar[i], bottom=Num_cumul, tick_label=kit_name)
            Num_cumul+=OTU_num_ar[i]
        
        x_coords.append(j)
        j+=1
        
    
    #ax.set_xlim([-0.5,0.5])
    ax.tick_params(axis='x', which='major', labelsize=12)
    ax.set_xticks(x_coords, labels=kit_types_list, minor=False, size=12, rotation=90)
    ax.set_ylabel('Number of OTUs', size=12)
    #ax.set_xlabel('DNA-extraction kits', size=12)
    ax.set_title(Sample_type, size=15)
    ax.spines['top'].set_color(None)
    ax.spines['right'].set_color(None)
    
    #Place legend outside of a graph. Taken from: https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
    box=ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(fontsize=12, ncol=1, loc='upper left', frameon=False, markerscale=1, handlelength=0.7, handletextpad=0.3, columnspacing=0.7, bbox_to_anchor=(1, 0.95))
    plt.tight_layout(rect=[0,0,0.7,1])
    plt.show()
    plt.savefig(f'{path_out}{Sample_type}_{Dataset_name}_OTUs_types_barplot.png', dpi=300, figsize=(2, 5))
    plt.savefig(f'{path_out}{Sample_type}_{Dataset_name}_OTUs_types_barplot.svg', dpi=300, figsize=(2, 5))    
    
    return

## Get unique OTUs, provided by just one kit.

sample_type_list=['soil', 'water', 'organism']
Color_ar=['#de7b42', '#5dc1e3', '#e880aa']
OTU_types_list=['Universal OTUs', 'Other OTUs', 'Unique OTUs']

for sample_type in sample_type_list:
    OTU_num_dict={}
    for kit_name_sel in Kit_types_list_short:
        Selected_comb=f'{sample_type}_{kit_name_sel}'
        Other_comb=[]
        All_OTUs=[]
        for kit_name_oth in Kit_types_list_short:
            All_OTUs+=Replicated_OTUs_dict_abs[f'{sample_type}_{kit_name_oth}']
            if kit_name_sel!=kit_name_oth:
                Other_comb.append(f'{sample_type}_{kit_name_oth}')
        
        Other_kits_OTU=[]
        for other_kit_name in Other_comb:
            Other_kits_OTU+=Replicated_OTUs_dict_abs[other_kit_name]
        Other_kits_OTU_nr=list(set(Other_kits_OTU))
        Selected_kit_OTU_nr=list(set(Replicated_OTUs_dict_abs[Selected_comb]))
        
        Unique_OTU=[]
        for sel_kit_OTU in Selected_kit_OTU_nr:
            if sel_kit_OTU not in Other_kits_OTU_nr:
                Unique_OTU.append(sel_kit_OTU)
                
        print(f'{sample_type} {kit_name_sel}')
        print("Unique OTUs num/Selected kit's OTU num/Other kits OTU/Total num of OTU")
        print(f'{len(Unique_OTU)}/{len(Selected_kit_OTU_nr)}/{len(Other_kits_OTU_nr)}/{len(list(set(All_OTUs)))}')
        
        Number_of_universal_OTUs=Universal_OTU_num_dict[sample_type] 
        Number_of_unique_OTUs=len(Unique_OTU)
        Number_of_other_OTUs=len(Selected_kit_OTU_nr)-Number_of_universal_OTUs-Number_of_unique_OTUs
        
        OTU_num_dict[kit_name_sel]=[Number_of_universal_OTUs, Number_of_other_OTUs, Number_of_unique_OTUs]
    
    barplot_OTUs_types(OTU_num_dict, OTU_types_list, Color_ar, Kit_types_list_short, sample_type, Dataset_name, Pathout)


        

## Investigate reproducibility between kits using Relaxed replication.

OTU_replication_numbers_by_st={}
sample_type_list=['soil', 'water', 'organism']
for sample_type in sample_type_list:
    OTU_replication_numbers={}
    for sk_comb, OTU_list in Replicated_OTUs_dict_rel.items():
        if sample_type in sk_comb:
            for OTU in OTU_list:
                if OTU in OTU_replication_numbers:
                    OTU_replication_numbers[OTU]+=1
                else:
                    OTU_replication_numbers[OTU]=1
    OTU_replication_numbers_by_st[sample_type]=OTU_replication_numbers


Rep_num_by_sample_type={}
for sample_type, st_data in OTU_replication_numbers_by_st.items():
    Rep_num_ar=[]
    for OTU_name, OTU_rep_num in st_data.items():
        Rep_num_ar.append(OTU_rep_num)
    Rep_num_by_sample_type[sample_type]=Rep_num_ar
    
    fig=plt.figure(figsize=(3, 3), dpi=100)
    plt.hist(Rep_num_ar, bins=range(1, 8 + 2, 1), facecolor='blue', edgecolor='black', alpha=0.5)
    plt.xticks(np.array(range(1, 8 + 1, 1))+0.5, range(1, 8 + 1, 1))
    plt.xlim([0,10])
    plt.xlabel('Number of kits')
    plt.ylabel('Number of OTUs')
    plt.title(f'{sample_type} samples')
    plt.tight_layout()
    plt.savefig(f'{Pathout}Relaxed_replication_for_{sample_type}.png', dpi=300, figsize=(4, 4))
    plt.close()           
    
#print(Rep_num_by_sample_type)

        


