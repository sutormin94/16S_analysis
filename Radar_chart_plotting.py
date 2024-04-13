#Import the packages
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

## Disclaimer
## A code was partially taken from https://curbal.com/curbal-learning-portal/radar-charts-in-matplotlib


#Read file with ranks.
ranks_dataframe=pd.read_excel('C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\Radar_plots\\Final_ranks.xlsx', index_col=0, sheet_name='Ranks')

#########
## Plot 1. Total kit ranks.
#########

#Output file.
Outpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\Radar_plots\\Total_ranks_radar_plots"

#Parameters (column names) to be visualized.
cat_init=['DNA amount, total', 'DIN, total', 'Presence of inhibitors, total', '18S/16S, total', 'Contamination level, total', 'Reproducibility level, total', 'Alpha diversity, total']

######
#Rescale dataframe to make all values between 0 and 1.
######

def rescale_ranks(dataframe, column_name_list):
    
    for param in column_name_list:
        init_data=dataframe[param]
        rescale_data=(init_data-max(init_data))/(min(init_data)-max(init_data))
        dataframe[param]=rescale_data
    
    print(dataframe)    
    
    return dataframe

ranks_dataframe=rescale_ranks(ranks_dataframe, cat_init)

######
#Set a descending order of plots.
######

def desc_order_of_plots(dataframe, column_name_list):
    
    Orig_kits_order=list(dataframe.index)
    kit_areas_list=[]
    for kit in Orig_kits_order:
        kit_area=0
        for param in column_name_list:
            kit_area+=dataframe.loc[kit][param]
        kit_areas_list.append(kit_area)  
    
    kit_areas_list_sort, Kits_order_sort=zip(*sorted(zip(kit_areas_list, Orig_kits_order), reverse=True)) 
    
    return Kits_order_sort

Kits_order_sort=desc_order_of_plots(ranks_dataframe, cat_init)

######
#Prepare axis data.
######

def prep_axis(cat_init):
    
    cat=[*cat_init, cat_init[0]] #to close the radar, duplicate the first column
    label_loc=np.linspace(start=0, stop=2*np.pi, num=len(cat))    
    
    return cat, label_loc

cat, label_loc=prep_axis(cat_init)
cat_polished=[r'$\nu$(DNA)', 'DIN', 'Inh', '18S/16S', 'Cont', 'Rep', r'$\alpha$', r'$\nu$(DNA)']

#Plot data.
Number_of_plots=4 
Number_of_series=2

def plot_radar_set(Number_of_series, Number_of_plots, Kits_order_sort, cat, cat_polished, ranks_dataframe, label_loc, pseudocount, Set_title, Outpath):
    
    fig, plots=plt.subplots(Number_of_series, Number_of_plots, subplot_kw={'projection': 'polar'}, figsize=(Number_of_plots*3, Number_of_series*3.5), dpi=100)
    rgb = plt.colormaps['rainbow']
    
    i=0
    for kit in Kits_order_sort:
        
        cust_color=rgb(((Number_of_plots*Number_of_series)-i)/(Number_of_plots*Number_of_series))
        
        kit_data=[]
        for param in cat:
            kit_point=ranks_dataframe.loc[kit][param]
            kit_data.append(kit_point+pseudocount)
        
        x=i//Number_of_plots
        y=i%Number_of_plots
        plots[x, y].plot(label_loc, kit_data, 'o-', linewidth=2, color=cust_color)
        plots[x, y].fill(label_loc, kit_data, alpha=0.45, color=cust_color)
        plots[x, y].set_theta_offset(np.pi/2)
        plots[x, y].set_theta_direction(-1)   
        plots[x, y].set_thetagrids(np.degrees(label_loc), cat_polished, size=13)
        plots[x, y].set_ylim(0, 1.1+pseudocount)
        plots[x, y].tick_params(axis='y', labelsize=0)
        plots[x, y].set_yticks([pseudocount, pseudocount+0.5, pseudocount+1], minor=False)
        plots[x, y].grid(color='#7f929c', linestyle='--', linewidth=1.5, which='major')
        plots[x, y].spines['polar'].set_color('white')
        plots[x, y].set_facecolor('white')
        plots[x, y].set_title(f'{kit}', y=1.08, size=14, fontweight="bold", pad=20)
        
        i+=1
    
    fig.suptitle(f'{Set_title}', size=16, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f'{Outpath}.png', dpi=300)
    plt.savefig(f'{Outpath}.svg', dpi=300)
    plt.show() 
    plt.close()
    
    return

Pseudocount=0.1
Set_title_total=""
plot_radar_set(Number_of_series, Number_of_plots, Kits_order_sort, cat, cat_polished, ranks_dataframe, label_loc, Pseudocount, Set_title_total, Outpath)



#########
## Plot 2. Kit ranks for soil, water, and guts samples.
#########

#Output file.
Outpath_total="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\Radar_plots\\Sample_specific_ranks_radar_plots"

#Rescale dataframe to make all values between 0 and 1.
cat_init_soil=['DNA amount, soil', 'DIN, soil', 'Presence of inhibitors, soil', '18S/16S, soil', 'Contamination level, soil', 'Reproducibility level, soil', 'Alpha diversity, soil']
ranks_dataframe=rescale_ranks(ranks_dataframe, cat_init_soil)
cat_init_water=['DNA amount, water', 'DIN, water', 'Presence of inhibitors, water', '18S/16S, water', 'Contamination level, water', 'Reproducibility level, water', 'Alpha diversity, water']
ranks_dataframe=rescale_ranks(ranks_dataframe, cat_init_water)
cat_init_guts=['DNA amount, guts', 'DIN, guts', 'Presence of inhibitors, guts', '18S/16S, guts', 'Contamination level, guts', 'Reproducibility level, guts', 'Alpha diversity, guts']
ranks_dataframe=rescale_ranks(ranks_dataframe, cat_init_guts)

#Set a descending order of plots.
Kits_order_sort_soil =desc_order_of_plots(ranks_dataframe, cat_init_soil)
Kits_order_sort_water=desc_order_of_plots(ranks_dataframe, cat_init_water)
Kits_order_sort_guts =desc_order_of_plots(ranks_dataframe, cat_init_guts)

#Prepare axis data.
cat_soil,  label_loc_soil =prep_axis(cat_init_soil)
cat_water, label_loc_water=prep_axis(cat_init_water)
cat_guts,  label_loc_guts =prep_axis(cat_init_guts)

#Plot data all together.
Number_of_plots_all=3
Number_of_series_all=len(Kits_order_sort_soil)

fig, plots=plt.subplots(Number_of_series_all, Number_of_plots_all, subplot_kw={'projection': 'polar'}, figsize=(Number_of_plots_all*3, Number_of_series_all*3.5), dpi=100)
rgb=plt.colormaps['rainbow']


def plot_all_radars(ranks_dataframe, Kits_order_sort, cat, cat_polished, label_loc, plots, y, sample_type, Number_of_series, Pseudocount):

    i=0
    for kit in Kits_order_sort:
        
        cust_color=rgb(((Number_of_series)-i)/(Number_of_series))
        
        kit_data=[]
        for param in cat:
            kit_point=ranks_dataframe.loc[kit][param]
            kit_data.append(kit_point+Pseudocount)
        
        x=i
        plots[x, y].plot(label_loc, kit_data, 'o-', linewidth=2, color=cust_color)
        plots[x, y].fill(label_loc, kit_data, alpha=0.45, color=cust_color)
        plots[x, y].set_theta_offset(np.pi/2)
        plots[x, y].set_theta_direction(-1)   
        plots[x, y].set_thetagrids(np.degrees(label_loc), cat_polished, size=13)
        plots[x, y].set_ylim(0, 1.1+Pseudocount)
        plots[x, y].tick_params(axis='y', labelsize=0)
        plots[x, y].grid(color='#7f929c', linestyle='--', linewidth=1.5)
        plots[x, y].spines['polar'].set_color('white')
        plots[x, y].set_facecolor('white')
        if x==0:
            plots[x, y].set_title(f'{sample_type}\n{kit}', y=1.08, size=14, fontweight="bold", pad=20)  
        else:
            plots[x, y].set_title(f'{kit}', y=1.08, size=14, fontweight="bold", pad=20)
        
        i+=1
    
    return

Pseudocount=0.1
y_pos=0
sample_type='Soil'
plot_all_radars(ranks_dataframe, Kits_order_sort_soil, cat_soil, cat_polished, label_loc_soil, plots, y_pos, sample_type, Number_of_series_all, Pseudocount)
y_pos=1
sample_type='Water'
plot_all_radars(ranks_dataframe, Kits_order_sort_water, cat_water, cat_polished, label_loc_water, plots, y_pos, sample_type, Number_of_series_all, Pseudocount)
y_pos=2
sample_type='Guts'
plot_all_radars(ranks_dataframe, Kits_order_sort_guts, cat_guts, cat_polished, label_loc_guts, plots, y_pos, sample_type, Number_of_series_all, Pseudocount)

plt.tight_layout()
plt.savefig(f'{Outpath_total}.png', dpi=300)
plt.savefig(f'{Outpath_total}.svg', dpi=300)
plt.show()
plt.close()

#Plot separate plots for diffrent sample types.
Outpath_soil="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\Radar_plots\\Soil_specific_ranks_radar_plots"
Set_title_soil="Kits' ranks for soil samples"
plot_radar_set(Number_of_series, Number_of_plots, Kits_order_sort_soil, cat_soil, cat_polished, ranks_dataframe, label_loc_soil, Pseudocount, Set_title_soil, Outpath_soil)
Outpath_water="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\Radar_plots\\Water_specific_ranks_radar_plots"
Set_title_water="Kits' ranks for water samples"
plot_radar_set(Number_of_series, Number_of_plots, Kits_order_sort_water, cat_water, cat_polished, ranks_dataframe, label_loc_water, Pseudocount, Set_title_water, Outpath_water)
Outpath_guts="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Bioinformatics\Kit_comparison_analysis\Radar_plots\\Guts_specific_ranks_radar_plots"
Set_title_guts="Kits' ranks for guts samples"
plot_radar_set(Number_of_series, Number_of_plots, Kits_order_sort_guts, cat_guts, cat_polished, ranks_dataframe, label_loc_guts, Pseudocount, Set_title_guts, Outpath_guts)