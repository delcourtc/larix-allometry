# -*- coding: utf-8 -*-
"""
@author: C. J. F. Delcourt
"""

#%% importing required packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# setting the working directory in which the excel files from
# Delcourt and Veraverbeke (2022) were saved
wdir = ""

# loading size-class-specific G and MSD values created from larix_fwd.R script
df_fwd = pd.read_csv(wdir+"outputs/table3_fwd_summary.csv")

#%% values of the multiplication factor M

# storing size-class-specific G and MSD values from from our study and those
# from other tree species in the Canadian Northwest Territories and
# Saskatchewan by diameter size class
boreal_dict = {'lcajYA': {'G': df_fwd.iloc[1:,3].to_list(),
                          'MSD': [round(x,2) for x in df_fwd.iloc[1:,1].to_list()]},
               
               'llarSK' : {'G' : [0.51, 0.51, 0.49, 0.55],
                           'MSD' : [0.475, 2.60, 14.1, 40.6]},
       
               'pmarSK' : {'G' : [0.56, 0.51, 0.49, 0.49],
                           'MSD' : [0.487, 3.49, 15.7, 33.8]},
       
                'pglaSK' : {'G' : [0.54, 0.54, 0.46, 0.41],
                            'MSD' : [0.528, 3.28, 15.8, 34.2]},
       
                'pmarNT' : {'G' : [0.62, 0.59, 0.55, 0.52],
                            'MSD' : [0.491, 3.573, 15.0, 34.7]},
       
                'pglaNT' : {'G' : [0.56, 0.54, 0.49, 0.45],
                            'MSD' : [0.498, 3.248, 15.5, 36.5]}}

# computing M values using Eq. (8) from our paper
for keys in boreal_dict.keys():
    boreal_dict[keys]['M'] = []
    for i in range(4):
        temp = (boreal_dict[keys]['G'][i]*1.13*boreal_dict[keys]['MSD'][i]*np.pi**2)/8
        boreal_dict[keys]['M'].append(temp)
del temp, i, keys

# computing differences between values of M from this study and those from
# other tree species and boreal regions for each size class
boreal_dict.keys()
diff_M = {}
for keys in list(boreal_dict.keys())[1:]:
    diff_M[keys] = ((np.array(boreal_dict['lcajYA']['M'])-np.array(boreal_dict[keys]['M']))/np.array(boreal_dict['lcajYA']['M']))*100    
del keys
    
print(np.mean(np.concatenate((diff_M['llarSK'], diff_M['pmarSK'],
                              diff_M['pglaSK'], diff_M['pmarNT'],
                              diff_M['pglaNT']))))

print(np.mean(np.abs(np.concatenate((diff_M['llarSK'], diff_M['pmarSK'],
                                     diff_M['pglaSK'], diff_M['pmarNT'],
                                     diff_M['pglaNT'])))))
 
print(np.max(np.concatenate((diff_M['llarSK'], diff_M['pmarSK'],
                             diff_M['pglaSK'], diff_M['pmarNT'],
                             diff_M['pglaNT']))))

print(np.min(np.concatenate((diff_M['llarSK'], diff_M['pmarSK'],
                             diff_M['pglaSK'], diff_M['pmarNT'],
                             diff_M['pglaNT']))))

#%% FWD biomass estimates per size class in 47 larch forest stands
    
# loading plot characteristics
df_plot = pd.read_excel(wdir+'YA2019_plots.xlsx', sheet_name='plots_summary')
# calculating slope correction factor (s) using Eq. (2) from our paper
df_plot['slope_corr'] = np.sqrt(1+(np.tan(np.radians(df_plot['slope'])))**2)
# storing s values per plot in a dictionary
slope_dict = dict(zip(df_plot.plotID, df_plot.slope_corr))

# loading FWD inventory data collected using the line-intersect method
df_count = pd.read_excel(wdir+"YA2019_fwd_transects.xlsx",
                         sheet_name="fwd_count")
# working only with larch FWD and pieces larger than 0.5 cm in diameter
df_count = df_count[(df_count.species=='LC') & (df_count.size_class!=1)]

# creating a function to calculate FWD biomass per size class and plot
# using Eq (4) from our paper
def fwd_pre (df,ref,tilt_corr):
    return ((np.pi**2)*df['count']*boreal_dict[ref]['G'][df['size_class']-2]*\
             boreal_dict[ref]['MSD'][df['size_class']-2]*tilt_corr*\
                 slope_dict[df['plotID']])/(8*30)

# deriving FWD biomass estimates using values of G and MSD from our study and
# those from other tree species and boreal regions
for ref in boreal_dict.keys():
    df_count['prefwd_'+ref] = df_count.apply(lambda row: fwd_pre(row,ref,1.13),
                                             axis=1)
    
#%% creating Table A1

table_a1 = df_count.set_index('size_class')
table_a1 = [pd.DataFrame(y).reindex([2,3,4,5]) for x, y in table_a1.groupby('plotID',
                                                                            as_index=False)]

for i in range(len(table_a1)):
    table_a1[i].plotID.fillna(method='bfill', inplace=True)
    table_a1[i].plotID.fillna(method='pad', inplace=True)
    table_a1[i].species.fillna(method='bfill', inplace=True)
    table_a1[i].species.fillna(method='pad', inplace=True)
    table_a1[i].iloc[:,2:] = table_a1[i].iloc[:,2:].fillna(0)
del i
    
table_a1 = pd.concat(table_a1)
table_a1.reset_index(drop=False,inplace=True)

class_count = table_a1[table_a1['count']!=0]
print(class_count.groupby('size_class')['plotID'].count())

table_a1 = table_a1.groupby('size_class').agg({'count':['mean', 'min', 'max'],
                                               'prefwd_llarSK': ['mean', 'std'],
                                               'prefwd_pmarSK': ['mean', 'std'],
                                               'prefwd_pglaSK': ['mean', 'std'],
                                               'prefwd_pmarNT': ['mean', 'std'],
                                               'prefwd_pglaNT': ['mean', 'std'],
                                               'prefwd_lcajYA': ['mean', 'std']})

df_all = df_count.groupby('plotID', as_index=False).agg({'count':sum,
                                                         'prefwd_llarSK':sum,
                                                         'prefwd_pmarSK':sum,
                                                         'prefwd_pglaSK':sum,
                                                         'prefwd_pmarNT':sum,
                                                         'prefwd_pglaNT':sum,
                                                         'prefwd_lcajYA':sum})

all_classes = df_all.agg({'count':['mean', 'min', 'max'],
                          'prefwd_llarSK': ['mean', 'std'],
                          'prefwd_pmarSK': ['mean', 'std'],
                          'prefwd_pglaSK': ['mean', 'std'],
                          'prefwd_pmarNT': ['mean', 'std'],
                          'prefwd_pglaNT': ['mean', 'std'],
                          'prefwd_lcajYA': ['mean', 'std']})
    
#%% differences in FWD biomass estimates

# calculating percentage difference in FWD biomass estimates in the 47 larch
# forest stands near Yakutsk using M factors derived for other species and
# boreal regions.
for ref in list(boreal_dict.keys())[1:]:
    df_all['diff_'+ref] = ((df_all['prefwd_lcajYA']-df_all['prefwd_'+ref])/df_all['prefwd_lcajYA'])*100 

data = [df_all.iloc[:,i].to_list() for i in np.arange(df_all.shape[1]-5,
                                                      df_all.shape[1])]
data = [data[x] for x in [0,2,1,4,3]]

#%% plotting Figure 3

plt.rcParams["font.family"] = "Arial"
plt.rcParams['figure.dpi'] = 300
cm = 1/2.54

spplabels = ['$\it{L. laricina}$, SK', '$\it{P. glauca}$, SK',
             '$\it{P. mariana}$, SK', '$\it{P. glauca}$, NT',
             '$\it{P. mariana}$, NT']
sppcolors = ['#285185', '#3669AC', '#6081D0', '#979AE6', '#B9B3F0']
xtickpos = [1,2,3,4,5]

fig1, ax1 = plt.subplots(figsize=(9*cm,9.85*cm))

bp = ax1.boxplot(data, patch_artist=True, showmeans=True)

# changing color of boxes
for patch, color in zip(bp['boxes'], sppcolors):
    patch.set_facecolor(color)
#changing linewidth of boxes
for box in bp['boxes']:
    box.set(linewidth=0.8)
#changing linewidth of whiskers
for whsk in bp['whiskers']:
    whsk.set(linewidth=0.8)
# changing color of medians
for median in bp['medians']:
    median.set(color ='w', linewidth=0.8)
# changing style and color of means
for mean in bp['means']:
    mean.set(marker="*",
             markersize=5,
             markerfacecolor = "white",
             markeredgecolor = "white")
# changing style and color of fliers
for flier, color in zip(bp['fliers'], sppcolors):
    flier.set_markeredgecolor(color)
    flier.set(marker='+',
              markersize=5,
              markeredgewidth=0.5)

ax1.set_ylabel("Percentage difference (%)", fontsize=7)
ax1.set_xticks(xtickpos)
ax1.set_xticklabels(spplabels, rotation = 45, fontsize=7)
ax1.tick_params(axis='y', which='major', direction='in',
                length=3, right=True, labelsize=7)
ax1.tick_params(axis='x', which='both', bottom=False, top=False)

fig1.tight_layout()
fig1.savefig(wdir+'outputs/Figures/figure3.png')
