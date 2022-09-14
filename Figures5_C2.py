# -*- coding: utf-8 -*-
"""
@author: C. J. F. Delcourt
"""

#%% importing required packages
import os
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# setting the working directory in wich the codes were saved
os.chdir('')
# importing dictionaries containing allometric coefficients from our nonlinear
# regressions and from existing allometric models
from Figures4_C1 import wls_dict, lrx_allo

# setting the working directory in which the excel files from
# Delcourt and Veraverbeke (2022) were saved
wdir = ''

#%% data preparation for plotting Figure 5 and Figure C2

# Larch aboveground biomass predicted for sites in Yakutia
# using different allometry models 

# adding coefficients from Magadan and Ust-Yansky sites in the lrx_allo dictionary
lrx_allo['MG'] = {'a': [wls_dict['MG']['a'][i] for i in [0,3,4,5]],
                  'b': [wls_dict['MG']['b'][i] for i in [0,3,4,5]]}
lrx_allo['UY'] = {'a': [wls_dict['UY']['a'][i] for i in [0,3,4,5]],
                  'b': [wls_dict['UY']['b'][i] for i in [0,3,4,5]]}

# loading field measurements collected within tree transects
df_trees = pd.read_excel(wdir+"YA2019_tree_transects.xlsx",
                         sheet_name="dbh")
# working only with larch trees
df_trees = df_trees[df_trees['species']=='LC']

# creating a dictionary in which we will store stem, branches, foliage, and
# total aboveground biomass values for each forest stand calculated using
# the different allometry models
dict_AGB = {}
bcomp = ["stem", "branches", "foliage", "above"]

for keys in lrx_allo.keys():
    
    tp = copy.deepcopy(df_trees)
    
    for i in range(4):
        
        if (keys in ["kajimoto_CK", "kajimoto_OM"]) & (i==3):  
            tp["W_"+bcomp[i]] = \
                lrx_allo[keys]["a"][i-3]*(tp["dbh"]**lrx_allo[keys]["b"][i-3]) +\
                    lrx_allo[keys]["a"][i-2]*(tp["dbh"]**lrx_allo[keys]["b"][i-2]) +\
                        lrx_allo[keys]["a"][i-1]*(tp["dbh"]**lrx_allo[keys]["b"][i-1])
        else:
            tp["W_"+bcomp[i]] = \
                lrx_allo[keys]["a"][i]*(tp["dbh"]**lrx_allo[keys]["b"][i])
    
    tp = tp.groupby(['plotID'],
                    as_index=False)['W_stem',
                                    'W_branches',
                                    'W_foliage',
                                    'W_above'].apply(lambda x : x.sum()/60)
    tp['model'] = np.repeat(keys, len(tp))
    dict_AGB[keys] = pd.melt(tp, id_vars=['plotID', 'model'],
                             value_vars=['W_stem','W_branches',
                                         'W_foliage','W_above'])
del i, keys, tp

frames = list(dict_AGB.values())
df_AGB = pd.concat(frames)

#%% plotting parameters

plt.rcParams["font.family"] = "Arial"
plt.rcParams['figure.dpi'] = 300
cm = 1/2.54

xtickpos = [1,2,3,4,5]
modcolors = ["#020305", "#3B6797", "#8B2705", "#DA9652", "#FACCFA"]
modlabels = ["Magadan site", "Ust-Yansky site", "Alexander (2012), CK",
             "Kajimoto (2006), CK", "Kajimoto (2006), OM"]

complabels = ["stem", "branches", "foliage", "aboveground"]
figlabels = ["(a)", "(b)", "(c)"]
meancomplabels = [["a","a","b","bc","c"],
                  ["a","b","b","c","d"],
                  ["ab","c","a","b","d"],
                  ["a","a","b","b","b"]]

meancompva = [[0.4]*6,
              [0.03]*6,
              [0.008,0.008,0.055,0.008,0.03],
              [0.4]*6]

frames = [frames[x] for x in [3,4,0,1,2]]

#%% plotting Figure 5

data  = [frames[x][frames[x]["variable"]=="W_above"]["value"].to_list() for x in range(len(frames))]

fig1, ax1 = plt.subplots(figsize=(9*cm,9.85*cm))

bp = ax1.boxplot(data, patch_artist=True, showmeans=True)

# changing color of boxes
for patch, color in zip(bp['boxes'], modcolors):
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
for flier, color in zip(bp['fliers'], modcolors):
    flier.set_markeredgecolor(color)
    flier.set(marker='+',
              markersize=5,
              markeredgewidth=0.5)
    
ax1.set_ylabel("Predicted larch aboveground biomass (kg m$^{-2}$)",
               fontsize=7)
ax1.set_xticks(xtickpos)
ax1.set_xticklabels(modlabels, rotation = 45, fontsize=7)
ax1.tick_params(axis='y', which='major', direction='in',
                length=3, right=True, labelsize=7)
ax1.tick_params(axis='x', which='both', bottom=False, top=False)

for j in xtickpos:
    ax1.text(j, bp['whiskers'][(2*j)-1].get_ydata()[1]+meancompva[3][j-1],
             meancomplabels[3][j-1], fontsize = 8, fontstyle = 'italic')

fig1.tight_layout()

fig1.savefig(wdir+'outputs/Figures/figure5.png')

#%% plotting Figure C2

import matplotlib.gridspec as gridspec

fig2 = plt.figure(figsize=(15.9*cm,18.5*cm), constrained_layout = True)
gs = gridspec.GridSpec(2, 4, figure=fig2)
gs.update(wspace=0.5)

#ax2 ##########################################################################
ax2 = plt.subplot(gs[0, :2], )
bp2 = ax2.boxplot([frames[x][frames[x]["variable"]=="W_stem"]["value"].to_list() for x in range(len(frames))],
                  patch_artist=True, showmeans=True)

# changing color of boxes
for patch, color in zip(bp2['boxes'], modcolors):
    patch.set_facecolor(color)
#changing linewidth of boxes
for box in bp2['boxes']:
    box.set(linewidth=0.8)
#changing linewidth of whiskers
for whsk in bp2['whiskers']:
    whsk.set(linewidth=0.8)
# changing color of medians
for median in bp2['medians']:
    median.set(color ='w', linewidth=0.8)
# changing style and color of means
for mean in bp2['means']:
    mean.set(marker="*",
             markersize=5,
             markerfacecolor = "white",
             markeredgecolor = "white")
# changing style and color of fliers
for flier, color in zip(bp2['fliers'], modcolors):
    flier.set_markeredgecolor(color)
    flier.set(marker='+',
              markersize=5,
              markeredgewidth=0.5)

ax2.set_ylabel("Predicted larch stem biomass (kg m$^{-2}$)",
               fontsize=7)
ax2.set_xticks(xtickpos)
ax2.set_xticklabels(modlabels, rotation = 45, fontsize=7)
ax2.tick_params(axis='y', which='major', direction='in',
                length=3, right=True, labelsize=7)
ax2.tick_params(axis='x', which='both', bottom=False, top=False)

for j in xtickpos:
    ax2.text(j, bp2['whiskers'][(2*j)-1].get_ydata()[1]+meancompva[0][j-1],
             meancomplabels[0][j-1], fontsize = 8, fontstyle = 'italic')
    
ax2.text(4.7, 0.89*ax2.get_ylim()[1], "(a)",
         fontsize = 12, fontweight = "bold")

#ax3 ##########################################################################
ax3 = plt.subplot(gs[0, 2:])
bp3 = ax3.boxplot([frames[x][frames[x]["variable"]=="W_branches"]["value"].to_list() for x in range(len(frames))],
                  patch_artist=True, showmeans=True)

# changing color of boxes
for patch, color in zip(bp3['boxes'], modcolors):
    patch.set_facecolor(color)
#changing linewidth of boxes
for box in bp3['boxes']:
    box.set(linewidth=0.8)
#changing linewidth of whiskers
for whsk in bp3['whiskers']:
    whsk.set(linewidth=0.8)
# changing color of medians
for median in bp3['medians']:
    median.set(color ='w', linewidth=0.8)
# changing style and color of means
for mean in bp3['means']:
    mean.set(marker="*",
             markersize=5,
             markerfacecolor = "white",
             markeredgecolor = "white")
# changing style and color of fliers
for flier, color in zip(bp3['fliers'], modcolors):
    flier.set_markeredgecolor(color)
    flier.set(marker='+',
              markersize=5,
              markeredgewidth=0.5)

ax3.set_ylabel("Predicted larch branches biomass (kg m$^{-2}$)",
               fontsize=7)
ax3.set_xticks(xtickpos)
ax3.set_xticklabels(modlabels, rotation = 45, fontsize=7)
ax3.tick_params(axis='y', which='major', direction='in',
                length=3, right=True, labelsize=7)
ax3.tick_params(axis='x', which='both', bottom=False, top=False)

for j in xtickpos:
    ax3.text(j, bp3['whiskers'][(2*j)-1].get_ydata()[1]+meancompva[1][j-1],
             meancomplabels[1][j-1], fontsize = 8, fontstyle = 'italic')
    
ax3.text(4.7, 0.89*ax3.get_ylim()[1], "(b)",
         fontsize = 12, fontweight = "bold")

#ax4 ##########################################################################
ax4 = plt.subplot(gs[1, 1:3])
bp4 = ax4.boxplot([frames[x][frames[x]["variable"]=="W_foliage"]["value"].to_list() for x in range(len(frames))],
                  patch_artist=True, showmeans=True)

# changing color of boxes
for patch, color in zip(bp4['boxes'], modcolors):
    patch.set_facecolor(color)
#changing linewidth of boxes
for box in bp4['boxes']:
    box.set(linewidth=0.8)
#changing linewidth of whiskers
for whsk in bp4['whiskers']:
    whsk.set(linewidth=0.8)
# changing color of medians
for median in bp4['medians']:
    median.set(color ='w', linewidth=0.8)
# changing style and color of means
for mean in bp4['means']:
    mean.set(marker="*",
             markersize=5,
             markerfacecolor = "white",
             markeredgecolor = "white")
# changing style and color of fliers
for flier, color in zip(bp4['fliers'], modcolors):
    flier.set_markeredgecolor(color)
    flier.set(marker='+',
              markersize=5,
              markeredgewidth=0.5)
    
ax4.set_ylabel("Predicted larch foliage biomass (kg m$^{-2}$)",
               fontsize=7)
ax4.set_xticks(xtickpos)
ax4.set_xticklabels(modlabels, rotation = 45, fontsize=7)
ax4.tick_params(axis='y', which='major', direction='in',
                length=3, right=True, labelsize=7)
ax4.tick_params(axis='x', which='both', bottom=False, top=False)

for j in xtickpos:
    ax4.text(j, bp4['whiskers'][(2*j)-1].get_ydata()[1]+meancompva[2][j-1],
             meancomplabels[2][j-1], fontsize = 8, fontstyle = 'italic')
    
ax4.text(4.7, 0.89*ax4.get_ylim()[1], "(c)",
         fontsize = 12, fontweight = "bold")

fig2.savefig(wdir+'outputs/Figures/figureC2.png')
