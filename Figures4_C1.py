# -*- coding: utf-8 -*-
"""
@author: C. J. F. Delcourt
"""

#%% importing required packages
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.ticker import AutoMinorLocator
from matplotlib.lines import Line2D

## setting the working directory in which the outputs from the
## larix_allometry.R script were saved 
wdir = ""

#%% data preparation for plotting Figure C1

## importing outputs from the larix_allometry.R script
# linear regression coefficients
ols = pd.read_csv(wdir+"tableC1_ols_summary.csv")
# nonlinear regression coefficients
wls = pd.read_csv(wdir+"table5_wls_summary.csv")
# weightings of the nonlinear regressions (c values)
cvals = pd.read_csv(wdir+"tableB1_wls_weightings.csv")
# tree biomass and size measurements at Magadan and Ust-Yansky sites
tree_db = pd.read_csv(wdir+"Larch_tree_DB.csv")

## storing linear regression coefficients for each site in a dictionnary
ols_dict = {"MG": {"a": [ols.iloc[5*i,1] for i in range(6)],
                   "b": [ols.iloc[5*i+1,1] for i in range(6)],
                   "rmse": [ols.iloc[5*i+3,1] for i in range(6)]},
            "UY": {"a": [ols.iloc[5*i,2] for i in range(6)],
                   "b": [ols.iloc[5*i+1,2] for i in range(6)],
                   "rmse": [ols.iloc[5*i+3,2] for i in range(6)]}}

## storing nonlinear regression coefficients for each site in a dictionnary
wls_dict = {"MG": {"a": wls.iloc[:6,5].to_list(),
                   "b": wls.iloc[:6,7].to_list(),
                   "rmse": wls.iloc[:6,9].to_list(),
                   "c": cvals["MG.c"].to_list()},
            "UY": {"a": wls.iloc[6:,5].to_list(),
                   "b": wls.iloc[6:,7].to_list(),
                   "rmse": wls.iloc[6:,9].to_list(),
                   "c": cvals["UY.c"].to_list()}}

#%% plotting Figure C1

# site-specific allometry models for Larix cajanderi developed using linear
# regressions (OLS) and weighted nonlinear regressions (WLS).

plt.rcParams["font.family"] = "Arial"
plt.rcParams['figure.dpi'] = 300
cm = 1/2.54

x = np.linspace(0, 55, 500)
ylims = [1850, 1850, 225, 175, 35, 1850]
ylabels = ["Stem", "Stem wood", "Stem bark", "Branches",
           "Foliage", "Total aboveground"]
yfiles = ["stem", "stemwood", "stembark", "branches",
          "foliage", "aboveground"]
figlabels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]

custom_lines = [Line2D([0], [0], color="#020305",
                       linestyle="dotted", lw=2),
                Line2D([0], [0], color="#020305",
                       linestyle="solid", lw=2),
                Line2D([0], [0], marker='.', color='w',
                       markerfacecolor="#020305",
                       markeredgecolor="black", markersize=7),
                Line2D([0], [0], color="w", lw=2),
                Line2D([0], [0], color="#3B6797",
                       linestyle=(0, (2, 1)), lw=2),
                Line2D([0], [0], color="#3B6797",
                       linestyle="solid", lw=2),
                Line2D([0], [0], marker='D', color='w',
                       markerfacecolor="#3B6797",
                       markeredgecolor="#000000",
                       markeredgewidth = 0.5, markersize=3.5)]

lgd_lines = [Line2D([0], [0], color="#020305",
                    linestyle="dotted", lw=2.5),
             Line2D([0], [0], color="#020305",
                    linestyle="solid", lw=2.5),
             Line2D([0], [0], marker='.', color='w',
                    markerfacecolor="#020305",
                    markeredgecolor="black", markersize=10),
             Line2D([0], [0], color="#3B6797",
                    linestyle=(0, (2, 1)), lw=2.5),
             Line2D([0], [0], color="#3B6797",
                    linestyle="solid", lw=2.5),
             Line2D([0], [0], marker='D', color='w',
                    markerfacecolor="#3B6797",
                    markeredgecolor="#000000",
                    markeredgewidth = 0.5, markersize=5)]

figlgd_labels = ["OLS (Magadan site)",
                 "WLS (Magadan site)",
                 "Biomass measurements (Magadan site)",
                 "OLS (Ust-Yansky site)",
                 "WLS (Ust-Yansky site)",
                 "Biomass measurements (Ust-Yansky site)"]

fig, axs = plt.subplots(3, 2, figsize=(16*cm,18.7*cm))

lgd2 = fig.legend(lgd_lines, figlgd_labels, ncol = 2, fontsize = 7.5,
                  loc = "lower center", bbox_to_anchor = (0.53, 0),
                  frameon = False, handletextpad = 1.4,
                  handlelength = 2.75)

plt.gca().add_artist(lgd2)

for ax, i in zip(axs.ravel(), range(6)):
    
    # OLS (log-log tansfromation)
    y1 = ols_dict["MG"]["a"][i]*(x**ols_dict["MG"]["b"][i])
    y2 = ols_dict["UY"]["a"][i]*(x**ols_dict["UY"]["b"][i])
    
    # WLS (DBH-2c)
    y3 = wls_dict["MG"]["a"][i]*(x**wls_dict["MG"]["b"][i])
    y4 = wls_dict["UY"]["a"][i]*(x**wls_dict["UY"]["b"][i])
    
    
    ax.plot(x, y1, label = "OLS (Magadan site)",
            linestyle = "dotted", lw = 1.5, color = "#020305", alpha = 0.5)
    ax.plot(x, y2, label = "OLS (Ust-Yansky site)",
            linestyle = (0, (5, 5)), lw = 1.5, color = "#3B6797", alpha = 0.5)

    ax.plot(x, y3, label = "WLS (Magadan site)",
            linestyle = "solid", lw = 1.5, color = "#020305")
    ax.plot(x, y4, label = "WLS (Ust-Yansky site)",
            linestyle = "solid", lw = 1.5, color = "#3B6797")

    ax.plot(tree_db[tree_db["Location"]=="MG"]["DBH"],
            tree_db[tree_db["Location"]=="MG"].iloc[:,i],
            '.', markersize=8, markerfacecolor="#020305",
            markeredgecolor="black", markeredgewidth=0.5)
    ax.plot(tree_db[tree_db["Location"]=="UY"]["DBH"],
            tree_db[tree_db["Location"]=="UY"].iloc[:,i],
            'D', markersize=4, markerfacecolor="#3B6797",
            markeredgecolor="black", markeredgewidth=0.5)
    
    ax.set_xlim(0,55)
    ax.set_ylim(0, ylims[i])
    ax.set_xlabel("DBH (cm)", fontsize = 7)
    ax.set_ylabel(ylabels[i]+' biomass (kg)', fontsize = 7)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', which='minor', direction='in',
                   length=2, top=True)
    ax.tick_params(axis='both', which='major', direction='in',
                   length=3, top=True, labelsize=7)
    
    ax.text(46, 0.86*ylims[i], figlabels[i],
            fontsize = 12, fontweight = "bold")
    
    lgd_labels = ["RMSE = "+\
                  str(round(ols_dict["MG"]["rmse"][i], 2))+" kg",
                  "RMSE = "+\
                  str(round(wls_dict["MG"]["rmse"][i], 2))+" kg",
                  "($\it{n}$ ="+\
                      str(cvals["MG.n"][i])+")",
                  "",
                  "RMSE = "+\
                  str(round(ols_dict["UY"]["rmse"][i], 2))+" kg",
                  "RMSE = "+\
                  str(round(wls_dict["UY"]["rmse"][i], 2))+" kg",
                  "($\it{n}$ ="+\
                      str(cvals["UY.n"][i])+")",
                  ""]
    
    ax.legend(custom_lines, lgd_labels, fontsize = 5)


fig.tight_layout()
fig.subplots_adjust(bottom = 0.13, wspace = 0.37)
    
fig.savefig(wdir+'Figures/figureC1.png')

del ax, axs, custom_lines, fig, figlabels, figlgd_labels, cm, ylabels
del i, lgd2, lgd_labels, lgd_lines, x, y1, y2, y3, y4, ylims, yfiles

#%% importing confidence intervals predictions from the larix_allometry.R script

## setting the directory in which the confidence intervals were saved
os.chdir(wdir+'/Figures/Figure4_CIs/')
    
## specifying an empty dictionnary for content
cis_dict = {'MG': [], 'UY': []}
for ref in cis_dict.keys():
    for file in [[x for x in os.listdir() if x.startswith('CIs_'+ref)][i] for i in [5,3,4,0]]:
        ## reading content into data frame
        df = pd.read_csv(file)
        cis_dict[ref].append(df)
del(df, file, ref)

#%% existing allometric models for Larix cajanderi in northeast Siberia

# for allometric relationships developed in Chersky (CK) and Oymyakon (OM)
# areas, regression coefficients can be found in Kajimoto et al. (2006)
# and Alexander et al. (2012)
lrx_allo = {"alexander_CK": {"a": [0.08142,0.06966,0.0405,0.1792],
                             "b": [2.1,1.99,1.41,2.01]},
            "kajimoto_CK": {"a": [0.113,0.0396,0.0106],
                            "b": [2.04,1.91,1.94]},
            "kajimoto_OM": {"a": [0.14,0.0596,0.0186],
                            "b": [1.96,1.48,1.31]}}

## selecting only stem, branches, foliage and total aboveground biomass for
## comparisons with existing biomass equations
tree_db = tree_db[['Pst', 'Pbr', 'Pf', 'Pabo', 'DBH', 'Location']]

#%% plotting Figure 4

plt.rcParams["font.family"] = "Arial"
plt.rcParams['figure.dpi'] = 300
plt.rcParams["legend.title_fontsize"] = 20
cm = 1/2.54

x = np.linspace(0,30,500)
idx_comp = [0,3,4,5]
ylims = [375, 65, 20, 420]
ylabels = ["Stem", "Branches", "Foliage", "Total aboveground"]
yfiles = ["stem", "branches", "foliage", "aboveground"]
figlabels = ["(a)", "(b)", "(c)", "(d)"]
colors = ["#020305", "#3B6797", "#8B2705", "#DA9652", "#FACCFA"]

fillcolors = {'MG': matplotlib.colors.to_rgb(colors[0])+(0.2,),
              'UY': matplotlib.colors.to_rgb(colors[1])+(0.2,)}

custom_lines = [Line2D([0], [0], color="#020305",
                       linestyle = "solid", lw=8),
                Line2D([0], [0], color="#3B6797",
                       linestyle = (0, (3.5, 2)), lw=8),
                Line2D([0], [0], color="#8B2705",
                       linestyle="dotted", lw=8),
                Line2D([0], [0], color="#DA9652",
                       linestyle = (0, (2, 2, 0.75, 2)), lw=8),
                Line2D([0], [0], color="#FACCFA",
                       linestyle = (0, (1.75, 1.25, 0.75, 1.25, 0.75, 1.25)),
                       lw=8)]

lgd1_labels = ["Magadan site ($\it{n}$ = 43, DBH: 3.9-52.8 cm)",
               "Ust-Yansky site ($\it{n}$ = 20, DBH: 1.8-18.9 cm)",
               "Alexander et al. (2012), CK ($\it{n}$ = 32, DBH: 0.08-29.3 cm)",
               "Kajimoto et al. (2006), CK ($\it{n}$ = 7, max. DBH: 18.6 cm)",
               "Kajimoto et al. (2006), OM ($\it{n}$ = 6, max. DBH: 16.1 cm)"]

custom_markers = [Line2D([0], [0], marker = '.', color="w",
                         markerfacecolor = "#020305",
                         markeredgecolor = 'black', markersize = 30),
                  Line2D([0], [0], marker = 'D', color="w",
                         markerfacecolor = "#3B6797",
                         markeredgecolor = 'black', markersize = 14)]

lgd2_labels = ["Magadan site", "Ust-Yansky site"]

fig, axs = plt.subplots(2, 2, figsize=(16, 16.5))

for ax, i in zip(axs.ravel(), range(4)):
  
    y1 = wls_dict["MG"]["a"][idx_comp[i]]*(x**wls_dict["MG"]["b"][idx_comp[i]])
    y2 = wls_dict["UY"]["a"][idx_comp[i]]*(x**wls_dict["UY"]["b"][idx_comp[i]])
    y3 = lrx_allo["alexander_CK"]["a"][i]*(x**lrx_allo["alexander_CK"]["b"][i])
    
    if i == 3:
        y4 = lrx_allo["kajimoto_CK"]["a"][i-3]*(x**lrx_allo["kajimoto_CK"]["b"][i-3]) +\
            lrx_allo["kajimoto_CK"]["a"][i-2]*(x**lrx_allo["kajimoto_CK"]["b"][i-2]) +\
                lrx_allo["kajimoto_CK"]["a"][i-1]*(x**lrx_allo["kajimoto_CK"]["b"][i-1])
            
        y5 = lrx_allo["kajimoto_OM"]["a"][i-3]*(x**lrx_allo["kajimoto_OM"]["b"][i-3]) +\
            lrx_allo["kajimoto_OM"]["a"][i-2]*(x**lrx_allo["kajimoto_OM"]["b"][i-2]) +\
                lrx_allo["kajimoto_OM"]["a"][i-1]*(x**lrx_allo["kajimoto_OM"]["b"][i-1])
            
    else:
        y4 = lrx_allo["kajimoto_CK"]["a"][i]*(x**lrx_allo["kajimoto_CK"]["b"][i])
        y5 = lrx_allo["kajimoto_OM"]["a"][i]*(x**lrx_allo["kajimoto_OM"]["b"][i])
        
    ax.plot(x, y1, label = "Magadan site",
            linestyle = 'solid', lw = 3.5, color = colors[0])
    ax.plot(x, y2, label = "Ust-Yansky site",
            linestyle = (0, (5, 5)), lw = 3.5, color = colors[1])
    ax.plot(x, y4, label = "Kajimoto et al. (2006) - Chersky",
            linestyle = (0, (3, 5, 1, 5)), lw = 3.5, color = colors[3])
    ax.plot(x, y5, label = "Kajimoto et al. (2006) - Oymyakon",
            linestyle = (0, (3, 5, 1, 5, 1, 5)), lw = 3.5, color = colors[4])
    ax.plot(x, y3, label = "Alexander et al. (2012)",
            linestyle = "dotted", lw = 3.5, color = colors[2])
    ax.plot(tree_db[tree_db["Location"]=="MG"]["DBH"],
            tree_db[tree_db["Location"]=="MG"].iloc[:,i],
            '.', markersize=16, markerfacecolor=colors[0],
            markeredgecolor='black', markeredgewidth=1)
    ax.plot(tree_db[tree_db["Location"]=="UY"]["DBH"],
            tree_db[tree_db["Location"]=="UY"].iloc[:,i],
            'D', markersize=8, markerfacecolor=colors[1],
            markeredgecolor='black', markeredgewidth=1)
    ax.fill_between(cis_dict['UY'][i]['DBH'],
                    cis_dict['UY'][i]['conf_lower'],
                    cis_dict['UY'][i]['conf_upper'],
                    facecolor=fillcolors['UY'], edgecolor=colors[1],
                    linewidth=1.5, linestyle = (0, (5, 5)))
    ax.fill_between(cis_dict['MG'][i]['DBH'],
                    cis_dict['MG'][i]['conf_lower'],
                    cis_dict['MG'][i]['conf_upper'],
                    facecolor=fillcolors['MG'], edgecolor=colors[0],
                    linewidth=1.5)
    
    ax.set_xlim(0, 30)
    ax.set_ylim(0, ylims[i])
    ax.set_xlabel("DBH (cm)", fontsize = 18, labelpad = 10)
    ax.set_ylabel(ylabels[i]+' biomass (kg)', fontsize = 18, labelpad = 10)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='both', which='minor', direction='in',
                   length=8, width = 2, top=True)
    ax.tick_params(axis='both', which='major', direction='in',
                   length=11, width = 2, top=True, labelsize=18, pad = 10)
    ax.text(2.2, 0.86*ylims[i], figlabels[i], fontsize = 34,
            fontweight = "bold")
    
    [x.set_linewidth(2) for x in ax.spines.values()]

lgd2 = plt.legend(custom_markers, lgd2_labels, fontsize = 20,
                  loc = "center left", bbox_to_anchor = (0, -0.361),
                  borderaxespad = -26, frameon = False,
                  title = "in situ biomass measurements",
                  labelspacing = 1, handletextpad = 1.4)

title2 = lgd2.get_title()

lgd2._legend_box.align = "left"
title2.set_weight("bold")
ax = plt.gca().add_artist(lgd2)

lgd1 = plt.legend(custom_lines, lgd1_labels, fontsize = 20,
                  loc = "center left", bbox_to_anchor = (0,-0.53),
                  borderaxespad = -9, frameon = False,
                  title = "Allometric relationships",
                  labelspacing = 1, handlelength = 3.5,
                  handletextpad = 1.4)
title1 = lgd1.get_title()

lgd1._legend_box.align = "left"
title1.set_weight("bold")
    
fig.tight_layout()
fig.subplots_adjust(hspace = .2, wspace = .2)

fig.savefig(wdir+'Figures/figure4.png')
