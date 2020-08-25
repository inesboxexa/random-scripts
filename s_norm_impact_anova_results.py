# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 09:51:25 2020

@author: tm6
"""

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

data_path = '//datasvr1/MALDI_AMBIENT_DATA/Beatson/PI3K study (drug swiss)/DPO/DESI Neg (Synapt)/anova/tissue only/no norm'
data_folder = Path(data_path)
data_file = data_folder / 'anova day block vehicle 2014 8186 6244.txt' # name of the data file
data_table = pd.read_csv(data_file, delimiter='\t')

table1 = data_table.copy()
table1['norm'] = "no norm"
table1.drop_duplicates('meas mz', inplace=True)

data_path = '//datasvr1/MALDI_AMBIENT_DATA/Beatson/PI3K study (drug swiss)/DPO/DESI Neg (Synapt)/anova/tissue only/RMS'
data_folder = Path(data_path)
data_file = data_folder / 'anova day block vehicle 2014 8186 6244.txt' # name of the data file
data_table = pd.read_csv(data_file, delimiter='\t')

table2 = data_table.copy()
table2['norm'] = "RMS"
table2.drop_duplicates('meas mz', inplace=True)

data_path = '//datasvr1/MALDI_AMBIENT_DATA/Beatson/PI3K study (drug swiss)/DPO/DESI Neg (Synapt)/anova/tissue only/pqn median'
data_folder = Path(data_path)
data_file = data_folder / 'anova day block vehicle 2014 8186 6244.txt' # name of the data file
data_table = pd.read_csv(data_file, delimiter='\t')

table3 = data_table.copy()
table3['norm'] = "pqn median"
table3.drop_duplicates('meas mz', inplace=True)

data_path = '//datasvr1/MALDI_AMBIENT_DATA/Beatson/PI3K study (drug swiss)/DPO/DESI Neg (Synapt)/anova/tissue only/zscore'
data_folder = Path(data_path)
data_file = data_folder / 'anova day block vehicle 2014 8186 6244.txt' # name of the data file
data_table = pd.read_csv(data_file, delimiter='\t')

table4 = data_table.copy()
table4['norm'] = "zscore"
table4.drop_duplicates('meas mz', inplace=True)

df = pd.concat([table1, table2, table3, table4], axis=0)

df_ = pd.DataFrame()
i = 0
for p in [0, 0.0001, 0.001, 0.005, 0.01, 0.015, 0.025, 0.05 ]: #, 0.1, 0.15, 0.2]:
    
    for norm in ["no norm", "RMS", "pqn median", "zscore"]:
        
        df_.at[i,'day'] = np.sum((df['p value for day effect (mean)']<=p) & (df['norm']==norm))
        df_.at[i,'block'] = np.sum((df['p value for block effect (mean)']<=p) & (df['norm']==norm))
        df_.at[i,'vehicle'] = np.sum((df['p value for vehicle effect (mean)']<=p) & (df['norm']==norm))
        df_.at[i,'2014'] = np.sum((df['p value for 2014 effect (mean)']<=p) & (df['norm']==norm))
        df_.at[i,'8186'] = np.sum((df['p value for 8186 effect (mean)']<=p) & (df['norm']==norm))
        df_.at[i,'6244'] = np.sum((df['p value for 6244 effect (mean)']<=p) & (df['norm']==norm))
        df_.at[i,'p'] = p
        df_.at[i,'norm'] = norm
        
        i=i+1
        
#%%

sns.set(style='whitegrid', font_scale=1.35)

colors = [ 'grey', 'green', 'blue', 'red' ]

f, ax = plt.subplots(figsize=(12, 4))
sns.despine(f, left=True, bottom=True)
sns.lineplot(data=df_,
                x='p',
                y='day',
                hue='norm',
                palette = sns.xkcd_palette(colors),
                linewidth=1,
                ax=ax,
                #s=60,
                legend=False,
                zorder=1
                )
sns.scatterplot(data=df_,
                x='p',
                y='day',
                hue='norm',
                palette = sns.xkcd_palette(colors),
                linewidth=1,
                ax=ax,
                s=60,
#                legend=False,
                zorder=1
                )
ax.set_xlabel('p value')
ax.set_ylabel('# peaks')
ax.set_xlim(0,0.01)
plt.title('date of acquisition effect')

sns.set(style='whitegrid', font_scale=1.35)

f, ax = plt.subplots(figsize=(12, 4))
sns.despine(f, left=True, bottom=True)
sns.lineplot(data=df_,
                x='p',
                y='block',
                hue='norm',
                palette = sns.xkcd_palette(colors),
                linewidth=1,
                ax=ax,
                #s=60,
                legend=False,
                zorder=1
                )
sns.scatterplot(data=df_,
                x='p',
                y='block',
                hue='norm',
                palette = sns.xkcd_palette(colors),
                linewidth=1,
                ax=ax,
                s=60,
#                legend=False,
                zorder=1
                )
ax.set_xlabel('p value')
ax.set_ylabel('# peaks')
ax.set_xlim(0,0.01)
plt.title('block effect')

sns.set(style='whitegrid', font_scale=1.35)

f, ax = plt.subplots(figsize=(8, 4))
sns.despine(f, left=True, bottom=True)
sns.lineplot(data=df_,
                x='p',
                y='vehicle',
                hue='norm',
                palette = sns.xkcd_palette(colors),
                linewidth=1,
                ax=ax,
                #s=60,
                legend=False,
                zorder=1
                )
sns.scatterplot(data=df_,
                x='p',
                y='vehicle',
                hue='norm',
                palette = sns.xkcd_palette(colors),
                linewidth=1,
                ax=ax,
                s=60,
#                legend=False,
                zorder=1
                )
ax.set_xlabel('p value')
ax.set_ylabel('# peaks')
ax.set_xlim(0, 0.0505)
#ax.set_ylim(0, 200)
plt.title('vehicle effect')

sns.set(style='whitegrid', font_scale=1.35)

f, ax = plt.subplots(figsize=(8, 4))
sns.despine(f, left=True, bottom=True)
sns.lineplot(data=df_,
                x='p',
                y='2014',
                hue='norm',
                palette = sns.xkcd_palette(colors),
                linewidth=1,
                ax=ax,
                #s=60,
                legend=False,
                zorder=1
                )
sns.scatterplot(data=df_,
                x='p',
                y='2014',
                hue='norm',
                palette = sns.xkcd_palette(colors),
                linewidth=1,
                ax=ax,
                s=60,
#                legend=False,
                zorder=1
                )
ax.set_xlabel('p value')
ax.set_ylabel('# peaks')
ax.set_xlim(0, 0.0505)
#ax.set_ylim(0, 4000)
plt.title('2014 effect')

sns.set(style='whitegrid', font_scale=1.35)

f, ax = plt.subplots(figsize=(8, 4))
sns.despine(f, left=True, bottom=True)
sns.lineplot(data=df_,
                x='p',
                y='8186',
                hue='norm',
                palette = sns.xkcd_palette(colors),
                linewidth=1,
                ax=ax,
                #s=60,
                legend=False,
                zorder=1
                )
sns.scatterplot(data=df_,
                x='p',
                y='8186',
                hue='norm',
                palette = sns.xkcd_palette(colors),
                linewidth=1,
                ax=ax,
                s=60,
#                legend=False,
                zorder=1
                )
ax.set_xlabel('p value')
ax.set_ylabel('# peaks')
ax.set_xlim(0, 0.0505)
#ax.set_ylim(0, 500)
plt.title('8186 effect')

sns.set(style='whitegrid', font_scale=1.35)

f, ax = plt.subplots(figsize=(8, 4))
sns.despine(f, left=True, bottom=True)
sns.lineplot(data=df_,
                x='p',
                y='6244',
                hue='norm',
                palette = sns.xkcd_palette(colors),
                linewidth=1,
                ax=ax,
                #s=60,
                legend=False,
                zorder=1
                )
sns.scatterplot(data=df_,
                x='p',
                y='6244',
                hue='norm',
                palette = sns.xkcd_palette(colors),
                linewidth=1,
                ax=ax,
                s=60,
#                legend=False,
                zorder=1
                )
ax.set_xlabel('p value')
ax.set_ylabel('# peaks')
ax.set_xlim(0, 0.0505)
#ax.set_ylim(0, 500)
plt.title('6244 effect')