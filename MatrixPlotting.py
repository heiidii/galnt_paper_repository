import os
import sys
import numpy as np
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
from datetime import date

oldcolor='red'
jeffcolor='black'
listres = ['L','V','A','I','M','P','G','F','W','Y','S','T','N','Q','H','K','R','D','E']
#listres1 = listres
#listres2 = listres
assert(len(listres)==19)
N= len(listres)

def add_metadata(metdadata):
  print("adding meta data") 

#palette = ['lightcoral','lightblue']
def plotdataframesingle(df,outfile=None,show=True,cbar=False,cmap='Purples',palette=['crimson','lightcoral','lavender','cornflowerblue'],title='',fmt='0.1f',mask=np.full((19,19),True,dtype=bool),vlimits=None,metadata=None,label_colorbar =None,annotate=True,listres1=listres,listres2=listres,cbar_kws_fmt=None):
  if cmap is None:
    cmap = colors.ListedColormap(palette)
  print(listres1)
  df_corr = df 
  if len(listres1)==19:
    figsize=(12,12)
    cbar_shrink = 0.5
  else:
    figsize=(12,8)
    cbar_shrink = 0.3
  fig = plt.figure(figsize=figsize)
  sns.set(rc={'xtick.labelsize': 30, 'ytick.labelsize': 30,"font.weight":'bold', 'xtick.color':jeffcolor,'ytick.color':jeffcolor,'xtick.top':False, 'xtick.labeltop':True, 'xtick.bottom':False, 'xtick.labelbottom':False})
  curaxes = plt.subplot(111)

  curaxes.set_xticklabels(listres1)#,fontsize=22,fontweight=50)
  curaxes.set_yticklabels(listres2)#,fontsize=22,fontweight=50)
  curaxes.axhline(y=0, color='k',linewidth=5)
  curaxes.axhline(y=df_corr.shape[0], color='k',linewidth=5)
  curaxes.axvline(x=0, color='k',linewidth=5)
  curaxes.axvline(x=df_corr.shape[1], color='k',linewidth=5)
  curaxes.xaxis.set_label_position('top')
  print(cmap)
  cbar_kws = {"shrink": cbar_shrink}
  if not cbar_kws_fmt is None:
    cbar_kws = {"shrink": cbar_shrink, 'format':cbar_kws_fmt}
  if vlimits is None:
      imaxes = sns.heatmap(df_corr,cmap = cmap,linewidth=0.05, square=True, ax=curaxes, cbar_kws=cbar_kws,xticklabels=1,yticklabels=1,annot=annotate,fmt=fmt,mask=mask)
  else:
      imaxes = sns.heatmap(df_corr,cmap = cmap,linewidth=0.05, square=True, ax=curaxes, cbar_kws=cbar_kws,xticklabels=1,yticklabels=1,annot=annotate,fmt=fmt,mask=mask, vmin=vlimits[0],vmax=vlimits[1])
  #IMPORTANT: needed with matplotlib 3.1.1 orw hatever is the latest version
  plt.yticks(rotation=0)
  plt.xticks(rotation=0)
  bottom, top = imaxes.get_ylim()
  imaxes.set_ylim(bottom+0.5,top-0.5)
  plt.yticks(rotation=0)
  plt.xticks(rotation=0)

  fig.suptitle(title)
  plt.tight_layout()

  if not outfile is None:
          if not metadata is None:
              add_metadata(metadata)
          print(outfile)
          plt.savefig(outfile,transparent=True,dpi=800)
  if show:
          plt.show()
  plt.close('all')
