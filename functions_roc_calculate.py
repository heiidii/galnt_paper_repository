import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle

from sklearn import svm, datasets
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from scipy import interp

from matplotlib import rc

rc={'xtick.labelsize': 20, 'ytick.labelsize': 20,"font.weight":'bold', 'xtick.color':'black','ytick.color':'black','xtick.top':False, 'xtick.labeltop':False, 'xtick.bottom':True, 'xtick.labelbottom':True}

def plot_roc(fpr,tpr,filename='testroc.png',title='',str_score='',show=False):
  plt.rcParams.update(**rc)
  fig,ax = plt.subplots(figsize=(12,12))
  ax.spines['bottom'].set_color('black')
  ax.spines['top'].set_color('black')
  ax.spines['right'].set_color('black')
  ax.spines['left'].set_color('black')
  fig.suptitle(title,fontsize=25)
  plt.plot([0.0,1.0],[0.0,1.0],color='black',linestyle='dashed',linewidth=3)
  plt.plot(fpr,tpr,color='mediumvioletred',linewidth=7)
  plt.text(0.5,0.1,str_score,size=42,color='black')
  plt.ylim(0.0,1.0)
  plt.xlim(0.0,1.0)
  axis_font = { 'size':48}
  plt.xlabel("False Positive Rate (FPR)",**axis_font)
  plt.ylabel("True Positive Rate (TPR)",**axis_font)
  plt.tick_params(axis='both',labelsize=35)
  plt.tight_layout()
  plt.subplots_adjust(top=0.85)
  plt.savefig(filename,transparent=True,dpi=800)
  if show:
    plt.show()

def plot_rocs(fprs,tprs,filename='testroc.png',title='',str_scores=[],colors=[],show=False):
  plt.rcParams.update(**rc)
  fig,ax = plt.subplots(figsize=(12,12))
  ax.spines['bottom'].set_color('black')
  ax.spines['top'].set_color('black')
  ax.spines['right'].set_color('black')
  ax.spines['left'].set_color('black')
  fig.suptitle(title,fontsize=25)
  plt.plot([0.0,1.0],[0.0,1.0],color='black',linestyle='dashed',linewidth=3)
  for i,fpr in enumerate(fprs):
    if len(colors) == len(fprs):
      plt.plot(fprs[i],tprs[i],color=colors[i],linewidth=7)
    else:
      plt.plot(fprs[i],tprs[i],color='black',linewidth=7)
    if len(str_scores) == len(fprs):
      plt.text(0.5,0.05 + 0.09*i,str_scores[i],size=35,color=colors[i])
  plt.ylim(0.0,1.0)
  plt.xlim(0.0,1.0)
  axis_font = { 'size':48}
  plt.xlabel("False Positive Rate (FPR)",**axis_font)
  plt.ylabel("True Positive Rate (TPR)",**axis_font)
  plt.tick_params(axis='both',labelsize=35)
  plt.tight_layout()
  plt.subplots_adjust(top=0.85)
  plt.savefig(filename,transparent=True,dpi=800)
  if show:
    plt.show()
  plt.close()

def simplistic_roc_curve(true_labels,scores,pos_labels=1): #pos_labels 0 and 1
  '''
  Parameters:
    parameters are exactly the same as parameters for roc curve
  Returns:
    fpr , tpr, thresholds
  '''
  fpr , tpr, thresholds = roc_curve(true_labels,scores,pos_labels)  
  auc_score = roc_auc_score(true_labels,scores)
  return fpr,tpr,thresholds,auc_score
