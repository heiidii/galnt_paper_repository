import os
import sys
import pandas as pd
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import numpy as np
import numpy.ma as ma
import seaborn as sns
import pickle
from MatrixPlotting import *
from datetime import date
import glob
from functions_scores_to_matrix_3 import makematrixfordict,makedataframefordict
from MultivariablePlottingGlobal import *
import math 
d = date.today()
datestring = d.isoformat()
res_dict = {498:'pos -1',499:'pos 0',500:'pos +1'}
key_energies_res_res = 'dataframe_pairwise_energies'
key_energies_res = 'dataframe_total_res_energies'

def setupaxes(ax1):
  ax1.tick_params(axis='x',labelsize=35)#, grid_linewidth=0.5, grid_linestyle='-', grid_color='black',grid_alpha=0.6)
  ax1.tick_params(axis='y', labelsize=35)#labelcolor='black',labelsize=25, grid_linewidth=0.5, grid_linestyle='-', grid_color='black',grid_alpha=0.6)
  #ax1.grid(linestyle='-', linewidth='0.5', color='gray')
  ax1.xaxis.set_label_text("")


def get_filtered_dataframe(df,criteria, rmsd_cutoff):
  if criteria=='rmsd_greater_than':
    df_filter = df[(df[labels_df[info['xfield']]] > rmsd_cutoff)]
    return df_filter
  if criteria=='rmsd_less_than':
    df_filter = df[(df[labels_df[info['xfield']]] < rmsd_cutoff)]
    return df_filter
  if criteria=='sc':
    df_filter = df[(df['sc_shapecomplementarity'] > rmsd_cutoff)]
    return df_filter

def get_stats(df,prop,key,dict_stats):
    
    dict_stats[key]['mean'] = df[prop].mean()
    dict_stats[key]['median'] = df[prop].median()
    dict_stats[key]['quartile_lower'] = df[prop].quantile(0.25)
    dict_stats[key]['quartile_upper'] = df[prop].quantile(0.75)
    dict_stats[key]['sd'] = df[prop].std()


def read_pickled_file(pfile):

    df_energies_res = None
    df_energies_res_res = None        
    if not os.path.exists(pfile): return
    data_ = pickle.load(open(pfile,'rb'))
    df = data_['dataframe']
    
    if key_energies_res in data_:
      df_energies_res = data_[key_energies_res]
    if key_energies_res_res in data_:
      df_energies_res_res = data_[key_energies_res_res]

    return df, df_energies_res, df_energies_res_res

def get_dict_from_doubledict(dict_stats,skey):
    prim_key = dict_stats.keys()
    dicttemp={}
    listexp=[]
    listtemp=[]
    for key in prim_key:
      value = dict_stats[key][skey]
      if not math.isnan(value):
        dicttemp[key] = value#dict_stats[key][skey]
        listtemp.append(value)
        listexp.append(dict_stats[key]['exp2018'])
    return dicttemp,listtemp,listexp

def get_stat(temp_,prop,aa,dict_pos_,by='median'):
        if by=='median':
          dict_pos_[aa]= temp_[prop].median(skipna=True)
        if by=='quartile_lower':
          dict_pos_[aa]= temp_[prop].quantile(0.25)
        if by=='quartile_upper':
          dict_pos_[aa]= temp_[prop].quantile(0.75)
        if by=='min':
          dict_pos_[aa]= temp_[prop].min(skipna=True)
        if by=='sd':
          dict_pos_[aa]= temp_[prop].std(skipna=True)
        if by=='mean':
          dict_pos_[aa]= temp_[prop].mean(skipna=True)

def get_stat_periodic(temp_,prop,aa,dict_pos_,by='sd'):
        from scipy.stats import circstd
        vals = list(temp_[prop])
        if by=='sd':
          dict_pos_[aa]= circstd(temp_[prop],high=np.pi,low=-np.pi)

def get_sorted_df(pos,df,prop,by='median'):
    dict_pos_={}
    unique_aa = set(list(df[pos].values))
    print(unique_aa)
    for aa in unique_aa:
        temp_ = df[ df[pos].str.match(aa)]
        #if prop.find('phi')!=-1 or prop.find('phi')!=-1 or prop.find('omega')!=-1:
        #  dict_pos_ = get_stat_periodic(temp_,prop,aa,by=by)
        #else:
        get_stat(temp_,prop,aa,dict_pos_,by=by)
    print(dict_pos_)
    from collections import OrderedDict
    dd = OrderedDict(sorted (dict_pos_.items(), key= lambda x:x[1]))
    print(dd)
    sorter = dd.keys()
    df[pos]=df[pos].astype('category')
    df[pos].cat.set_categories(sorter,inplace=True)
    print(df[pos])
    tempdf = df.sort_values([pos])
    tempdf = tempdf.reset_index(drop=True)
    return tempdf

def get_positionwise_plots(dict_stats,skey,prop,criteria,rmsd_cutoff,N,basename,show=False,suffix=None,reverse=False, applyfraction=False,maxmin=None,annotate=False,pointplot=False,positions=['-1','+1'],by='median',boxplot=True):
    dicttemp,listtemp,listexp = get_dict_from_doubledict(dict_stats,skey)
    value_name = prop
    if skey=="sd":
      value_name= prop+'_'+skey
    from functions_matrices_from_clusters import plot_roc_auc_for_dict, get_default_value, cmap_for_field
    default_value = get_default_value(value_name)
    df =  makedataframefordict(dicttemp,prop)
    tempdata = np.array(list(dicttemp.values()))
    name_suffix = '%s_%s_%f_%s_top%d' %(prop,criteria,rmsd_cutoff,skey,N)
    if not suffix is None:
      name_suffix += '_'+suffix
    if applyfraction:
      name_suffix += '_fractioncutoff'
    if reverse:
      name_suffix += 'reverse'
    title = '%s\n %s' %(datestring,name_suffix)
    if pointplot:
      base='pointplot'
      positions = ['+1']
    elif boxplot:
      base='boxplot'
      positions = positions
    else:
      base='swarmplot'
      positions = ['+1']
    for pos in positions:
      print(pos)
      outfile = '%s/%s_%s_pos%s_by%s.png' %(basename,base,name_suffix,pos,by)
      print( df)
      print(df.columns.values)
      #exit()
      tempdf = get_sorted_df(pos,df,prop,by=by)
      print(tempdf)
      fig,ax = plt.subplots(figsize=[10,7])
      if pointplot:
        ax = sns.pointplot(y=prop,x=pos,data=tempdf,ax=ax,color='gainsboro',size=12) #mediumseagreen
      elif boxplot:
        ax = sns.boxplot(y=prop,x=pos,data=tempdf,ax=ax,color='gainsboro') #mediumseagreen
        ax =sns.swarmplot(y=prop,x=pos,data=tempdf,color=".2",size=8)
      else:
        ax =sns.swarmplot(y=prop,x=pos,data=tempdf,color=".2",size=8,palette='Set2')
      setupaxes(ax)
      plt.tight_layout()
      print(outfile)
      plt.savefig(outfile,transparent=True,dpi=1200)
      if show:
        plt.show()
      plt.close()

def get_positionwise_plots_2(df,prop,criteria,rmsd_cutoff,N,basename,show=False,suffix=None,reverse=False, applyfraction=False,maxmin=None,annotate=False,write_pdbs=True,plot=True,by='median',pos='+1'):
    #print(df)
    value_name = prop
    name_suffix = '%s_%s_%f_top%d' %(prop,criteria,rmsd_cutoff,N)
    if not suffix is None:
      name_suffix += '_'+suffix
    if applyfraction:
      name_suffix += '_fractioncutoff'
    if reverse:
      name_suffix += 'reverse'
    title = '%s\n %s' %(datestring,name_suffix)
    base='boxplot2'
    outfile = '%s/%s_%s_pos%s.png' %(basename,base,name_suffix,pos)
    if pos=='+1':
      df[pos]=[t[3] for t in list(df['key'])]
    if pos=='-1':
      df[pos]=[t[1] for t in list(df['key'])]
    print('dfpos ',pos,df[pos].values)

    tempdf = get_sorted_df(pos,df,prop,by=by)
    #print(tempdf)
    if write_pdbs:
      outfile_pdb = '%s/pdbfiles_%s_pos%s.txt' %(basename,name_suffix,pos)
      outfp=open(outfile_pdb,'w')
      for i,entry in enumerate(list(tempdf['description'])):
        outfp.write(str(entry)+'\t'+str(tempdf[prop][i])+'\t'+str(tempdf['rmsd'][i])+'\t'+str(tempdf['IE'][i])+'\n')
      outfp.close()  
    if plot:
      outfile = '%s/%s_%s_pos%s.png' %(basename,base,name_suffix,pos)
      fig,ax = plt.subplots(figsize=[12,8])
      ax = sns.boxplot(y=prop,x=pos,data=tempdf,ax=ax,color='lightseagreen')
      ax =sns.swarmplot(y=prop,x=pos,data=tempdf,color=".2",size=5)
      setupaxes(ax)
      if not maxmin is None:
        ax.set(ylim=maxmin)
      plt.tight_layout()
      plt.savefig(outfile,transparent=True,dpi=1200)
      if show:
        plt.show()
      plt.close()

def setup_and_plot_heatmap_from_dict(dict_stats,skey,prop,criteria,rmsd_cutoff,N,basename,show=False,suffix=None,reverse=False, applyfraction=False,maxmin=None,annotate=False,listres1=aalist):

    dicttemp,listtemp,listexp = get_dict_from_doubledict(dict_stats,skey)
    value_name = prop
    if skey=="sd":
      value_name= prop+'_'+skey
    from functions_matrices_from_clusters import plot_roc_auc_for_dict, get_default_value, cmap_for_field
    default_value = get_default_value(value_name)
    mat = np.full((len(listres1),19), default_value ,dtype=float)
    tempdf,mat =  makematrixfordict(dicttemp,mat,listres_1=listres1) 
    #df,mat =  makematrixfordict(dictmat,mat,listres_1=listres1)
    if len(list(dicttemp.values()))<5: return
    tempdata = np.array(list(dicttemp.values()))
    mask = (mat == default_value)
    name_suffix = '%s_%s_%f_%s_top%d' %(prop,criteria,rmsd_cutoff,skey,N)
    if not suffix is None:
      name_suffix += '_'+suffix
    if applyfraction:
      name_suffix += '_fractioncutoff'
    if reverse:
      name_suffix += 'reverse'
    title = '%s\n %s' %(datestring,name_suffix)
    outfile = '%s/heatmap_%s.png' %(basename,name_suffix)
    if maxmin is None:
      vmin = np.quantile(tempdata,0.05)
      vmax = np.quantile(tempdata,0.95)
    else:
      vmin = maxmin[0]
      vmax = maxmin[1]
    #plot_overalltrends(None,np.array(listtemp),np.array(listexp),show=True)

    plotdataframesingle( tempdf ,cmap = cmap_for_field(value_name) , cbar=True , title=title , mask=mask, vlimits=[float(vmin),float(vmax)], outfile =outfile ,show=show,annotate=annotate,fmt='0.2f',listres1=listres1)

def get_unique_full_list(listofdfs,column='residue_j',filterres=None):
  longlist = []
  for df in listofdfs:
    df_ = df
    if not filterres is None:
      df_ = df [ df.residue_i.eq(filterres) ]
    for entry in list(df_[column]):
      longlist.append(entry)
  return list(set(longlist))

def get_scoretype_dict_for_residue_pair(listofdfs_energies_res_res,listofdfs,res_i,res_j,st,filter_by_pdbs=True,prop = 'score_value',pairwise=True):
        dict_stats = {}
        for idf,df in enumerate(listofdfs_energies_res_res):
          if pairwise:
            temp_df_f = df [ df.residue_i.eq(res_i) & df.residue_j.eq(res_j) & df.score_type.eq(st)]
          else:
            temp_df_f = df [ df.residue.eq(res_i) &  df.score_type.eq(st)]
          temp_df_f.reset_index()
          if temp_df_f.shape[0] < 1:
            #print('not enough samples in energy dataframe 1 ',df.shape,temp_df_f.shape)
            continue
          parent_df = listofdfs[idf] 
          listpdbs_cutoff = list(set(parent_df['description'])) #unique list
          if len(listpdbs_cutoff)==0:
            #print('empty filtered parent dataframe ')
            continue
          if filter_by_pdbs:
            temp_df = temp_df_f [ temp_df_f['file'].isin( listpdbs_cutoff ) ]
            temp_df.reset_index()
          else:
            temp_df = temp_df_f

          if temp_df is None:
            continue
          if temp_df.shape[0] < 1:
            #print('not enough samples in energy dataframe 2')
            continue

          listkey = list(parent_df['key'])
          key = listkey[0]
          listexp = list(parent_df['exp2018'])
          exp2018 = listexp[0]
          dict_stats[key] ={}
          dict_stats[key]['exp2018']=exp2018
          get_stats(temp_df_f[:N],prop,key,dict_stats)
        
        return dict_stats

def get_dict_stat_for_res_pair(pickledfiles='PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_Energies*.p',basename='results/heatmaps_energies_pairwise_filtered/',N=20,rmsd_cutoff=1.0,criteria='rmsd_less_than',res_i = 498,res_j=291,st='fa_atr',skey='median',off=0.10,plotdist=False,reverse=False,applyfraction=False,prop = 'score_value',pairwise=True):
  listofdfs =[]
  listofdfs_nofilter = []
  listofdfs_energies_res = []
  listofdfs_energies_res_res = []
  pfiles = glob.glob(pickledfiles)
  dict_stats = {}
  dict_stats_size = {}
  pfiles.sort()
  for pfile in pfiles:
    print(pfile)
    df, df_energies_res, df_energies_res_res = read_pickled_file(pfile)

    key = df['key'][0]
    exp = df['exp2018'][0]

    df_filter = get_filtered_dataframe(df,criteria, rmsd_cutoff)
    if df_filter is None:
      continue
    total_size = float(len(df))
    Nsize = len(df_filter)/float(total_size)
    if applyfraction:
      if Nsize < 0.10:
        continue
    if not reverse:
      df_filter_top = df_filter[:N]
    else:
      df_top = df[:N]
      df_filter_top = get_filtered_dataframe(df,criteria, rmsd_cutoff)

    if df_filter_top is None:
        continue

    listofdfs.append(df_filter_top)
    listofdfs_nofilter.append(df)
    listofdfs_energies_res.append(df_energies_res)
    listofdfs_energies_res_res.append(df_energies_res_res)

  dict_stats = get_scoretype_dict_for_residue_pair(listofdfs_energies_res_res,listofdfs,res_i,res_j,st,pairwise=pairwise)
  #setup_and_plot_heatmap_from_dict(dict_stats,skey,prop,criteria,rmsd_cutoff,N,basename,show=True,suffix=None,reverse=False, applyfraction=False,maxmin=None,annotate=False)

  dicttemp,_,_ =  get_dict_from_doubledict(dict_stats,skey)
  return dicttemp

def process_res_energies_to_heatmaps(listofdfs,listofdfs_energies_res,criteria='',cutoff=100.0,N=10,basename='./',show=False,filter_by_pdbs=True,pepres=None,e_terms = None,plotheatmaps=True):
  df_0 = listofdfs_energies_res[0]
  if pepres is None:
    list_residue_i = list(set(list(df_0['residue'])))
  else:
    list_residue_i = pepres
  if e_terms is None:
    score_types = list(set(list(df_0['score_type'])))
  else:
    score_types = e_terms

  #ierate over i-j pairs and score terms
  prop = 'score_value'
  for res_i in list_residue_i:
    print(res_i)
    for st in score_types:
        print(st)
        dict_stats = get_scoretype_dict_for_residue_pair(listofdfs_energies_res,listofdfs,res_i,288,st,pairwise=False)
        prim_key = list(dict_stats.keys())
        print(prim_key, len(prim_key))
        if len(prim_key)<5:
          print("empty dict")
          continue
        key = prim_key[0]
        sec_keys =  dict_stats[key].keys()
        for skey in sec_keys:
          if skey == 'exp2018': continue
          suffix = '%d_%s' %(res_i,st)
          if plotheatmaps:
            setup_and_plot_heatmap_from_dict(dict_stats,skey,prop,criteria,cutoff,N,basename,show=show,suffix=suffix)


def process_pairwise_energies_to_heatmaps(listofdfs,listofdfs_energies_res_res,criteria='',cutoff=100.0,N=10,basename='./',show=False,filter_by_pdbs=True,pepres=[498],enzymeres=[291],e_terms = None,plotheatmaps=True,listres1=aalist):
 
  #Use 0th dataframe to determine scoretypes and residue_j
  df_0 = listofdfs_energies_res_res[0]

  if pepres is None:
    list_residue_i = list(set(list(df_0['residue_i'])))
  else:
    list_residue_i = pepres

  if e_terms is None:
    score_types = list(set(list(df_0['score_type']))) 
  else:
    score_types = e_terms

  #ierate over i-j pairs and score terms
  prop = 'score_value'
  for res_i in list_residue_i:
    print(res_i)
    if enzymeres is None:
      list_residue_j = get_unique_full_list(listofdfs_energies_res_res,filterres = res_i)
    else:
      list_residue_j = enzymeres
    print('residue js: ',list_residue_j)
    for res_j in list_residue_j:
      for st in score_types:
        print(st)
        dict_stats = get_scoretype_dict_for_residue_pair(listofdfs_energies_res_res,listofdfs,res_i,res_j,st,pairwise=True)
        prim_key = list(dict_stats.keys())
        print(prim_key, len(prim_key))
        if len(prim_key)<5:
          print("empty dict")
          continue
        key = prim_key[0]
        sec_keys =  dict_stats[key].keys()
        for skey in sec_keys:
          if skey == 'exp2018': continue
          suffix = '%d-%d_%s' %(res_i,res_j,st)
          if plotheatmaps:
            setup_and_plot_heatmap_from_dict(dict_stats,skey,prop,criteria,cutoff,N,basename,show=show,suffix=suffix,listres1=listres1)

def plot_heatmaps_from_pickled_dataframe_energies(pickledfiles='PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_Energies*.p',basename='results/heatmaps_energies_pairwise_filtered/',N=5,rmsd_cutoff=1.0,criteria='rmsd_less_than',cmap=None,off=0.00,plotdist=False,fitstuff=False,reverse=False,maxmin=None,annotate=False,suffix=None,applyfraction=False,show=False,donotapplytopNtoparent=True,pairwise=True,listres1=aalist):

  os.system('mkdir -p %s' %basename)
  if cmap is None:
    cmap = 'Greens'
  listofdfs =[]
  listofdfs_nofilter = []
  listofdfs_energies_res = []
  listofdfs_energies_res_res = []
  pfiles = glob.glob(pickledfiles)
  dict_stats = {}
  dict_stats_size = {}
  pfiles.sort()
  for pfile in pfiles:
    print(pfile)
    df, df_energies_res, df_energies_res_res = read_pickled_file(pfile)

    key = df['key'][0]
    exp = df['exp2018'][0]

    df_filter = get_filtered_dataframe(df,criteria, rmsd_cutoff)
    if df_filter is None:
      continue
    total_size = float(len(df))
    Nsize = len(df_filter)/float(total_size)
    if applyfraction:
      if Nsize < 0.10:
        continue
    if not reverse:
      df_filter_top = df_filter[:N]
    else:
      df_top = df[:N]
      df_filter_top = get_filtered_dataframe(df,criteria, rmsd_cutoff)

    if df_filter_top is None:
        continue
    if donotapplytopNtoparent:
      df_filter_top = df_filter

    listofdfs.append(df_filter_top)
    listofdfs_nofilter.append(df)
    listofdfs_energies_res.append(df_energies_res)
    listofdfs_energies_res_res.append(df_energies_res_res)

  if pairwise:
    process_pairwise_energies_to_heatmaps(listofdfs,listofdfs_energies_res_res,criteria=criteria,cutoff=rmsd_cutoff,N=N,basename=basename,show=show,listres1=listres1)
  else:
    process_res_energies_to_heatmaps(listofdfs,listofdfs_energies_res,criteria=criteria,cutoff=rmsd_cutoff,N=N,basename=basename,show=show)

def write_metrics(keys_true_positives,keys_true_negatives,keys_false_positives,keys_false_negatives,file_handle=None):
  tp=len(keys_true_positives)
  tn=len(keys_true_negatives) #2 cases where the data is not entered correctly in the matrix
  fp=len(keys_false_positives)
  fn=len(keys_false_negatives)
  acc = float(tp+tn)/(tp+tn+fp+fn)
  prec=float(tp)/(tp+fp)
  rec=float(tp)/(tp+fn)
  tpr=float(tp)/(tp+fn)
  tnr=float(tn)/(tn+fp)
  bal_acc=0.5*(tpr+tnr)
  fpr=float(fp)/(tn+fp)
  infd=tpr+tnr-1.0
  if not file_handle is None:
     file_handle.write('TP/TN/FP/FN %d/%d/%d/%d\n' %(tp,tn,fp,fn))
     file_handle.write('Acc/Prec/tpr/tnr/fpr %f/%f/%f/%f/%f\n'%(acc,prec,tpr,tnr,fpr))
     file_handle.write('bal_acc/infd %f/%f\n' %(bal_acc,infd))

def plot_heatmaps_from_pickled_dataframe(pickledfiles='PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_UpdateDist_plusone_*.p',basename='results/heatmaps_rmsd_filtered/',N=10,rmsd_cutoff=1.0,criteria='rmsd_less_than',prop = 'dist_498OG-288O',cmap=None,off=0.10,plotdist=False,fitstuff=False,reverse=False,maxmin=None,annotate=False,suffix=None,applyfraction=False,show=True,sc_cutoff=None,rocs=False,heatmaps=True,off_list=[0.0,0.05,0.10,0.20,0.30,0.40,0.50,0.55],positionwise=True,pointplot=False,by='median',pos=['+1'],listres1=aalist,boxplot=True,pick_fpr=None, pick_tpr=0.92,pick_threshold=None):
  from functions_matrices_from_clusters import get_roc_auc_for_dict
  if cmap is None:
    cmap = 'Greens'
  listofdfs =[]
  listofdfs_nofilter = []
  pfiles = glob.glob(pickledfiles)
  dict_stats = {}
  dict_stats_size = {}
  pfiles.sort()
  print(len(pfiles))
  for pfile in pfiles:
    #print(pfile)
    df, df_energies_res, df_energies_res_res = read_pickled_file(pfile)
    key = df['key'][0]
    exp = df['exp2018'][0]
    total_size = float(len(df))
    df_filter = get_filtered_dataframe(df,criteria, rmsd_cutoff)
    if df_filter is None:
      continue
    
    #print(Nsize,len(df_filter),total_size)
    Nsize = len(df_filter)/float(total_size)
    print(Nsize,len(df_filter),total_size)
    if applyfraction:
      if Nsize < 0.10:
        continue
    if not reverse:
      df_filter_top = df_filter[:N]
    else:
      df_top = df[:N]
      df_filter_top = get_filtered_dataframe(df,criteria, rmsd_cutoff)
      
      if df_filter_top is None:
        continue
    if not sc_cutoff is None:
      if df_filter_top['sc_shapecomplementarity'].mean() < sc_cutoff:
        continue
    dict_stats[key] ={}
    get_stats(df_filter_top,prop,key,dict_stats)


    dict_stats_size[key]=Nsize
    dict_stats[key]['exp2018'] = exp
    listofdfs.append(df_filter_top)
    listofdfs_nofilter.append(df)
  df_concat = pd.concat(listofdfs,ignore_index=True)    
  df_concat_nofilter = pd.concat(listofdfs_nofilter,ignore_index=True)
  tempdata_size = np.array(list(dict_stats_size.values()))
  prim_key = dict_stats.keys()
  key = list(prim_key)[0]
  sec_keys =  dict_stats[key].keys()
  if pointplot:
    get_positionwise_plots_2(df_concat,prop,criteria,rmsd_cutoff,N,basename,show=show,suffix=suffix,by=by,maxmin=maxmin,pos=pos)
  if rocs:
    suffix='tpr%f' %pick_tpr
    if not pick_threshold is None:
       suffix='threshold%f' %pick_threshold
    if not pick_fpr is None:
       suffix='threshold%f' %pick_fpr
    filename = 'results_aucscores/AUCscores_%s_%s_%f_top%d_%s.txt' %(prop,criteria,rmsd_cutoff,N,suffix)
    file_handle_auc = open(filename,'w')
    filename = 'results_aucscores/AccuracyMetrics_%s_%f_%s_top%d_%s.txt' %(criteria,rmsd_cutoff,prop,N,suffix)
    file_handle = open(filename,'w')
  for skey in sec_keys:
   if skey == 'exp2018': continue
   if positionwise:
    get_positionwise_plots(dict_stats,skey,prop,criteria,rmsd_cutoff,N,basename,show=show,suffix=suffix,by=by,boxplot=boxplot,pointplot=pointplot)
   #if get_pdbs:
   # get_rep_pdb(df_concat,dict_stats,skey,prop,criteria,rmsd_cutoff,N,basename,show=show)
   if heatmaps:
    print()
    #setup_and_plot_heatmap_from_dict(dict_stats,skey,prop,criteria,rmsd_cutoff,N,basename,show=show,suffix=suffix,listres1=listres1) 
   if rocs:
    dicttemp ,_,_ = get_dict_from_doubledict(dict_stats,skey)
    value_name = prop
    if prop=='rmsd':
      value_name = 'substrate_ca_no_super_rmsd'
    if skey=='sd':
      value_name = prop +'_%s' %skey
    from functions_matrices_from_clusters import get_false_positives_and_false_negatives,heatmap_topN_sequon, heatmap
    for off in off_list:
      #_,_,_,score = get_roc_auc_for_dict(dicttemp,value_name,off=off)
      #file_handle.write('%f,%s,%s,%f\n'%(off,prop,skey,score)) 
      prefix= 'results_aucscores/FullAccuracyMetrics_%s_%f_%s_%s_off%f_top%d' %(criteria,rmsd_cutoff,prop,skey,off,N)
      suffix='tpr%f' %pick_tpr
      if not pick_threshold is None:
       suffix='threshold%f' %pick_threshold
      if not pick_fpr is None:
       suffix='threshold%f' %pick_fpr

      outfile=prefix+'_'+suffix+'.out'
      file_handle.write('%s %s off %f\n' %(prop,skey,off))
      file_handle_auc.write('%s %s off %f\n' %(prop,skey,off))
      print(outfile)
      if len(list(dicttemp.values()))<10: continue
      dictmat_base,dictmat_compare,list_fps,list_fns,list_tps,list_tns = get_false_positives_and_false_negatives(dicttemp,value_name,off=off,pick_tpr=pick_tpr,outfile=outfile,pick_threshold=pick_threshold,pick_fpr=pick_fpr,file_handle_auc=file_handle_auc)
      #write_metrics(list_tps,list_tns,list_fps,list_fns,file_handle=file_handle)

      suffix = '%s_' %prop + '_%s' %skey + '_off%f_pick_tpr%f' %(off,pick_tpr)
      if not pick_fpr is None:
        suffix = '%s_' %prop + '_%s' %skey + '_off%f_fpr%f' %(off,pick_fpr)
      if not pick_threshold is None:
        suffix = '%s_' %prop + '_%s' %skey + '_off%f_threshold%f' %(off,pick_threshold)
      #heatmap_topN_sequon(dictmat_compare , heatmap_type=heatmap.compare, criteria=criteria , filterTopN=True ,suffix=suffix,annotate=False)
      
  if rocs:
    file_handle.close()
    file_handle_auc.close()

if __name__=='__main__':
  #for aa in ['S','T','P','A','G']:
    #for by in ['min']:
    #plot_heatmaps_from_pickled_dataframe(pickledfiles='PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_UpdateDist_plusone_A%sT*APRC.p' %aa,prop = 'sc_shapecomplementarity',N=5,rmsd_cutoff=1.0,criteria='rmsd_less_than',reverse=False,rocs=False,heatmaps=False,positionwise=False,show=True,suffix=aa,pointplot=True,by='median',maxmin=[0.65,0.80])
  for prop in ['psi_498','psi_499','phi_499','psi_500','phi_500']:
      plot_heatmaps_from_pickled_dataframe(pickledfiles='PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_UpdateDist_plusone_A[AVTPGS]T*APRC.p' ,prop = prop,N=5,rmsd_cutoff=1.0,criteria='rmsd_less_than',reverse=False,rocs=False,heatmaps=False,positionwise=False,show=True,suffix=aa,pointplot=True,by='sd',pos='-1')

  #plot_heatmaps_from_pickled_dataframe(pickledfiles='PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_UpdateDist_plusone_AST*APRC.p',prop = 'dist_500CA-287CB',N=5,rmsd_cutoff=1.0,criteria='rmsd_less_than',reverse=False,rocs=False,heatmaps=False,positionwise=False,show=False,suffix='S',pointplot=True)
