import os
import sys
import pandas as pd
import numpy as np
from MultivariablePlottingGlobal import info, labels_df, fields, categories
import compare_simulations_to_experiments as CompSE
from HelperPlotting import *


def getdatafromfile(curfile,use_fields=None):
    f = open(curfile,'r')
    lines = f.readlines()
    f.close()

    data = dict()
    from functions_filterpdbforparams import getfiltereddatafromlines
    #print('getdatafromfile ',fields)
    if use_fields is None:
      fields_local = fields
      data, fields_temp = getfiltereddatafromlines(lines,fields_local)
    if use_fields=='all':
      data, fields_temp = getfiltereddatafromlines(lines)
    #print('getdatafromfile ',data)
    return data

def getdata_xyzc(inputfile,serialize=True,json=True):

    def get_full_pdbname(curfile,description):
      fullname = curfile.split('/')[0]+'/'+description+'.pdb.gz'
      #print(description,fullname)
      return fullname

    dictdata = dict()
    for curfield in ['xfield','yfield','zfield','cfield','z2field','description','glyc_yfield','glyc_z2field','z3field']:
      dictdata[curfield] = dict()
    data = getdatafromfile(inputfile)
    if data is None: return None

    for curfield in ['xfield','yfield','zfield','cfield','z2field','description','glyc_yfield','glyc_z2field','z3field']:
      if info[curfield] in data:
        dictdata[curfield]["vals"] = data[info[curfield]]

    #fix pdbfile names    
    if 'description' in dictdata:
      for i,val in enumerate(dictdata['description']['vals']):
        dictdata['description']['vals'][i] = get_full_pdbname(inputfile,val)

    return dictdata

def get_df_from_dict(dictdata,entries=None):
    cleandict = {}
    if entries==None:
      entries = ['xfield','yfield','zfield','cfield','glyc_yfield','z3field','description']
    for curfield in entries:
      if curfield in dictdata:
        #print(curfield)
        if len(dictdata[curfield]['vals'])<1: continue
        cleandict[labels_df[info[curfield]]] = dictdata[curfield]['vals']

    df = pd.DataFrame(cleandict)
    return df

def get_dfs_from_dicts(listofdicts,keys):
  listofdfs = []
  for inpdict,key in zip(listofdicts,keys):
    df = get_df_from_dict(inpdict)
    listofdfs.append(df)
  return get_dfs_from_dicts


def get_combined_df(listofdicts,keys):
  listofdfs = []
  for inpdict,key in zip(listofdicts,keys):
    df = get_df_from_dict(inpdict) 
    Ndf = df.shape[0]
    #print(df.shape , Ndf)
    df['sequon']=[key for _ in range(0,Ndf)]
    matrix_key = get_matrix_key(key)
    df['key']=[matrix_key for _ in range(0,Ndf)]
    key1, key2 = get_keys_from_key(key)

    for cat in categories:
      catlist = categories[cat]
      if key2 in catlist:
        df['category']=[cat for _ in range(0,Ndf)]
        break

    exp = CompSE.get_experimental_data_for_sequon(key1,key2)
    df['exp2018']=[exp for _ in range(0,Ndf)]
    val = 'MEDIUM'
    if exp <= 0.10:
      val = 'OFF'
    elif exp > 0.70:
      val = 'HIGH'
      
    df['exp2018_category'] = [val for _ in range(0,Ndf)]

    #print(df.shape , Ndf)
    listofdfs.append(df)
  return pd.concat(listofdfs,ignore_index=True)


def filter_by_expdata(df,expdata):
    if expdata=='on':
      df_filter = df[df['exp2018']>0]
    elif expdata=='off':
      df_filter = df[df['exp2018']<0.10]
    elif expdata=='high':
      df_filter = df[df['exp2018']>=0.70]
    elif expdata=='medium':
      temp_ = df[ (df['exp2018']<0.70)]
      df_filter = temp_[temp_['exp2018']>0.10 ]
      del temp_
    del df
    return df_filter

def filter_interaction_energy(df,filterbyfield):
      if filterbyfield[key]=='high':
        #print("filtering ",key,filterbyfield[key])
        df_filter =        df[ ( df [labels_df[key] ] <= -34 ) ]
      elif filterbyfield[key]=='medium':
        temp_ =        df[ ( df [labels_df[key] ] <= -15 ) ]
        df_filter =        temp_[ ( temp_ [labels_df[key] ] > -34 ) ]
        del temp_
      elif filterbyfield[key]=='low':
        temp_ =        df[ ( df [labels_df[key] ] <= -5 ) ]
        df_filter =        temp_[ ( temp_ [labels_df[key] ] > -15 ) ]
        del temp_
      del df
      df = df_filter
      return df

def filter_by_field(df,filterbyfield):
  for key in filterbyfield:
      if key=='interaction_energy':
        return filter_interaction_energy(df,filterbyfield)
        

def get_filtered_combined_dataframe(files,category='all',expdata='all',filterbyfield=None,filterTopN=None):
  listofdicts = []
  keys=[]
  for curfile in files:
    listofdicts.append(getdata_xyzc(curfile))
    keys.append(get_key(curfile))
  #print("KEYS",keys)
  #print('Dicts',listofdicts)
  combined_df = get_combined_df(listofdicts,keys)


  combined_df_filter = None
  if not category=='all':
    #filter
    #print('filtering by',category)
    combined_df_filter = combined_df[combined_df['category']==category]
    del combined_df
    combined_df = combined_df_filter

  if not expdata=='all':
    #print('filtering by exp data ',expdata)
    combined_df_filter = filter_by_expdata(combined_df,expdata)
    combined_df = combined_df_filter

  if not filterbyfield is None:
    combined_df_filter = filter_by_field(combined_df,filterbyfield)
    combined_df = combined_df_filter

  if not filterTopN is None:
    #print("Sorting by interaction energies")
    combined_df.sort_values(by=labels_df['interaction_energy'],inplace=True)
    combined_df.reset_index(drop=True,inplace=True)
    #print(combined_df[labels_df['interaction_energy']][:20])

  return combined_df

def get_clusters_from_combined_dataframe(files,category='all',expdata='all',filterbyfield=None,k_means=False,filterTopN=None,TopN=200,clus2d=False):

  combined_df = get_filtered_combined_dataframe(files,category=category,expdata=expdata,filterbyfield=filterbyfield,filterTopN=filterTopN)
  if combined_df is None: return
  if combined_df.shape[0]<5: return
  if not filterTopN is None:
    #print("Taking TopN values",combined_df.shape)
    combined_df_TopN = combined_df[:TopN]
    del combined_df
    combined_df = combined_df_TopN
    #print("Taking TopN values",combined_df.shape)
    #print("sequon rmsd ", combined_df[labels_df[info['z3field']]].values) 
  if clus2d:
    listfromdf = np.array([ combined_df[labels_df[info['xfield']]], combined_df[labels_df[info['yfield']]] ])
  else:
    listfromdf = np.array([ combined_df[labels_df[info['xfield']]] , combined_df[labels_df[info['yfield']]], combined_df[labels_df[info['zfield']]] ])
  from ClusterPoints import get_clusterpoints_from_list
  centroids , sizedict, sizedict_normalized,labels,n_clusters,core_samples_ = get_clusterpoints_from_list(listfromdf,k_means=k_means,core_samples=True)
  return combined_df,centroids , sizedict, sizedict_normalized,labels,n_clusters,core_samples_
