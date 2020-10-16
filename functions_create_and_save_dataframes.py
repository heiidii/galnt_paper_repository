
import pickle
from functions_matrices_from_clusters_criteria import *
from functions_matrices_from_clusters import *
import pandas as pd

def apply_criterion(criteria,heatmap_type,field_name):

  if heatmap_type == heatmap.value:
    dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name])

  if heatmap_type == heatmap.size:
    dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type)

  return dictmat

def generatelist():
  from functions_matrices_from_clusters import get_files
  from DataSortingFiltering import get_clusters_from_combined_dataframe
  pattern = 'PackandRelax_498*/unglycosylated/score.sc'
  files = get_files(pattern )

  dictmat={}
  for curfile in files:
    print(curfile)
    suffix = curfile.split('/')[0]
    df, centroids, sizedict, sizedict_normalized,labels,n_clusters,core_samples_ =  get_clusters_from_combined_dataframe([curfile],filterTopN=True)
    outf = open('list_%s' %suffix, 'w')
    for curi in df.index.values:
      print(curi,df.loc[curi,'description'])
      outf.write(df.loc[curi,'description']+'\n')
    outf.close()

def generatelist_from_dfs(pattern='',topN=1,prefix='../',filter_by_rmsd=False,suffix=''):
  files = glob.glob(pattern)
  outf = open('list_pdbfiles_topN_%d%s' %(topN,suffix), 'w')
  outf_v = open('list_pdbfiles_topN_%d_verify%s' %(topN,suffix), 'w')
  for pfile in files:
    if not os.path.exists(pfile): return
    data_ = pickle.load(open(pfile,'rb'))
    df = data_['dataframe']
    key = df['key'][0]
    if filter_by_rmsd:
      df_filter = df [ df['rmsd']<1.0 ]
      del df
      df = df_filter
    for curi in (list(df.index.values))[:topN]:
      outf.write(prefix+df.loc[curi,'description']+'\n')
      outf_v.write('%d %s %5.2f\n' %(curi,df.loc[curi,'description'],df.loc[curi,'IE']))
  outf.close()
  outf_v.close()

def add_prop_columns_to_df(df,prop='distance'):
    from functions_AAprops import getdistancemetricforfile, get_sc_for_file, getdihedrals_for_file,get_fnat_and_fnotnat,getdistanceforfile_OHsidechain_Obackbone_bond,getdistanceforfile_plusone
    dict_array = {}
    for curi in df.index.values:
      pdbfile = df.loc[curi,'description']
      print(pdbfile)
      dict_ ={}
      if prop == 'sc':
        dict_ = get_sc_for_file(pdbfile)
      if prop=='distance':
        dict_ = getdistancemetricforfile(pdbfile)
      if prop=='dih':
        dict_ = getdihedrals_for_file(pdbfile)
      if prop=='fnat':
        dict_ = get_fnat_and_fnotnat(pdbfile)
      if prop=='distance_thr_ser_bk':
        dict_ = getdistanceforfile_OHsidechain_Obackbone_bond(pdbfile)
      if prop=='distance_plusone_F287':
        dict_ = getdistanceforfile_plusone(pdbfile)
      
      for key in dict_:
        if not key in dict_array:
          dict_array[key] = []
        dict_array[key].append(dict_[key])

    for key in dict_array:
      df[key] = pd.Series(dict_array[key])
    return df


def get_energy_df_for_df(df,residues=[498,499,500],etype = 'res-res'):
    '''
    This is when for each pdbfile you want to generate a new dataframe
    and add that to pickle file
    '''
    from functions_resenergies import PairwiseEnergyMapForPdbFile, EnergyMapForPdbFile
    dict_array = {}
    for curi in df.index.values:
      pdbfile = df.loc[curi,'description']
      sequonkey = df.loc[curi,'key']
      #print(pdbfile)
      dict_ ={}
      #NOTE: these functions return dict of array so when combining then each value must be appended
      if etype=='res-res':
        dict_ = PairwiseEnergyMapForPdbFile(pdbfile,tag=sequonkey,residues = residues)  
      if etype=='res-total':
        dict_ = EnergyMapForPdbFile(pdbfile,tag=sequonkey,residues=residues)
      for key in dict_:
        if not key in dict_array:
          dict_array[key] = []
        for value in dict_[key]:
          dict_array[key].append(value)

    new_df = pd.DataFrame()

    for key in dict_array:
      new_df[key] = pd.Series(dict_array[key])
    return new_df


def add_distance_columns_to_df(df):
    from functions_AAprops import getdistancemetricforfile
    dict_distances_array = {}
    for curi in df.index.values:
      print(curi,df.loc[curi,'description'])
      pdbfile = df.loc[curi,'description']
      for key in dict_distances:
        if not key in dict_distances_array:
          dict_distances_array[key] = []
        dict_distances_array[key].append(dict_distances[key])
    for key in dict_distances_array:
      df[key] = pd.Series(dict_distances_array[key])
    return df

def add_distance_columns_to_df_for_files(files=None):
  from functions_matrices_from_clusters import get_files
  from DataSortingFiltering import get_clusters_from_combined_dataframe
  if files is None:
    pattern = 'PackandRelax_498*/unglycosylated/score.sc'
    files = get_files(pattern )
  for curfile in files:
    print(curfile)
    suffix = curfile.split('/')[0]
    df, centroids, sizedict, sizedict_normalized,labels,n_clusters,core_samples_ =  get_clusters_from_combined_dataframe([curfile],filterTopN=True)
    df = add_distance_columns_to_df(df)
    print(df)

def generate_distance_df_for_files(files,reslist1,reslist2,tags):
  from functions_AAprops import getdistancemetricforfile
  dictdistarray = {}
  dictdistarray['description']=[]
  dictdistarray['tag']=[]
  dictdistarray['distance']=[]
  dictdistarray['distance_type']=[]
  dictdistarray['distance_type_pep']=[]
  dictdistarray['distance_type_pepAtom']=[]
  dictdistarray['distance_type_udpAtom']=[]
  for j,fname in enumerate(files):
    dictdist,tupledict = getdistancemetricforfile(fname,reslist1[j],reslist2[j])
    for key in dictdist:
      dictdistarray['distance'].append(dictdist[key])
      dictdistarray['distance_type'].append(key)
      dictdistarray['description'].append(fname)
      dictdistarray['tag'].append(tags[j])
      dictdistarray['distance_type_pep'].append(tupledict[key][0])
      dictdistarray['distance_type_pepAtom'].append(tupledict[key][1])
      dictdistarray['distance_type_udpAtom'].append(tupledict[key][2])
  df = pd.DataFrame.from_dict(dictdistarray)
  print(df)
  return df


def serialize_processed_data_nofilter():
  from functions_matrices_from_clusters import get_files
  from DataSortingFiltering import get_filtered_combined_dataframe
  pattern = 'PackandRelax_498*/unglycosylated/score.sc'
  files = get_files(pattern )

  for curfile in files:
    print(curfile)
    df_nofilter  = get_filtered_combined_dataframe([curfile])
    key = df_nofilter['key'][0]
    pickle.dump(df_nofilter,open('PickleFiles/20191124/UnglycosylatedPeptideData_NoFilter_%s.p' %key,'wb'))

def serialize_processed_data_for_file(curfile,clus2d=False,addDistance=True,addSc=True,write=True,verbose=False,path='',filename='UnglycosylatedPeptideData_Filtered_ClusteredDBScan'):

    from DataSortingFiltering import get_clusters_from_combined_dataframe
    from functions_matrices_from_clusters import get_sinks_for_columns
    from HelperPlotting import get_matrix_key_from_file
    #print(curfile)
    key = get_matrix_key_from_file(curfile)
    if write:
      if not os.path.exists(path):
        os.system(path)
      if clus2d:
        outfilename = '%s/%s2D_%s.p' %(path,filename,key)
      else:
        outfilename = '%s/%s_%s.p' %(path,filename,key)
      if os.path.exists(outfilename):
        #print('file %s exists' %outfilename)
        return
      else:
        print(curfile)
        print('file %s does not exist' %outfilename)
        #return
    df, centroids, sizedict, sizedict_normalized,labels,n_clusters,core_samples_ =  get_clusters_from_combined_dataframe([curfile],filterTopN=True,clus2d=clus2d)
    #if addDistance:
    #  df = add_prop_columns_to_df(df,prop='distance')
    if addSc:
      df = add_prop_columns_to_df(df,prop='sc') 
    
    sorted_d = dict([(k, sizedict_normalized[k]) for k in sorted(sizedict_normalized, key=sizedict_normalized.get, reverse=True)])
    sinks = get_sinks_for_columns(df,labels,type_value='sinks')
    sizedict_sorted={}
    sizedict_normalized_sorted={}
    centroids_sorted=centroids
    sinks_sorted = sinks
    list_sorted_clusterids=[]
    i=0
    for curkey in sorted_d:
      if curkey == -1: continue
      centroids_sorted[i,:] = centroids[curkey,:]
      sizedict_sorted[i] = sizedict[curkey]
      sizedict_normalized_sorted[i] = sizedict_normalized[curkey]
      sinks_sorted[i,:] = sinks[curkey,:]
      list_sorted_clusterids.append(curkey)
      i+=1

    labels_sorted =[]
    list_sorted_labels = list(sorted_d.keys())
    if verbose:
      print('list_sorted_labels ',list_sorted_labels)
      print('list corrected sorted labels ',list_sorted_clusterids)
    # [1,0,-1]
    for label in labels:
      if label != -1:
        labels_sorted.append(list_sorted_clusterids.index(label))
      else:
        labels_sorted.append(-1)
    #add dummy row for each label
    for label in [0,1,2]:
      if label in list_sorted_labels: continue
      dummy_list = [ np.nan for _ in df.columns.values ]
      dummy_row = pd.Series(dummy_list, index=df.columns.values)

    #key = df['key'][0]
    assert(key == df['key'][0])

    data ={}
    data['dataframe'] = df
    data['source'] = curfile
    data['centroids_sorted'] = centroids_sorted
    data['sizedict_sorted'] = sizedict_sorted
    data['sizedict_normalized_sorted']=sizedict_normalized_sorted
    data['labels_sorted']=labels_sorted
    data['n_clusters'] = n_clusters
    #data['core_samples_sorted']=core_samples_
    data['sinks_sorted']=sinks_sorted
    if write:
      pickle.dump(data,open(outfilename,'wb'))
      print(outfilename)

def serialize_processed_data(pattern='PackandRelax_498*/unglycosylated/score.sc',clus2d=False,addDistance=False,addSc=False,write=True,verbose=False,path=''):
  from functions_matrices_from_clusters import get_files
  files = get_files(pattern )
  for curfile in files:
    serialize_processed_data_for_file(curfile,clus2d=clus2d,addDistance=addDistance,addSc=addSc,write=write,verbose=verbose,path=path)



def update_stored_dataframe_distances(pfile,outfile=None):
    if not os.path.exists(pfile): return
    data_ = pickle.load(open(pfile,'rb'))
    df = data_['dataframe']
    key = df['key'][0]
    df = add_distance_columns_to_df(df)
    data_['dataframe'] = df
    if outfile is None:
      pickle.dump(data_,open('PickleFiles/20191203/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_UpdateDistance%s_.p' %key,'wb'))
  

def update_stored_dataframe_shapecomplementarity(pfile,outfile=None):
    if not os.path.exists(pfile): return
    data_ = pickle.load(open(pfile,'rb'))
    df = data_['dataframe']
    key = df['key'][0]
    df = add_prop_columns_to_df(df,prop='sc')
    data_['dataframe'] = df
    if outfile is None:
      pickle.dump(data_,open('PickleFiles/20191210/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_UpdateSC_%s.p' %key,'wb'))

def update_stored_dataframe_dihedrals(pfile,outfile=None):
    if not os.path.exists(pfile): return
    print(pfile)
    data_ = pickle.load(open(pfile,'rb'))
    df = data_['dataframe']
    key = df['key'][0]
    df = add_prop_columns_to_df(df,prop='dih')
    data_['dataframe'] = df
    if outfile is None:
      pickle.dump(data_,open('PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_UpdateDih_%s.p' %key,'wb'))

def update_stored_dataframe_fnat(pfile,outfile=None):
    if not os.path.exists(pfile): return
    print(pfile)
    data_ = pickle.load(open(pfile,'rb'))
    df = data_['dataframe']
    key = df['key'][0]
    df = add_prop_columns_to_df(df,prop='fnat')
    data_['dataframe'] = df
    if outfile is None:
      pickle.dump(data_,open('PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_UpdateFnat_%s.p' %key,'wb'))

def update_stored_dataframe_distance_thr_bk(pfile,outfile=None):
    if not os.path.exists(pfile): return
    print(pfile)
    data_ = pickle.load(open(pfile,'rb'))
    df = data_['dataframe']
    key = df['key'][0]
    df = add_prop_columns_to_df(df,prop='distance_thr_ser_bk')
    data_['dataframe'] = df
    if outfile is None:
      pickle.dump(data_,open('PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_UpdateDist_Thr_Ser_Bk_%s.p' %key,'wb'))

def update_stored_dataframe_distance_plusone_F287(pfile,outfile=None):
    if not os.path.exists(pfile): return
    print(pfile)
    data_ = pickle.load(open(pfile,'rb'))
    df = data_['dataframe']
    key = df['key'][0]
    df = add_prop_columns_to_df(df,prop='distance_plusone_F287')
    data_['dataframe'] = df
    if outfile is None:
      pickle.dump(data_,open('PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_UpdateDist_plusone_%s.p' %key,'wb'))

def update_stored_dataframe_pairwise_energies(pfile,outfile=None):
    if not os.path.exists(pfile): return
    print(pfile)
    data_ = pickle.load(open(pfile,'rb'))
    df = data_['dataframe']
    key = df['key'][0]
    new_df = get_energy_df_for_df(df[:50])
    data_['dataframe_pairwise_energies'] = new_df
    new_df_2 = get_energy_df_for_df(df[:50],etype='res-total')
    data_['dataframe_total_res_energies'] = new_df_2
    if outfile is None:
      outfile = 'PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_Energies_%s.p' %key
      if not os.path.exists(outfile):
        pickle.dump(data_,open(outfile,'wb'))

def generate_and_store_processed_dataframe():
    from run_matrices_from_cluster import get_deltas
    import compare_simulations_to_experiments as CompSE
    import pandas as pd
    import numpy as np
    dictmat_exp = CompSE.getexperimentaldatafrompickle_dict()
    df_PTX = pd.DataFrame()
    df_PTX['experimental'] = pd.Series(dictmat_exp)

    #step1 - Interaction PTX
    stepwise_dictmat = {}
    criteria =  hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
    heatmap_type = heatmap.value
    field_name = 'interaction_energy'
    stepwise_dictmat[field_name] = apply_criterion(criteria,heatmap_type,field_name)
    dictmat_PTX = stepwise_dictmat[field_name]

    df_PTX['interaction_energy_PTX']=pd.Series(dictmat_PTX)

    print("0. ",df_PTX)
    #step2 - Cluster size PTX
    criteria =  hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
    heatmap_type = heatmap.size
    field_name = 'cluster_size'
    stepwise_dictmat[field_name] = apply_criterion(criteria,heatmap_type,field_name)
    dictmat_size = stepwise_dictmat[field_name]
    df_PTX['cluster_size']=pd.Series(dictmat_size)
    df_PTX['cluster_size'] = df_PTX['cluster_size'].fillna(0.0)
    print("1. ",df_PTX)

    criteria =  hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
    heatmap_type = heatmap.value
    field_name = 'substrate_ca_no_super_rmsd'
    stepwise_dictmat[field_name] = apply_criterion(criteria,heatmap_type,field_name)
    df_PTX['substrate_ca_no_super_rmsd']=pd.Series(stepwise_dictmat[field_name])

    #step3 - GTX
    criteria =  hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_rmsd
    heatmap_type = heatmap.value
    field_name = 'interaction_energy'
    dictmat_GTX = apply_criterion(criteria,heatmap_type,field_name)
    df_PTX['interaction_energy_GTX']=pd.Series(dictmat_GTX)

    dictmat_delta_PTX_GTX = get_deltas(dictmat_PTX,dictmat_GTX)
    default_value = -2
    heatmap_type = heatmap.size
    field_name = 'cluster_size'
    stepwise_dictmat[field_name] = apply_criterion(criteria,heatmap_type,field_name)
    dictmat_size_GTX = stepwise_dictmat[field_name]
    df_PTX['cluster_size_GTX']=pd.Series(dictmat_size_GTX)
    df_PTX['cluster_size_GTX'] = df_PTX['cluster_size_GTX'].fillna(0.0)
    df_PTX['interaction_energy_delta_GTX']=pd.Series(dictmat_delta_PTX_GTX)
    df_PTX['interaction_energy_delta_GTX'] = df_PTX['interaction_energy_delta_GTX'].fillna(-2.0)


    #step4 - TTQ
    criteria =  hb_bond_criteria.significant_clusters_sinks_TTY_cutoff
    heatmap_type = heatmap.value
    field_name = 'interaction_energy'
    dictmat_TTQ = apply_criterion(criteria,heatmap_type,field_name)
    df_PTX['interaction_energy_TTQ']=pd.Series(dictmat_TTQ)

    dictmat_delta_PTX_TTQ = get_deltas(dictmat_PTX,dictmat_TTQ)
    field_name = 'cluster_size'
    heatmap_type = heatmap.size
    dictmat_size_TTQ = apply_criterion(criteria,heatmap_type,field_name)
    df_PTX['cluster_size_TTQ'] = pd.Series(dictmat_size_TTQ)
    df_PTX['cluster_size_TTQ'] = df_PTX['cluster_size_TTQ'].fillna(0.0)

    df_PTX['interaction_energy_delta_TTQ'] = pd.Series(dictmat_delta_PTX_TTQ)
    df_PTX['interaction_energy_delta_TTQ'] = df_PTX['interaction_energy_delta_TTQ'].fillna(-2.0)

    criteria =  hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
    heatmap_type = heatmap.value
    field_name = 'substrate_sequon_ca_no_super_rmsd'
    stepwise_dictmat[field_name] = apply_criterion(criteria,heatmap_type,field_name)
    df_PTX['substrate_sequon_ca_no_super_rmsd']=pd.Series(stepwise_dictmat[field_name])
    print('2. ',df_PTX)
    df_PTX.to_pickle("PickleFiles/20191124/clustered_top10percent_PTX_GTX_TTQ_rmsd_dataframe.p")
    dictcombined = df_PTX.to_dict('dict')
    
    pickle.dump(dictcombined,open("PickleFiles/20191124/clustered_top10percent_PTX_GTX_TTQ_rmsd_dict.p",'wb'))


def get_conserved_distances_df():
  files =['PackandRelax_T12_YY_correct_531THR_533PRO/T12_YYITP_unglycosylated_correct_YY_1_531THR_533PRO.pk.pdb',
          'PackandRelax_498THR_500PRO/T2_with_UDPGalNAc_withS_from_4d0z_ATTAAPRC_498THR_500PRO.pk.pdb',
          'T4_6h0b_chainBandcatalyticF_HWUandMnfromT12_0001.pdb']
  pepres = [[531,532,533],[498,499,500],[525,526,527]]
  udpres = [527,495,1]
  tags = ['T12','T2','T4']
  df = generate_distance_df_for_files(files,pepres,udpres,tags)
  pickle.dump(df,open("PickleFiles/ConservedDistances/conserveddistances_dataframe.p",'wb'))
  import matplotlib.pyplot as plt
  import seaborn as sns
  df_sorted_1 = df.sort_values(['distance_type_pepAtom','distance_type_udpAtom','distance_type_pep','distance'])
  ax1 = sns.pointplot(data=df_sorted_1,x='distance_type',y='distance',hue='tag',palette='Set2')
  ax1.tick_params(axis='x',labelsize=9, grid_linewidth=0.5, grid_linestyle='-', grid_color='black',grid_alpha=0.6,labelrotation=70)
  ax1.tick_params(axis='y', labelcolor='gray',labelsize=15, grid_linewidth=0.5, grid_linestyle='-', grid_color='black',grid_alpha=0.6)
  ax1.xaxis.set_label_text("")
  ax1.yaxis.set_label_text("")
  ax1.grid(linestyle='-', linewidth='0.5', color='gray')
  plt.tight_layout()
  #plt.savefig('PickleFiles/ConservedDistances/conserveddistances_all.png',transparent=True,dpi=1200)
  #plt.show()
  plt.close()
  df_filter = df[ df['distance']<=15.0]
  df_sorted = df_filter.sort_values(['distance_type_pep','distance_type_pepAtom','distance_type_udpAtom'])
  #df_filter
  fig, ax1 = plt.subplots(figsize=(15,8))
  ax1 = sns.pointplot(data=df_sorted,x='distance_type',y='distance',hue='tag',palette='Set2',ax=ax1)
  ax1.tick_params(axis='x',labelsize=10, grid_linewidth=0.5, grid_linestyle='-', grid_color='black',grid_alpha=0.6,labelrotation=70)
  ax1.tick_params(axis='y', labelcolor='gray',labelsize=15, grid_linewidth=0.5, grid_linestyle='-', grid_color='black',grid_alpha=0.6)
  ax1.xaxis.set_label_text("")
  ax1.yaxis.set_label_text("")
  ax1.grid(linestyle='-', linewidth='0.5', color='gray')
  plt.tight_layout()
  plt.savefig('PickleFiles/ConservedDistances/conserveddistances.png',transparent=True,dpi=1200)
  plt.show()
  return df

def parallelrun(pattern,runfunc='dihedrals',n_pool=8):
  from functions_matrices_from_clusters import get_files
  files = get_files(pattern ) 
  import multiprocessing as mp
  if pattern.find('.sc') == -1:
    files = glob.glob(pattern)
  print(len(files))
  pool = mp.Pool(min(n_pool,len(files)))
  #pool.map(serialize_processed_data_for_file,files)
  if runfunc=='dihedrals': 
    pool.map(update_stored_dataframe_dihedrals,files)
  if runfunc=='energies':
    pool.map(update_stored_dataframe_pairwise_energies,files)
  if runfunc=='fnat':
    pool.map(update_stored_dataframe_fnat,files)
  if runfunc=='distance_thr_ser_bk':
    pool.map(update_stored_dataframe_distance_thr_bk,files)
  if runfunc=='distance_plusone_F287':
    pool.map(update_stored_dataframe_distance_plusone_F287,files)

if __name__ == '__main__':    
  #generate_and_store_processed_dataframe()
  #serialize_processed_data(pattern='PackandRelax_498PHE_500GLN/unglycosylated/score.sc',clus2d=False,addDistance=True,addSc=True,path='PickleFiles/20200207/')
  serialize_processed_data(pattern='PackandRelax_498*/unglycosylated/score.sc',clus2d=False,addDistance=True,addSc=True,path='PickleFiles/20200207/')
  #pattern = 'PackandRelax_498*/unglycosylated/score.sc'
  #pattern = 'PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_UpdateFnat_A*T*APRC.p'
  #parallelrun(pattern,runfunc='distance_plusone_F287',n_pool=4)
  #generatelist_from_dfs(pattern=pattern,topN=1,prefix='../',filter_by_rmsd=True,suffix='filtered_rmsd')
  #get_conserved_distances_df() 
