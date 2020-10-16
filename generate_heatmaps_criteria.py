from functions_matrices_from_clusters import *
import copy
off=0.10

def get_deltas(dict1,dict2):
      deltadict = {}
      keys_all = list(dict1.keys()) + list(dict2.keys())
      unique_keys = set(keys_all)
      print(unique_keys)
      for sequon in unique_keys:
        value1 = None
        value2 = None
        delta = None
        if sequon in dict1:
          value1 = dict1[sequon]
        if sequon in dict2:
          value2 = dict2[sequon]

        if (not value1 is None) and (not value2 is None):
          delta = value1 - value2
        if not delta is None:
          deltadict[sequon] = delta
      return deltadict


def get_heatmaps_for_criteria_for_all_clusters(criteria,show=True,roc=False,suffix=None,clus2d=False,annotate=False,pattern='PackandRelax_498*/unglycosylated/score.sc'):

  list_of_clusters = [0,1] #only matters for signicant clusters being written separately
  for clusid in list_of_clusters:
    fields=['substrate_ca_no_super_rmsd','distance_catalysis_HWUO1B-THR7N','interaction_energy','substrate_sequon_ca_no_super_rmsd','distance_catalysis']
    for field_name in fields:

      heatmap_type = heatmap.value
      dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],cluster_size_sorted_id=clusid,clus2d=clus2d,pattern=pattern)
      print(dictmat)
      heatmap_topN_sequon(dictmat , heatmap_type=heatmap_type, criteria=criteria , column_name=labels_df[field_name], filterTopN=True , show=show,roc=roc,cluster_size_sorted_id=clusid,suffix=suffix,clus2d=clus2d,annotate=annotate)

    heatmap_type = heatmap.size
    dictmat = get_topN_sequon_unglycosylated( criteria = criteria ,  heatmap_type =  heatmap_type,cluster_size_sorted_id=clusid,clus2d=clus2d,pattern=pattern)
    heatmap_topN_sequon(dictmat , heatmap_type=heatmap_type, criteria=criteria , filterTopN=True,show=show,roc=True,cluster_size_sorted_id=clusid,suffix=suffix,clus2d=clus2d,annotate=annotate)

   

def print_tprs_and_fprs(criteria,heatmap_type=heatmap.value,fields=[],pick_tpr=0.999,off=0.10,pick_fpr=None,pick_threshold=None,file_handle=open('accmetrics.out','w'),file_handle_auc=open('aucmetrics.out','w')):

  if heatmap_type==heatmap.value:

   for field_name in fields:

    suffix = '%s_' %field_name + '_off%f_pick_tpr%f' %(off,pick_tpr)
    if not pick_fpr is None:
      suffix = '%s_' %field_name + '_off%f_fpr%f' %(off,pick_fpr)
    if not pick_threshold is None:
      suffix = '%s_' %field_name + '_off%f_threshold%f' %(off,pick_threshold)

    prefix= 'results_aucscores/FullAccuracyMetrics_%s' %criteria
    outfile_fullmetrics = '%s_%s.out' %(prefix,suffix)
    dictmat_base,dictmat_compare,list_fps,list_fns,list_tps,list_tns = get_false_positives_and_false_negatives_for_criteria(criteria,heatmap_type,field_name,off=off,pick_tpr=pick_tpr,pick_threshold=pick_threshold,pick_fpr=pick_fpr,file_handle_auc=file_handle_auc,outfile=outfile_fullmetrics)
     
  else:
    heatmap_type = heatmap.size
    field_name = 'cluster_size'
    suffix = '%s_' %field_name + '_off%f_pick_tpr%f' %(off,pick_tpr)
    if not pick_fpr is None:
      suffix = '%s_' %field_name + '_off%f_fpr%f' %(off,pick_fpr)
    if not pick_threshold is None:
      suffix = '%s_' %field_name + '_off%f_threshold%f' %(off,pick_threshold)

    prefix= 'results_aucscores/FullAccuracyMetrics_%s' %criteria
    outfile_fullmetrics = '%s_%s.out' %(prefix,suffix)
    dictmat_base,dictmat_compare,list_fps,list_fns,list_tps,list_tns = get_false_positives_and_false_negatives_for_criteria(criteria,heatmap_type,field_name,off,pick_tpr=pick_tpr,pick_fpr=pick_fpr,pick_threshold=pick_threshold,file_handle_auc=file_handle_auc,outfile=outfile_fullmetrics)
    
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
  fpr=float(fp)/(tn+fp)
  bal_acc=0.5*(tpr+tnr)
  infd=tpr+tnr-1.0
  if not file_handle is None:
     file_handle.write('TP/TN/FP/FN %d/%d/%d/%d\n' %(tp,tn,fp,fn))
     file_handle.write('Acc/Prec/tpr/tnr/fpr %f/%f/%f/%f/%f\n'%(acc,prec,tpr,tnr,fpr))
     file_handle.write('bal_acc/infd %f/%f\n' %(bal_acc,infd))

def adjust_predictions(dictnew,key,list_tps,list_tns,list_fps,list_fns):
        #new prediction says non-glycosylatable
        if key in list_tns:
          #was correct still correct
          dictnew[key]=0.7 #correct #tns are 0.7
        if key in list_tps:
          #was correct now incorrect
          dictnew[key]=0.0 #incorrect #fns are 0.0
          list_tps.remove(key)
          list_fns.append(key)
        if key in list_fps:
          #was incorrect now correct
          dictnew[key]=0.7 #correct #tns are 0.7
          list_fps.remove(key)
          list_tns.append(key)
        if key in list_fns:
          #was incorrect still incorrect
          dictnew[key]=0.0 #incorrect #fns are 0.0
        return dictnew,list_tps,list_tns,list_fps,list_fns

def get_compare_heatmaps_GTX_TTY_states(pick_tpr=0.90,pick_fpr=None,pick_threshold=None):

  suffix = 'off%f_pick_tpr%f' %(off,pick_tpr)
  if not pick_fpr is None:
      suffix = 'off%f_fpr%f' %(off,pick_fpr)
  if not pick_threshold is None:
      suffix = 'off%f_threshold%f' %(off,pick_threshold)
  outfile_all = 'results/compare/compare_states_%s.txt' %suffix
  outf=open(outfile_all,'w')
  outf.write('#PTX_size\t')# tpr=%f\n' %pick_tpr)
  outfile = 'results/compare/compare_PTXsize_%s.txt' %(suffix)
  criteria=hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
  dictmat_PTX = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap.value , column_name = labels_df['interaction_energy'])
  dictmat_PTX_size = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap.size)
  _,dictmat_compare_PTX_size,list_fps,list_fns,list_tps,list_tns =  get_false_positives_and_false_negatives(dictmat_PTX_size,'cluster_size',off=0.10,pick_tpr=pick_tpr,pick_fpr=pick_fpr, pick_threshold = pick_threshold, outfile=outfile)
  heatmap_topN_sequon(dictmat_compare_PTX_size , heatmap_type=heatmap.compare, criteria=hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only , filterTopN=True ,suffix=suffix,annotate=False)

  write_metrics(list_tps,list_tns,list_fps,list_fns,file_handle=outf)
  assert((len(list_tps)+len(list_tns)+len(list_fps)+len(list_fns)) == 361)
  criteria = hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_rmsd
  dictmat_GTX = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap.value , column_name = labels_df['interaction_energy'])
  dictmat_GTX_size = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap.size)
  dictmat_GTX_delta = get_deltas(dictmat_PTX,dictmat_GTX)

  dictmat_compare_PTX_size_minusGTX = copy.deepcopy(dictmat_compare_PTX_size)
  for key in dictmat_PTX_size:
    if key in dictmat_GTX_delta:
      if (dictmat_GTX_size[key]>0.05) and (dictmat_GTX_delta[key]>1.5):
        dictmat_compare_PTX_size_minusGTX,list_tps,list_tns,list_fps,list_fns =adjust_predictions(dictmat_compare_PTX_size_minusGTX,key,list_tps,list_tns,list_fps,list_fns)

  assert((len(list_tps)+len(list_tns)+len(list_fps)+len(list_fns)) == 361)
  outf.write('GTX\t') #%d\t%d\t%d\t%d\n' %(len(list_tps),len(list_tns),len(list_fps),len(list_fns)))
  write_metrics(list_tps,list_tns,list_fps,list_fns,file_handle=outf)
  outf_1 = open('results/compare/Summary_classification_sequonlists_%s_%s_GTX' %(hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only,heatmap.size),'w')
  for lis,lis_name in zip([list_fps,list_fns],['FPs','FNs']):
    outf_1.write(lis_name+'\n')
    list_trim = [t[1:4] for t in lis]
    list_trim.sort()
    outf_1.write('\t'.join(list_trim)+'\n')
  outf_1.close()

  assert((len(list_tps)+len(list_tns)+len(list_fps)+len(list_fns)) == 361)
  print('minusGTX TPs,TNs,FPs,FNs:',len(list_tps),len(list_tns),len(list_fps),len(list_fns))
  suffix = 'PTXsize_minusGTX'
  heatmap_topN_sequon(dictmat_compare_PTX_size_minusGTX , heatmap_type=heatmap.compare, criteria=hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only , filterTopN=True ,suffix=suffix,annotate=False)
  
  criteria = hb_bond_criteria.significant_clusters_sinks_TTY_cutoff
  dictmat_TTY = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap.value , column_name = labels_df['interaction_energy'])
  dictmat_TTY_delta = get_deltas(dictmat_PTX,dictmat_TTY)
  dictmat_TTY_size = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap.size)
 
  dictmat_compare_PTX_size_minusTTY=copy.deepcopy(dictmat_compare_PTX_size)
  for key in dictmat_PTX_size:
    if key in dictmat_TTY_delta:
      if (dictmat_TTY_size[key]>0.05) and (dictmat_TTY_delta[key]>1.5):
        dictmat_compare_PTX_size_minusTTY,list_tps,list_tns,list_fps,list_fns =adjust_predictions(dictmat_compare_PTX_size_minusTTY,key,list_tps,list_tns,list_fps,list_fns)

  print('minusTTY TPs,TNs,FPs,FNs:',len(list_tps),len(list_tns),len(list_fps),len(list_fns))
  assert((len(list_tps)+len(list_tns)+len(list_fps)+len(list_fns)) == 361)
  outf.write('TTQ\t')#%d\t%d\t%d\t%d\n' %(len(list_tps),len(list_tns),len(list_fps),len(list_fns)))
  write_metrics(list_tps,list_tns,list_fps,list_fns,file_handle=outf)
  outf_1 = open('results/compare/Summary_classification_sequonlists_%s_%s_TTY' %(hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only,heatmap.size),'w')
  for lis,lis_name in zip([list_fps,list_fns],['FPs','FNs']):
    outf_1.write(lis_name+'\n')
    list_trim = [t[1:4] for t in lis]
    list_trim.sort()
    outf_1.write('\t'.join(list_trim)+'\n')
  outf_1.close()

  suffix = 'PTXsize_minusTTY' 
  heatmap_topN_sequon(dictmat_compare_PTX_size_minusTTY , heatmap_type=heatmap.compare, criteria=hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only , filterTopN=True ,suffix=suffix,annotate=False)

  dictmat_compare_PTX_size_minusGTX_minusTTY = copy.deepcopy(dictmat_compare_PTX_size_minusGTX)

  for key in dictmat_PTX_size:
    if key in dictmat_TTY_delta:
      if (dictmat_TTY_size[key]>0.05) and (dictmat_TTY_delta[key]>1.5):
        dictmat_compare_PTX_size_minusGTX_minusTTY,list_tps,list_tns,list_fps,list_fns =adjust_predictions(dictmat_compare_PTX_size_minusGTX_minusTTY,key,list_tps,list_tns,list_fps,list_fns)
  assert((len(list_tps)+len(list_tns)+len(list_fps)+len(list_fns)) == 361)
  outf.write('GTX+TTQ\t')#%d\t%d\t%d\t%d\n' %(len(list_tps),len(list_tns),len(list_fps),len(list_fns)))
  write_metrics(list_tps,list_tns,list_fps,list_fns,file_handle=outf)
  outf.close()

  outf_1 = open('results/compare/Summary_classification_sequonlists_%s_%s_GTX_TTY' %(hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only,heatmap.size),'w')
  for lis,lis_name in zip([list_fps,list_fns],['FPs','FNs']):
    outf_1.write(lis_name+'\n')
    list_trim = [t[1:4] for t in lis]
    list_trim.sort()
    outf_1.write('\t'.join(list_trim)+'\n')
  outf_1.close()
  assert((len(list_tps)+len(list_tns)+len(list_fps)+len(list_fns)) == 361)
  print('minusGTX_minusTTY TPs,TNs,FPs,FNs:',len(list_tps),len(list_tns),len(list_fps),len(list_fns))
  suffix = 'PTXsize_minusGTX_minusTTY'
  heatmap_topN_sequon(dictmat_compare_PTX_size_minusGTX_minusTTY, heatmap_type=heatmap.compare, criteria=hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only , filterTopN=True ,suffix=suffix,annotate=False)

def print_rocs(clus2d=False,off=0.1,pick_tpr=0.999,threshold_GTX=0.45,threshold_TTX=0.45):
  criteria = hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
  heatmap_type = heatmap.size
  field_name = 'cluster_size'
  dictmat_base,list_fps_1,list_fns_1,list_tps_1 = get_false_positives_and_false_negatives(criteria,heatmap_type,field_name,off,pick_tpr)
 
  criteria = hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_rmsd
  heatmap_type = heatmap.size
  field_name = 'cluster_size'
  dictmat = get_topN_sequon_unglycosylated( criteria = criteria ,  heatmap_type =  heatmap_type)
  list_keep = []
  list_remove = []
  for key in list_fps_1:
    if dictmat[key] < threshold_GTX:
      list_keep.append(key)
    else:
      list_remove.append(key)
      del dictmat_base[key]
  print('GTX remove: ',len(list_remove),'keep: ',len(list_keep))
  list_keep_ = []
  list_remove_ = []
  for key in list_tps_1:
    if dictmat[key] < threshold_GTX:
      list_keep_.append(key)
    else:
      list_remove_.append(key)
      print(key,dictmat[key])
      del dictmat_base[key]
  print('GTX remove true positives: ',len(list_remove_))

  criteria = hb_bond_criteria.significant_clusters_sinks_TTY_cutoff
  heatmap_type = heatmap.size
  field_name = 'cluster_size'
  dictmat = get_topN_sequon_unglycosylated( criteria = criteria ,  heatmap_type =  heatmap_type)
  list_keep_1 = []
  list_remove_1 = []
  for key in list_keep:
    if dictmat[key] < threshold_TTX:
      list_keep_1.append(key)
      print(key,dictmat[key])
    else:
      list_remove_1.append(key)
      del dictmat_base[key]
  print('TTX remove: ',len(list_remove_1),'keep: ',len(list_keep_1))
  list_keep_ = []
  list_remove_ = []
  for key in list_tps_1:
    if dictmat[key] < threshold_TTX:
      list_keep_.append(key)
    else:
      list_remove_.append(key)
      del dictmat_base[key]
      print(key,dictmat[key])
    
  print('TTX remove true positives: ',len(list_remove_))
  criteria = hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
  heatmap_topN_sequon(dictmat_base , heatmap_type=heatmap_type, criteria= criteria, filterTopN=True,show=True,roc=True,suffix='rmsd_GTX_TTX',annotate=False,off=off)


def get_heatmaps_for_criteria(criteria,show=True,roc=False,suffix=None,clus2d=False,annotate=False,topN=200,pattern=None,listres1=None):

    for field_name in ['interaction_energy','substrate_ca_no_super_rmsd','distance_catalysis_HWUO1B-THR7N']:#,'distance_catalysis','interaction_energy','substrate_sequon_ca_no_super_rmsd']:#labels_df:

      heatmap_type = heatmap.value
      if pattern is None:
       dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],clus2d=clus2d,topN=topN)
      else:
       dictmat = get_topN_sequon_unglycosylated( pattern =pattern , criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],clus2d=clus2d,topN=topN)
      print(dictmat)
      if listres1 is None:
        heatmap_topN_sequon(dictmat , heatmap_type=heatmap_type, criteria=criteria , column_name=labels_df[field_name], filterTopN=True , show=show,roc=roc,suffix=suffix,clus2d=clus2d,annotate=annotate)
      else:
        heatmap_topN_sequon(dictmat , heatmap_type=heatmap_type, criteria=criteria , column_name=labels_df[field_name], filterTopN=True , show=show,roc=roc,suffix=suffix,clus2d=clus2d,annotate=annotate,listres1=listres1)

def write_arrays_for_criteria_positionwise(criteria,show=True,roc=False,suffix=None,clus2d=False,annotate=False,topN=200):
    os.system('mkdir -p results_positionwise')
    filename = 'results_positionwise/aucscores_%s_top%d.txt' %(criteria,topN)
    outf=open(filename,'w')
    for field_name in ['interaction_energy','substrate_ca_no_super_rmsd','distance_catalysis_HWUO1B-THR7N','distance_catalysis','substrate_sequon_ca_no_super_rmsd']:#labels_df:
      heatmap_type = heatmap.value
      dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],clus2d=clus2d,topN=topN)
      get_positionwise_roc_auc_for_dict(dictmat,value_name=field_name,file_handle=outf)
    field_name='cluster_size'
    dictmat = get_topN_sequon_unglycosylated( criteria = criteria ,  heatmap_type =  heatmap.size,clus2d=clus2d,topN=topN)
    get_positionwise_roc_auc_for_dict(dictmat,value_name=field_name,file_handle=outf)

    outf.close()


def write_arrays_for_criteria(criteria,show=True,roc=False,suffix=None,clus2d=False,annotate=False,topN=200):
    os.system('mkdir -p results_aucscores')
    filename = 'results_aucscores/aucscores_%s_top%d.txt' %(criteria,topN)
    outf=open(filename,'w')
    outf.write('#off,field_name,auc\n')
    for field_name in ['interaction_energy','substrate_ca_no_super_rmsd','distance_catalysis_HWUO1B-THR7N','distance_catalysis','substrate_sequon_ca_no_super_rmsd']:#labels_df:
      heatmap_type = heatmap.value
      dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],clus2d=clus2d,topN=topN)

      write_roc_auc_for_dict(dictmat,value_name=field_name,file_handle=outf)
    if criteria== hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_rmsd:
      field_name='cluster_size_GTX'
    else:
      field_name='cluster_size'
    dictmat = get_topN_sequon_unglycosylated( criteria = criteria ,  heatmap_type =  heatmap.size,clus2d=clus2d,topN=topN)
    write_roc_auc_for_dict(dictmat,value_name=field_name,file_handle=outf)

    outf.close()

def write_scores_for_criteria():
  criteria = hb_bond_criteria.significant_clusters_sinks_cutoff
  write_arrays_for_criteria(criteria)
  criteria = hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
  write_arrays_for_criteria(criteria)
  get_heatmaps_for_criteria_cutoff(criteria=criteria,heatmap_type=heatmap.size)
  criteria = hb_bond_criteria.significant_clusters_sinks
  print_tprs_and_fprs(criteria)
  write_arrays_for_criteria(criteria)

def get_heatmaps_for_criteria_cutoff(show=True,criteria=hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only,heatmap_type = heatmap.size,suffix=None,clus2d=False,topN=200,pattern=None,listres1=None,roc=False,annotate=False):

  if pattern is None:
   dictmat = get_topN_sequon_unglycosylated( criteria = criteria ,  heatmap_type =  heatmap_type,clus2d=clus2d,topN=topN)
  else:
   dictmat = get_topN_sequon_unglycosylated(pattern=pattern, criteria = criteria ,  heatmap_type =  heatmap_type,clus2d=clus2d,topN=topN)
  if listres1 is None:
   heatmap_topN_sequon(dictmat , heatmap_type=heatmap_type, criteria=criteria , filterTopN=True , show=show,suffix=suffix,clus2d=clus2d,annotate=annotate,roc=roc)
  else:
   heatmap_topN_sequon(dictmat , heatmap_type=heatmap_type, criteria=criteria , filterTopN=True , show=show,suffix=suffix,clus2d=clus2d,annotate=annotate,listres1=listres1,roc=roc)

def get_heatmaps_for_PTX_GTX_criteria(show=False,roc=False,annotate=False,pattern=None,listres1=None):

    for field_name in ['interaction_energy']:#,'distance_catalysis_HWUO1B-THR7N','substrate_ca_no_super_rmsd','substrate_sequon_ca_no_super_rmsd']:#labels_df:

      heatmap_type = heatmap.value
      criteria = hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_rmsd
      if pattern is None:
         dictmat_GTX = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name])
      else:
         dictmat_GTX = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],pattern=pattern)
      #print('dict GTX',dictmat_GTX)
      #heatmap_topN_sequon(dictmat_GTX , heatmap_type=heatmap_type, criteria=criteria, column_name=labels_df[field_name],filterTopN=True , show=show,roc=roc,suffix='GTXonly')

      criteria = hb_bond_criteria.significant_clusters_sinks_TTY_cutoff
      dictmat_TTY = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name])
      #heatmap_topN_sequon(dictmat_TTY , heatmap_type=heatmap_type, criteria=criteria, column_name=labels_df[field_name],filterTopN=True , show=show,roc=roc,suffix='GTXonly')

      if field_name == 'interaction_energy':
        criteria = hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
        dictmat_PTX = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name])
        heatmap_type = heatmap.delta
        dictmat_GTX = get_deltas(dictmat_PTX,dictmat_GTX)
        dictmat_TTY = get_deltas(dictmat_PTX,dictmat_TTY)
        if listres1 is None:
          heatmap_topN_sequon(dictmat_GTX , heatmap_type=heatmap_type, criteria=None, filterTopN=True , show=show,roc=False,suffix='delta_PTX_GTX',annotate=annotate)
          heatmap_topN_sequon(dictmat_TTY , heatmap_type=heatmap_type, criteria=None, filterTopN=True , show=show,roc=False,suffix='delta_PTX_TTY',annotate=annotate)
        else:
          heatmap_topN_sequon(dictmat_GTX , heatmap_type=heatmap_type, criteria=None, filterTopN=True , show=show,roc=False,suffix='delta_PTX_GTX',annotate=annotate,listres1=listres1)
          heatmap_topN_sequon(dictmat_TTY , heatmap_type=heatmap_type, criteria=None, filterTopN=True , show=show,roc=False,suffix='delta_PTX_TTY',annotate=annotate,listres1=listres1)


def get_experimental_heatmap(show=False,annotate=False):
  import compare_simulations_to_experiments as CompSE
  exp_dictmat = CompSE.getexperimentaldatafrompickle_dict(percent=True)
  heatmap_type=heatmap.exp
  heatmap_topN_sequon(exp_dictmat , heatmap_type=heatmap_type, criteria=None, filterTopN=False , show=show,roc=False,suffix='exp2018',annotate=annotate)


def generate_heatmaps_for_all_criteria(annotate=False,clus2d=False,roc=False):
  all_criteria = [hb_bond_criteria.significant_clusters, #centroid of cluster
              hb_bond_criteria.significant_clusters_sinks, #lowest scoring decoy of cluster
              hb_bond_criteria.significant_clusters_sinks_TopKCal, #average over lowest scoring decoys within 1REU of lowest score
              hb_bond_criteria.significant_clusters_sinks_cutoff, #dHB criterion
              hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only #rmsd criterion
              ]
  for criteria in all_criteria:
    get_heatmaps_for_criteria_for_all_clusters(criteria,show=False,annotate=annotate)

  pattern = []
  listres1 = ['A','G','S','T']
  for aa in ['ALA','GLY','SER','THR']:
    pattern += ['PackandRelax_498%s*/unglycosylated/score.sc' %aa]
  get_heatmaps_for_PTX_GTX_criteria(show=False,roc=roc,annotate=annotate,pattern=pattern,listres1=listres1)
  for criteria in all_criteria:
    get_heatmaps_for_criteria(criteria, show = False , roc=roc,clus2d=clus2d,annotate=annotate,pattern=pattern,listres1=listres1)
    get_heatmaps_for_criteria_cutoff(criteria=criteria,heatmap_type=heatmap.size,pattern=pattern,listres1=listres1,roc=roc)

  get_experimental_heatmap(annotate=annotate)

def accuracy_metrics(off=0.10):
  criteria=hb_bond_criteria.significant_clusters_sinks
  file_handle_auc=open('results_aucscores/AUCScores_%s.txt' %(criteria),'w')
  fields=['distance_catalysis_HWUO1B-THR7N','substrate_ca_no_super_rmsd','substrate_sequon_ca_no_super_rmsd','interaction_energy']
  thresholds=[-4.0,-1.0,-0.9,34.0]
  for field,pick_threshold in zip(fields,thresholds):
    print_tprs_and_fprs(criteria,fields=[field],off=off,pick_threshold=pick_threshold,file_handle_auc=file_handle_auc)
  file_handle_auc.close()


  criterias=[hb_bond_criteria.significant_clusters_sinks,
            hb_bond_criteria.significant_clusters_sinks_cutoff,
            hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only]
  for criteria in criterias:
    file_handle_auc=open('results_aucscores/AUCScores_%s_IEonly.txt' %(criteria),'w')
    field='interaction_energy'
    for pick_threshold in thresholds:
      print_tprs_and_fprs(criteria,fields=[field],off=off,pick_threshold=pick_threshold,file_handle_auc=file_handle_auc)
  file_handle_auc.close()

  criterias=[hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only]
  fprs=[0.075,0.08,0.09,0.10,0.11,0.12]
  for criteria in criterias:
    field='cluster_size'
    file_handle_auc=open('results_aucscores/AUCScores_%s_fprs_%s.txt' %(criteria,field),'w')
    off=0.10
    for pick_fpr in fprs:
      print_tprs_and_fprs(criteria,heatmap_type=heatmap.size,off=off,pick_fpr=pick_fpr,file_handle_auc=file_handle_auc)
    file_handle_auc.close()
  criterias=[hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only]
  tprs=[0.89,0.90,0.91,0.92,0.93]
  for criteria in criterias:
    field='cluster_size'
    file_handle_auc=open('results_aucscores/AUCScores_%s_tprs_%s.txt' %(criteria,field),'w')
    off=0.10
    for pick_tpr in tprs:
      print_tprs_and_fprs(criteria,heatmap_type=heatmap.size,off=off,pick_tpr=pick_tpr,file_handle_auc=file_handle_auc)
    file_handle_auc.close()

  criterias=[hb_bond_criteria.significant_clusters_sinks_cutoff]
  for criteria in criterias:
    field='cluster_size'
    file_handle_auc=open('results_aucscores/AUCScores_%s_fprs_%s.txt' %(criteria,field),'w')
    fprs=[0.143,0.14,0.13,0.16,0.17,0.18]
    off=0.10
    for pick_fpr in fprs:
      print_tprs_and_fprs(criteria,heatmap_type=heatmap.size,off=off,pick_fpr=pick_fpr,file_handle_auc=file_handle_auc)
    file_handle_auc.close()

def cli():
  annotate=False
  clus2d=False
  roc=False
  generate_heatmaps_for_all_criteria(annotate=annotate,clus2d=clus2d,roc=roc)
  get_compare_heatmaps_GTX_TTY_states()
  write_scores_for_criteria() #write auc scores
  accuracy_metrics()

if __name__ == '__main__':
  #cli()
  #annotate=False
  #clus2d=False
  #roc=False
  #generate_heatmaps_for_all_criteria(annotate=annotate,clus2d=clus2d,roc=roc)
  #get_compare_heatmaps_GTX_TTY_states()
  write_scores_for_criteria() #write auc scores
  #accuracy_metrics()
