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

  list_of_clusters = [0]#,1] #only matters for signicant clusters being written separately
  for clusid in list_of_clusters:
    fields=['substrate_ca_no_super_rmsd']#,'distance_catalysis_HWUO1B-THR7N','interaction_energy','substrate_sequon_ca_no_super_rmsd','distance_catalysis']
    for field_name in fields:

      heatmap_type = heatmap.value
      dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],cluster_size_sorted_id=clusid,clus2d=clus2d,pattern=pattern)
      print(dictmat)
      heatmap_topN_sequon(dictmat , heatmap_type=heatmap_type, criteria=criteria , column_name=labels_df[field_name], filterTopN=True , show=show,roc=roc,cluster_size_sorted_id=clusid,suffix=suffix,clus2d=clus2d,annotate=annotate)

    heatmap_type = heatmap.size
    dictmat = get_topN_sequon_unglycosylated( criteria = criteria ,  heatmap_type =  heatmap_type,cluster_size_sorted_id=clusid,clus2d=clus2d,pattern=pattern)
    heatmap_topN_sequon(dictmat , heatmap_type=heatmap_type, criteria=criteria , filterTopN=True,show=show,roc=True,cluster_size_sorted_id=clusid,suffix=suffix,clus2d=clus2d,annotate=annotate)

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
    #write_metrics(list_tps,list_tns,list_fps,list_fns,file_handle=file_handle) 
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
    #write_metrics(list_tps,list_tns,list_fps,list_fns,file_handle=file_handle)


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


def get_compare_heatmaps_dhb():

  #criteria = hb_bond_criteria.significant_clusters_sinks_cutoff
  #dictmat_distance = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap.value , column_name = labels_df['distance_catalysis_HWUO1B-THR7N'])
  criteria = hb_bond_criteria.significant_clusters_sinks
  dictmat_distance = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap.value , column_name = labels_df['distance_catalysis_HWUO1B-THR7N'],cluster_size_sorted_id=0)
  pick_threshold=-4.0
  suffix = '_%s_threshold%f' %(criteria,pick_threshold)
  outfile = 'results/compare/compare%s.txt' %(suffix)
  _,dictmat_compare_distance,list_fps,list_fns,list_tps,list_tns =  get_false_positives_and_false_negatives(dictmat_distance,'distance_catalysis_HWUO1B-THR7N',off=0.10,pick_threshold=pick_threshold,outfile=outfile)

  outf_1 = open('results/compare/Summary_classification_sequonlists_%s_%s_threshold%f' %(criteria,heatmap.size,pick_threshold),'w')
  for lis,lis_name in zip([list_fps,list_fns],['FPs','FNs']):
       outf_1.write(lis_name+'\n')
       list_trim = [t[1:4] for t in lis]
       list_trim.sort()
       outf_1.write('\t'.join(list_trim)+'\n')
  outf_1.close()

  suffix = 'distance_catalysis_HWUO1B-THR7N_threshold%f' %pick_threshold
  heatmap_topN_sequon(dictmat_compare_distance , heatmap_type=heatmap.compare, criteria=criteria, filterTopN=True ,suffix=suffix,annotate=False)

  for pick_tpr in [0.91,0.93,0.95,0.97,0.99]:
     outfile = 'results/compare/compare%s.txt' %(suffix)
     _,dictmat_compare_distance,list_fps,list_fns,list_tps,list_tns =  get_false_positives_and_false_negatives(dictmat_distance,'distance_catalysis_HWUO1B-THR7N',off=0.10,pick_tpr=pick_tpr,outfile=outfile)

     outf_1 = open('results/compare/Summary_classification_sequonlists_%s_%s' %(criteria,heatmap.size),'w')
     for lis,lis_name in zip([list_fps,list_fns],['FPs','FNs']):
       outf_1.write(lis_name+'\n')
       list_trim = [t[1:4] for t in lis]
       list_trim.sort()
       outf_1.write('\t'.join(list_trim)+'\n')
     outf_1.close()

     suffix = '_distance_catalysis_HWUO1B-THR7N_tpr%f' %pick_tpr
     heatmap_topN_sequon(dictmat_compare_distance , heatmap_type=heatmap.compare, criteria=criteria, filterTopN=True ,suffix=suffix,annotate=False)

def get_compare_heatmaps_rmsd():
  criteria=hb_bond_criteria.significant_clusters_sinks
  dictmat_rmsd = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap.value , column_name = labels_df['substrate_ca_no_super_rmsd'],cluster_size_sorted_id=0)
  pick_threshold=-1.0
  suffix = '_%s_threshold%f' %(criteria,pick_threshold)
  outfile = 'results/compare/compare%s.txt' %(suffix)
  _,dictmat_compare_rmsd,list_fps,list_fns,list_tps,list_tns =  get_false_positives_and_false_negatives(dictmat_rmsd,'substrate_ca_no_super_rmsd',off=0.10,pick_threshold=pick_threshold,outfile=outfile)

  outf_1 = open('results/compare/Summary_classification_sequonlists_%s_%s_threshold%f' %(criteria,heatmap.size,pick_threshold),'w')
  for lis,lis_name in zip([list_fps,list_fns],['FPs','FNs']):
    outf_1.write(lis_name+'\n')
    list_trim = [t[1:4] for t in lis]
    list_trim.sort()
    outf_1.write('\t'.join(list_trim)+'\n')
  outf_1.close()

  suffix = '_substrate_ca_no_super_rmsd_threshold%f' %pick_threshold
  heatmap_topN_sequon(dictmat_compare_rmsd , heatmap_type=heatmap.compare, criteria=criteria, filterTopN=True ,suffix=suffix,annotate=False)

  for pick_tpr in [0.91,0.93,0.95,0.97,0.99]:
     suffix='_substrate_ca_no_super_rmsd_tpr%f' %pick_tpr
     outfile = 'results/compare/compare%s.txt' %(suffix)
     _,dictmat_compare_rmsd,list_fps,list_fns,list_tps,list_tns =  get_false_positives_and_false_negatives(dictmat_rmsd,'substrate_ca_no_super_rmsd_',off=0.10,pick_tpr=pick_tpr,outfile=outfile)

     outf_1 = open('results/compare/Summary_classification_sequonlists_%s_%s%s' %(criteria,heatmap.size,suffix),'w')
     for lis,lis_name in zip([list_fps,list_fns],['FPs','FNs']):
       outf_1.write(lis_name+'\n')
       list_trim = [t[1:4] for t in lis]
       list_trim.sort()
       outf_1.write('\t'.join(list_trim)+'\n')
     outf_1.close()
     heatmap_topN_sequon(dictmat_compare_rmsd , heatmap_type=heatmap.compare, criteria=criteria, filterTopN=True ,suffix=suffix,annotate=False)

def get_compare_interaction_energy():

  criteria=hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
  dictmat_PTX_size = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap.size)
  heatmap_type = heatmap.value
  field_name='interaction_energy'
  clusid=0
  dictmat_0 = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],cluster_size_sorted_id=clusid,clus2d=clus2d) 
  clusid=1
  dictmat_1 = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],cluster_size_sorted_id=clusid,clus2d=clus2d)
  keys = list(set(list(dictmat_0.keys())+list(dictmat_1.keys())))
  dictmat_1={}
  for key in keys:
     dictmat_lowest[key]=99.0
     if key in dictmat_0:
        dictmat_lowest[key]=dictmat_0[key]
     if key in dictmat_1:
        if dictmat_1[key]<dictmat_lowest[key]:
           dictmat_lowest[key]=dictmat_1[key]

  

def get_compare_heatmaps_other_states(pick_tpr=0.90,pick_fpr=None,pick_threshold=None):

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


def get_rocs_for_values(clus2d=False,combined_id=1):
  dictmats=[]
  value_names = []
  if combined_id==1:
    criteria=hb_bond_criteria.significant_clusters_sinks_cutoff
    field_name = 'distance_catalysis_HWUO1B-THR7N'
    heatmap_type = heatmap.value
    dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],clus2d=clus2d)
    dictmats.append(dictmat)
    value_names.append(field_name)
    heatmap_type = heatmap.size
    dictmat = get_topN_sequon_unglycosylated( criteria = criteria ,  heatmap_type =  heatmap_type)
    dictmats.append(dictmat)
    value_names.append('cluster_size')
    name_suffix = 'dhb_and_clustersize'

  if combined_id==2:
    criteria=hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd
    field_name = 'substrate_ca_no_super_rmsd'
    heatmap_type = heatmap.value
    dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],clus2d=clus2d)
    dictmats.append(dictmat)
    value_names.append(field_name)
    heatmap_type = heatmap.size
    dictmat = get_topN_sequon_unglycosylated( criteria = criteria ,  heatmap_type =  heatmap_type)
    dictmats.append(dictmat)
    value_names.append('cluster_size')
    name_suffix = 'rmsd_and_clustersize'

  if combined_id==3:
    criteria=hb_bond_criteria.significant_clusters_cutoff_rmsd
    field_name = 'substrate_ca_no_super_rmsd'
    heatmap_type = heatmap.value
    dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],clus2d=clus2d)
    dictmats.append(dictmat)
    value_names.append(field_name)
    heatmap_type = heatmap.size
    dictmat = get_topN_sequon_unglycosylated( criteria = criteria ,  heatmap_type =  heatmap_type)
    dictmats.append(dictmat)
    value_names.append('cluster_size')
    name_suffix = 'rmsd_and_clustersize'

  if combined_id==4:

    heatmap_type = heatmap.value
    criteria=hb_bond_criteria.significant_clusters_sinks
    field_name = 'interaction_energy'
    dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],cluster_size_sorted_id=0)

    dictmats.append(dictmat)
    value_names.append(field_name)

    criteria=hb_bond_criteria.significant_clusters_sinks_cutoff
    dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],clus2d=clus2d)

    dictmats.append(dictmat)
    value_names.append(field_name)

    criteria=hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
    dictmat = get_topN_sequon_unglycosylated( criteria = criteria ,  heatmap_type =  heatmap_type)

    dictmats.append(dictmat)
    value_names.append(field_name)
    name_suffix = 'interaction_energy_3criteria'

  if combined_id==5:

    heatmap_type = heatmap.value
    criteria=hb_bond_criteria.significant_clusters_sinks
    field_name = 'interaction_energy'
    dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],cluster_size_sorted_id=0)

    dictmats.append(dictmat)
    value_names.append(field_name)

    criteria=hb_bond_criteria.significant_clusters_sinks_cutoff
    dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name],clus2d=clus2d)

    dictmats.append(dictmat)
    value_names.append(field_name)

    criteria=hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
    dictmat = get_topN_sequon_unglycosylated( criteria = criteria ,  heatmap_type =  heatmap_type)

    dictmats.append(dictmat)
    value_names.append(field_name)

    #field_name = 'fa_atr'
    #from functions_additional_heatmaps_paper import get_dict_stat_for_res_pair
    #dictmat = get_dict_stat_for_res_pair(res_i = 498,res_j=291,st=field_name,skey='median',rmsd_cutoff=100.0)
    #dictmats.append(dictmat)
    #value_names.append(field_name)

    field_name = 'fa_atr'
    from functions_additional_heatmaps_paper import get_dict_stat_for_res_pair 
    dictmat = get_dict_stat_for_res_pair(res_i = 498,res_j=291,st=field_name,skey='median',N=10)
    dictmats.append(dictmat)
    value_names.append(field_name)
 
 
    name_suffix = 'interaction_energy_4criteria'

  plot_roc_auc_for_dicts(dictmats,name_suffix=name_suffix,value_names=value_names,clus2d=clus2d,show=False)

  


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

def get_features_for_criteria(criteria,show=True,roc=False,suffix=None):

    listofdicts = []
    list_names = []
    for field_name in ['interaction_energy']:#,'substrate_sequon_ca_no_super_rmsd']:#,'interaction_energy','substrate_ca_no_super_rmsd']:#labels_df:

      
      heatmap_type = heatmap.value
      dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name])
      print(dictmat)
      list_names.append(field_name)
      listofdicts.append(dictmat)

      '''
      criteria_GTX = hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_rmsd
      dictmat_GTX = get_topN_sequon_unglycosylated( criteria=criteria_GTX , heatmap_type=heatmap_type , column_name = labels_df[field_name])
      listofdicts.append(dictmat_GTX)
      list_names.append(field_name+'_GTX')

      criteria_TTY = hb_bond_criteria.significant_clusters_sinks_TTY_cutoff
      dictmat_TTY = get_topN_sequon_unglycosylated( criteria=criteria_TTY , heatmap_type=heatmap_type , column_name = labels_df[field_name])
      listofdicts.append(dictmat_TTY)
      list_names.append(field_name+'_TTY')
    
      '''

    field_name = 'interaction_energy'
    criteria_GTX = hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_rmsd
    dictmat_GTX = get_topN_sequon_unglycosylated( criteria=criteria_GTX , heatmap_type=heatmap_type , column_name = labels_df[field_name])
    listofdicts.append(dictmat_GTX)
    list_names.append(field_name+'_GTX')

    criteria_TTY = hb_bond_criteria.significant_clusters_sinks_TTY_cutoff
    dictmat_TTY = get_topN_sequon_unglycosylated( criteria=criteria_TTY , heatmap_type=heatmap_type , column_name = labels_df[field_name])
    listofdicts.append(dictmat_TTY)
    list_names.append(field_name+'_TTY')

    print(dictmat_TTY)
    #list_names.append('cluster_size')
    #heatmap_type = heatmap.size
    #dictmat = get_topN_sequon_unglycosylated( criteria = criteria ,  heatmap_type =  heatmap_type)
    #listofdicts.append(dictmat)
    feature_importance(listofdicts,list_value_names=list_names)

def get_heatmaps_for_criteria_cutoff(show=True,criteria=hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_rmsd,heatmap_type = heatmap.size,suffix=None,clus2d=False,topN=200,pattern=None,listres1=None,roc=False):

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


def get_heatmaps_glycosylated(show=True,roc=True):
  heatmap_type = heatmap.value
  for field_name in ['substrate_ca_no_super_rmsd','distance_catalysis_UDP2OPB-THR7N','interaction_energy','substrate_sequon_ca_no_super_rmsd']:
    dictmat = get_topN_sequon_glycosylated(column_name=labels_df[field_name])
    heatmap_topN_sequon(dictmat , heatmap_type=heatmap_type, column_name=labels_df[field_name], filterTopN=True , show=show, roc=roc)

def get_experimental_heatmap(show=False,annotate=False):
  import compare_simulations_to_experiments as CompSE
  exp_dictmat = CompSE.getexperimentaldatafrompickle_dict(percent=True)
  heatmap_type=heatmap.exp
  heatmap_topN_sequon(exp_dictmat , heatmap_type=heatmap_type, criteria=None, filterTopN=False , show=show,roc=False,suffix='exp2018',annotate=annotate)

def experimental_vs_simulations(criteria,heatmap_type,field_name=''):
  import compare_simulations_to_experiments as CompSE
  _,exp_mat, exp_mat_onoff = CompSE.getexperimentaldatafrompickle()
  dictmat = {}
  if heatmap_type==heatmap.value:
    dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name])  
  if heatmap_type==heatmap.size:
    dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type)
  
  _,mat = mat_for_dict(dictmat,field_name)
  suffix = '%s_%s' %(criteria,field_name)
  #CompSE.scatterplot_A_vs_B(mat,exp_mat,suffix=suffix) 
  X_ravel = mat.ravel()
  Y_ravel = exp_mat.ravel()
  fig,ax = plt.subplots(figsize=(16,16))
  ax.spines['bottom'].set_color('black')
  ax.spines['top'].set_color('black')
  ax.spines['right'].set_color('black')
  ax.spines['left'].set_color('black')
  fig.suptitle(datestring + '\n%s' %suffix,fontsize=10)
  plt.plot([0,1],[0,1],color='dimgrey',linestyle='dashed')
  plt.scatter(X_ravel, Y_ravel,color='xkcd:aquamarine',s=100,edgecolor='black',alpha=1)
  plt.xticks(fontsize=60)
  plt.yticks(fontsize=60)
  plt.ylim(-0.02,1.02)
  plt.xlim(-0.02,1.02)
  axis_font = { 'size':70}
  plt.xlabel('predictions',**axis_font)
  plt.ylabel('experiments',**axis_font)
  plt.tight_layout()
  plt.subplots_adjust(top=0.9)
  plt.savefig("results/compare/XvsY_"+ suffix + ".png",transparent=True,dpi=fig.dpi)
  plt.show()


annotate=False
clus2d=False
roc=False
criteria=hb_bond_criteria.significant_clusters_sinks_TopKCal
#get_heatmaps_for_criteria_for_all_clusters(criteria,show=False,annotate=False)
#get_heatmaps_for_criteria_for_all_clusters(criteria,show=False,annotate=True)
criteria=hb_bond_criteria.significant_clusters
#get_heatmaps_for_criteria_for_all_clusters(criteria,show=False,annotate=True)
#get_heatmaps_for_criteria(criteria, show = True , roc=True,clus2d=clus2d,annotate=annotate,topN=50)
#criteria = hb_bond_criteria.significant_clusters_sinks_TTY_cutoff
#criteria = hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_rmsd
#criteria = hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
#get_heatmaps_for_criteria(criteria, show = True , roc=True,clus2d=clus2d,annotate=annotate)
pattern = []
listres1 = ['A','G','S','T']
for aa in ['ALA','GLY','SER','THR']:
  pattern += ['PackandRelax_498%s*/unglycosylated/score.sc' %aa]
#get_heatmaps_for_PTX_GTX_criteria(show=False,roc=False,annotate=False,pattern=pattern,listres1=listres1)
#get_heatmaps_for_criteria(criteria, show = False , roc=roc,clus2d=clus2d,annotate=annotate,pattern=pattern,listres1=listres1)
#get_heatmaps_for_criteria_cutoff(criteria=criteria,heatmap_type=heatmap.size,pattern=pattern,listres1=listres1,roc=False)
#get_experimental_heatmap(annotate=True)
get_compare_heatmaps_other_states()
exit()
#get_heatmaps_for_criteria(criteria, show = False , roc=roc,clus2d=clus2d,annotate=annotate,pattern=pattern,listres1=listres1)
#get_heatmaps_for_criteria_cutoff(criteria=criteria,heatmap_type=heatmap.size,pattern=pattern,listres1=listres1,roc=False)
#criteria = hb_bond_criteria.significant_clusters_sinks_cutoff
#get_heatmaps_for_criteria(criteria, show = False , roc=roc,clus2d=clus2d,annotate=annotate)
#get_heatmaps_for_criteria_cutoff(criteria=criteria,heatmap_type=heatmap.size,roc=False)
#criteria = hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
#get_heatmaps_for_criteria(criteria, show = False , roc=roc,clus2d=clus2d,annotate=annotate)
#get_heatmaps_for_criteria_cutoff(criteria=criteria,heatmap_type=heatmap.size,roc=False)

criteria=hb_bond_criteria.significant_clusters_sinks
file_handle_auc=open('results_aucscores/AUCScores_%s.txt' %(criteria),'w')
fields=['distance_catalysis_HWUO1B-THR7N','distance_catalysis_HWUO1B-THR7N','substrate_ca_no_super_rmsd','substrate_sequon_ca_no_super_rmsd','interaction_energy']
thresholds=[-4.0,-3.7,-1.0,-0.9,34.0]
off=0.10
for field,pick_threshold in zip(fields,thresholds):
 print_tprs_and_fprs(criteria,fields=[field],off=off,pick_threshold=pick_threshold,file_handle_auc=file_handle_auc)
file_handle_auc.close()
'''
fields=['distance_catalysis_HWUO1B-THR7N','distance_catalysis_HWUO1B-THR7N','substrate_ca_no_super_rmsd','substrate_ca_no_super_rmsd','substrate_sequon_ca_no_super_rmsd','substrate_sequon_ca_no_super_rmsd']
thresholds=[-3.3,-3.2,-0.7,-0.68,-0.65,-0.6]
off=0.55
file_handle.write('\n\n')
for field,pick_threshold in zip(fields,thresholds):
 file_handle.write('%s off threshold %f %f\n' %(field,off,pick_threshold))
 print_tprs_and_fprs(criteria,fields=[field],off=off,pick_threshold=pick_threshold,file_handle=file_handle,file_handle_auc=file_handle_auc)
file_handle.close()
file_handle_auc.close()
'''
criterias=[hb_bond_criteria.significant_clusters_sinks,hb_bond_criteria.significant_clusters_sinks_cutoff,hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only]
for criteria in criterias:
 #file_handle=open('results_aucscores/AccuracyMetrics_%s_IEonly.txt' %(criteria),'w')
 file_handle_auc=open('results_aucscores/AUCScores_%s_IEonly.txt' %(criteria),'w')
 fields=['interaction_energy']
 thresholds=[33,33.5,34,34.5,35,35.5,36,36.5,37,37.5,38.0]
 off=0.10
 field='interaction_energy'
 for pick_threshold in thresholds:
    #file_handle.write('%s off threshold %f %f\n' %(field,off,pick_threshold))
    print('%s off threshold %f %f\n' %(field,off,pick_threshold))
    print_tprs_and_fprs(criteria,fields=[field],off=off,pick_threshold=pick_threshold,file_handle_auc=file_handle_auc)
  #thresholds=[0.05,0.10,0.15]
  #off=0.10
  #for pick_threshold in thresholds:
  # field='cluster_size'
  # file_handle.write('%s off threshold %f %f\n' %(field,off,pick_threshold))
  # print_tprs_and_fprs(criteria,heatmap_type=heatmap.size,off=off,pick_threshold=pick_threshold,file_handle=file_handle,file_handle_auc=file_handle_auc)
 #file_handle.close()
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

criterias=[hb_bond_criteria.significant_clusters_sinks]
for criteria in criterias:
 field='cluster_size'
 file_handle_auc=open('results_aucscores/AUCScores_%s_fprs_%s.txt' %(criteria,field),'w')
 fprs=[0.143,0.14,0.13,0.15,0.16,0.12,0.11,0.125]
 off=0.10
 for pick_fpr in fprs:
  print_tprs_and_fprs(criteria,heatmap_type=heatmap.size,off=off,pick_fpr=pick_fpr,file_handle_auc=file_handle_auc)
 file_handle_auc.close()


#write_arrays_for_criteria(criteria)
#criteria = hb_bond_criteria.significant_clusters_sinks_cutoff
#write_arrays_for_criteria(criteria)
#criteria = hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only
#write_arrays_for_criteria(criteria)
#get_heatmaps_for_criteria_cutoff(criteria=criteria,heatmap_type=heatmap.size)
#criteria = hb_bond_criteria.significant_clusters_sinks
#print_tprs_and_fprs(criteria,off=0.55)
#write_arrays_for_criteria(criteria)

