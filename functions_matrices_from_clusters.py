import os
import glob
import sys
import seaborn as sns
from MultivariablePlotting import *
from DataSortingFiltering import getdata_xyzc
from MultivariablePlottingGlobal import *
logging.basicConfig(filename='scatterplot.log',level=logging.DEBUG)
from MatrixPlotting import *
from functions_scores_to_matrix_3 import makematrixfordict,getkeys
from enum import Enum, unique , auto
from datetime import date
import statistics
from functions_matrices_from_clusters_criteria import *

d = date.today()
datestring = d.isoformat()
print( datestring)
hb_cutoff = 4.0
rmsd_cutoff = 1.0
topN = 200
top_average = 1

def get_default_value(value_name):
  if value_name == 'interaction_energy': return 99.0
  if value_name == 'distance_catalysis_HWUO1B-THR7N': return 99.0
  if value_name == 'distance_catalysis_HWUO1-THR7N': return 99.0
  if value_name == 'cluster_size' or value_name == 'cluster_size_filtered': return -1.0
  if value_name == 'substrate_ca_no_super_rmsd': return 99.0
  if value_name == 'substrate_sequon_ca_no_super_rmsd': return 99.0
  if value_name == 'distance_catalysis': return 99.0
  if value_name == 'experimental': return -1.0
  if value_name == 'interaction_energy_GTX': return 99.0
  if value_name == 'interaction_energy_TTY': return 99.0
  if value_name == 'sc_shapecomplementarity_sd': return 99.0
  if value_name == 'sc_shapecomplementarity': return -1
  if value_name == 'cluster_size_GTX' or value_name == 'cluster_size_TTY': return -1.0
  #if value_name == 'sc_shapecomplementarity_mean': return -1
  #if value_name == 'sc_shapecomplementarity_median': return -1
  #if value_name == 'sc_shapecomplementarity_quartile_lower': return -1
  #if value_name == 'sc_shapecomplementarity_quartile_upper': return -1
  if value_name == 'fa_atr': return 99.0
  if value_name == 'compare': return 1 #This is a special case
  

def cmap_for_field(value_name):
  if value_name == 'interaction_energy': return 'Purples_r'
  if value_name == 'distance_catalysis_HWUO1B-THR7N': return 'Blues_r'
  if value_name == 'distance_catalysis_HWUO1-THR7N': return 'Blues_r'
  if value_name == 'cluster_size' or value_name == 'cluster_size_filtered': return 'gist_yarg'
  if value_name == 'substrate_ca_no_super_rmsd': return 'RdPu_r'
  if value_name == 'substrate_sequon_ca_no_super_rmsd': return 'BuPu_r'
  if value_name == 'distance_catalysis': return 'Greens_r'
  if value_name == 'experimental': return 'Purples'
  if value_name == 'compare': return None # 'Paired'#'RdBu'
  if value_name == 'sc_shapecomplementarity': return 'Greys'
  if value_name == 'sc_shapecomplementarity_sd': return 'Greys_r'
  if value_name == 'fnat': return 'YlOrBr'
  if value_name == 'fnonnat': return 'YlOrBr'
  
  return 'Greys_r'

def color_for_field(value_name):
  if value_name == 'interaction_energy': return 'indigo'
  if value_name == 'distance_catalysis_HWUO1B-THR7N': return '#2976bb'
  if value_name == 'distance_catalysis_HWUO1-THR7N': return '#2976bb'
  if value_name == 'cluster_size' or value_name == 'cluster_size_filtered': return '#738595'
  if value_name == 'substrate_ca_no_super_rmsd': return 'mediumvioletred'
  if value_name == 'substrate_sequon_ca_no_super_rmsd': return 'turquoise'
  if value_name == 'distance_catalysis': return 'green'
  if value_name == 'experimental': return 'purple'
  if value_name.find('sc_shapecomplementarity') != -1: return 'green'
  if value_name == 'fa_atr': return 'maroon'
  return None

class heatmap(Enum):
  size = auto()
  value = auto()
  delta = auto()
  exp = auto()
  size_filtered = auto()
  compare = auto()


def get_colorbar_label(criteria):
  if criteria==hb_bond_criteria.significant_clusters_cutoff:
    label = 'd_hb < %2.1f' %hb_cutoff
    return label

def get_suffix(name_suffix,replace_this='.' , replace_by='-'):
  return name_suffix.replace(replace_this, replace_by)

def get_files(pattern,size=0.500):
  import glob
  raw_files = glob.glob(pattern)
  #print(len(raw_files))
  files=[]
  for rfile in raw_files:
      if os.path.exists(rfile):
        s = os.path.getsize(rfile) #s in bytes
        #print(rfile,s/(1000*1000.0),"mb")
        if s/(1000.0*1000.0) >size:
          files.append(rfile)
  #print( "Files:",'\n'.join(files))
  #print("TOTALFILES: ", len(files))
  return files

def get_df_indexlist_for_clus(df,column_name,labels,clus):
  indexlist = []
  for j,label in enumerate(labels):
    if label==clus: #Not just core_samples - "funnel" not be a core_sample and ( j in core_samples_):
      indexlist.append(j)
  return indexlist

def get_sink(df,column_name,labels,clus,type_value='sinks'):
  dfvalues = df[column_name].values
  indexlist = get_df_indexlist_for_clus(df,column_name,labels,clus)
  valuelist = [dfvalues[t] for t in indexlist]

  if type_value =='sinks_TopKCal':
    min_val_ie = df[labels_df['interaction_energy']][indexlist[0]] + 2.0*0.593
    temp_ = []
    for ik,index in enumerate(indexlist):
      cur_val_ie =  df[labels_df['interaction_energy']][index]
      if cur_val_ie <= min_val_ie:
        temp_.append(valuelist[ik])
        return np.mean(temp_)
  else:
    temparray = np.array(list(valuelist[:top_average]))
    return temparray.mean()


def filter_clus_by_criteria(df,columns,labels,clus,criteria,type_value='sinks'):

  [column_name1,column_name2] = columns
  print ('filter_clus_by_criteria',columns,column_name1,column_name2)
  if column_name1 is None:
    dfvalues1 = None
    value1 = None
  else:
    dfvalues1 = df[column_name1].values
    indexlist = get_df_indexlist_for_clus(df,column_name1,labels,clus)
  
  if column_name2 is None:
    dfvalues2 = None
    value2 = None
  else:
    dfvalues2 = df[column_name2].values
    indexlist = get_df_indexlist_for_clus(df,column_name2,labels,clus)
  
  #valuelist = [dfvalues[t] for t in indexlist]

  #print('Values of sinks', valuelist[0] , df[labels_df['interaction_energy']][indexlist[0]],df[labels_df['interaction_energy']][indexlist[1]] )

  count=0.0
  for index in indexlist:
    if not dfvalues1 is None:
      value1 = dfvalues1[index]
    if not dfvalues2 is None:
      value2 = dfvalues2[index]

    if check_criteria_decoy([value1,value2],criteria=criteria):
      print("filtered",value1,value2)
      count+=1.0

  print(count)
  return count

def get_sinks_for_columns(df,labels=None,type_value='sinks', columns=None):
  if columns is None:
    columns = [info['xfield'],info['yfield'],info['zfield'],info['z3field']] 

  if labels is None:
    labels = [ 0 for _ in range(0,df.shape[0]) ] #all belong to same cluster -since no clustering for fixed backbone case
  unique_labels = set(labels)
  for elem in unique_labels:
      if elem<0:
        unique_labels.remove(elem)
        break
  sinkarray = np.zeros((len(unique_labels) , len(columns)))
  for i,clus in enumerate((unique_labels)):
    if clus<0:
      continue
    for j,column in enumerate(columns):
      column_name = labels_df[column]
      sinkarray[clus,j] =  get_sink(df,column_name,labels,clus,type_value=type_value)
  return sinkarray


def get_matrix_for_dict(dictmat,value_name):
    default_value = get_default_value(value_name)
    mat = np.full((19,19), default_value ,dtype=float)
    df,mat =  makematrixfordict(dictmat,mat)
    return mat

def get_roc_auc_for_dict(dictmat,value_name,off=0.00):
    fpr = [0.0,1.0]
    tpr = fpr
    from compare_simulations_to_experiments import getexperimentaldatafrompickle,pfile_exp
    import sklearn.preprocessing as preprocess
    x_glyc, y_expdata,y_expdata_onoff = getexperimentaldatafrompickle(pfile_exp)
    data_y = np.ravel(y_expdata)
    #data_y_onoff = np.ravel(y_expdata_onoff)
    data_y_onoff = (preprocess.binarize(data_y.reshape(-1,1),threshold=off)).reshape(data_y.shape)
    filename = "pic.png"

    default_value = get_default_value(value_name)
    mat = np.full((19,19), default_value ,dtype=float)
    df,mat =  makematrixfordict(dictmat,mat)
    if value_name =="cluster_size" or value_name=="sc_shapecomplementarity":
      data_x =  np.ravel(mat)
    else:
      data_x =  -1.0 * np.ravel(mat)
    import functions_roc_calculate as roccalc
    return roccalc.simplistic_roc_curve(data_y_onoff,data_x)

def get_minmax_matrix(mat,default_value):
    vmax = round(np.amax(mat[mat != default_value ]))
    vmin = round(np.amin(mat[mat != default_value ]))

    for i in range(19):
      for j in range(19):
        if mat[i,j]==default:
          mat[i,j]=0.0
        else:
          mat[i,j] = round((mat[i,j]-vmin)/(vmax-vmin))

    return mat

def plot_roc_auc_for_dict(dictmat,name_suffix=None,show=False,title=None,metadata=None,value_name='interaction_energy',off=0.0,clus2d=False):
  str_score = ""
  fpr,tpr,thresholds,score = get_roc_auc_for_dict(dictmat,value_name,off)

  filename = 'testroc.png'
  basename = "results/rocs_off%f_ta%d" %(off,top_average)
  if not name_suffix is None:
      if not clus2d:
        basename = "results/rocs_off%f_ta%d" %(off,top_average)
        filename="%s/roccurve_%s_%s.png" %(basename,name_suffix,value_name)
      else:
        basename = "results_clus2d/rocs_off%f_ta%d" %(off,top_average)
        filename="%s/roccurve_%s_%s.png" %(basename,name_suffix,value_name)
  if not os.path.exists(basename):
    os.system('mkdir -p %s' %basename)

  str_score = "AUC = %4.3f" %score
  if not name_suffix is None:
    title = '%s\n%s\n%s OFF=%f' %(datestring,name_suffix,value_name,off) 
  import functions_roc_calculate as roccalc
  roccalc.plot_roc(fpr,tpr,filename=filename,title=title,str_score=str_score,show=show)

def get_positionwise_roc_auc_for_dict(dictmat,name_suffix=None,show=False,title=None,metadata=None,value_name='interaction_energy',off_list=[0.0,0.05,0.10,0.20,0.30],clus2d=False,positions=[-1,1],file_handle=open('temp.txt','w')):
 for off in off_list:
  for pos in positions:
    fpr = [0.0,1.0]
    tpr = fpr
    from compare_simulations_to_experiments import getexperimentaldatafrompickle,pfile_exp
    import sklearn.preprocessing as preprocess
    x_glyc, y_expdata,y_expdata_onoff = getexperimentaldatafrompickle(pfile_exp)

    default_value = get_default_value(value_name)
    mat = np.full((19,19), default_value ,dtype=float)
    df,mat =  makematrixfordict(dictmat,mat)
    mask = (mat == default_value)
    if pos==1:
      data_x = np.ravel(np.mean(mat,axis=0))
      data_y = np.ravel(np.mean(y_expdata,axis=0))
    if pos==-1:
      data_x = np.ravel(np.mean(mat,axis=1))
      data_y = np.ravel(np.mean(y_expdata,axis=1))

    data_y_onoff = (preprocess.binarize(data_y.reshape(-1,1),threshold=off)).reshape(data_y.shape)
    print(data_y_onoff,data_y_onoff.shape)
    nonzeros=np.count_nonzero(data_y_onoff)
    if nonzeros==0.0 or nonzeros==len(data_y_onoff):
      continue

    str_data_x = [str(v) for v in list(data_x)]
    vals = ','.join(str_data_x)    

    if value_name =="cluster_size" or value_name=="sc_shapecomplementarity":
      data_x =  data_x
    else:
      data_x =  -1.0 * data_x
    import functions_roc_calculate as roccalc
    fpr,tpr,thresholds,score= roccalc.simplistic_roc_curve(data_y_onoff,data_x)
    file_handle.write('%f,%s,%d,%f\n'%(off,value_name,pos,score))

def write_roc_auc_for_dict(dictmat,name_suffix=None,show=False,title=None,metadata=None,value_name='interaction_energy',off_list=[0.0,0.05,0.10,0.20,0.30,0.40,0.50,0.55],clus2d=False,file_handle=open('temp.txt','w')):
 for off in off_list:
    _,_,_,score = get_roc_auc_for_dict(dictmat,value_name,off=off)
    file_handle.write('%f,%s,%f\n'%(off,value_name,score))


def plot_roc_auc_for_dicts(dictmats,name_suffix=None,show=False,title=None,metadata=None,value_names=['interaction_energy'],off=0.00,clus2d=False):

  list_fpr = []
  list_tpr = []
  list_str_scores = []
  colors = []
  for i,dictmat in enumerate(dictmats):
    fpr,tpr,thresholds,score = get_roc_auc_for_dict(dictmat,value_names[i],off)
    list_fpr.append(fpr)
    list_tpr.append(tpr)
    list_str_scores.append("AUC = %4.3f" %score)
    colors.append(color_for_field(value_names[i]))
  
  value_name_combined = '-'.join(value_names)
  filename = 'testroc.png'
  basename = ''
  if not name_suffix is None:
      if not clus2d:
        basename = 'results/rocs_combined_off%f_ta%d' %(off,top_average)
        filename="%s/roccurve_%s_%s.png" %(basename,name_suffix,value_name_combined)
      else:
        basename = 'results_clus2d/rocs_combined_off%f_ta%d' %(off,top_average)
        filename="%s/roccurve_%s_%s.png" %(basename,name_suffix,value_name_combined)
  
  if not os.path.exists(basename):
    os.system('mkdir -p %s' %basename)
  if not name_suffix is None:
      title = '%s\n%s\n%s' %(datestring,name_suffix,value_name_combined)
  import functions_roc_calculate as roccalc
  roccalc.plot_rocs(list_fpr,list_tpr,filename=filename,title=title,str_scores=list_str_scores,colors=colors)

def get_false_positives_and_false_negatives(dictmat,field_name,off=0.10,pick_tpr=0.99,pick_fpr=None,pick_threshold=None,verbose=True,dict_roc={'tp':1.0,'tn':0.7,'fp':0.3,'fn':0.0},file_handle_auc=open('aucscores.txt','w'),outfile='out.out'):

  fpr,tpr,thresholds,score = get_roc_auc_for_dict(dictmat,field_name,off)
  file_handle_auc.write('%f,%s,%f\n'%(off,field_name,score)) 
  if verbose:
    for i,ifpr in enumerate(fpr):
      print(ifpr,tpr[i],thresholds[i])
    print(score)
  from compare_simulations_to_experiments import getexperimentaldatafrompickle_dict,pfile_exp
  import sklearn.preprocessing as preprocess
  expdict = getexperimentaldatafrompickle_dict(pfile_exp)
  def false_positive(d_,key,off):
    if d_[key] < off:
      return True
    else:
      return False

  def false_negative(d_,key,off):
    if d_[key] >= off:
      return True
    else:
      return False
  
  allkeys=getkeys()
  for key in allkeys:
    if not key in dictmat:
       dictmat[key]=get_default_value(field_name)

  sortedkeys = sorted(list(dictmat.keys()))
  print('Total keys ',len(sortedkeys))
  import numpy as np
  def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

  if (pick_threshold is None) and (pick_fpr is None):
   i_max = find_nearest(tpr,pick_tpr)
   threshold=thresholds[i_max]
   if verbose:
    print('i_max ',i_max)
    print('threshold ',threshold)
  if not pick_threshold is None:
   threshold=pick_threshold
  if not pick_fpr is None:
   i_max = find_nearest(fpr,pick_fpr)
   threshold=thresholds[i_max]
   if verbose:
    print('fpr i_max ',i_max)
    print('fpr threshold ',threshold)

  outf = open(outfile,'w')
  outf.write('threshold: %f\n'%threshold)

  if not (field_name =="cluster_size" or field_name=="sc_shapecomplementarity"):
      for key in dictmat:
         dictmat[key] = -1*dictmat[key]
  dictmat_compare = {}

  print('False Positives')
  outf.write('False Positives\n')
  keys_false_positives_size = []
  for key in sortedkeys:
    if dictmat[key]>=threshold and false_positive(expdict,key,off):
      #print(key,dictmat[key],expdict[key],'\n')
      outf.write('%s\t%f\t%f\n' %(key,dictmat[key],expdict[key]))
      keys_false_positives_size.append(key)
      dictmat_compare[key] = dict_roc['fp']
  print('# %d' %(len(keys_false_positives_size)))

  print('False Negatives')
  outf.write('False Negatives\n')
  keys_false_negatives_size = []
  for key in sortedkeys:
    if dictmat[key]<threshold and false_negative(expdict,key,off):
      #print(key,dictmat[key],threshold,expdict[key],'\n')
      outf.write('%s\t%f\t%f\n' %(key,dictmat[key],expdict[key]))
      keys_false_negatives_size.append(key)
      dictmat_compare[key] = dict_roc['fn']
  print('# %d' %(len(keys_false_negatives_size)))

  print('True Positives or True Negatives')
  outf.write('True Positives and Negatives\n')
  keys_true_positives_size = []
  keys_true_negatives_size = []
  for key in sortedkeys:
    if dictmat[key]>=threshold and (not false_positive(expdict,key,off)):
      #print(key,dictmat[key],threshold,expdict[key],'\n')
      outf.write('%s\t%f\t%f\n' %(key,dictmat[key],expdict[key]))
      keys_true_positives_size.append(key)
      dictmat_compare[key] = dict_roc['tp']
    if dictmat[key]<threshold and not false_negative(expdict,key,off):
      #print('tn ',key,dictmat[key],threshold,expdict[key],'\n')
      keys_true_negatives_size.append(key)
      dictmat_compare[key] = dict_roc['tn']

  for key in dictmat:
   if key not in dictmat_compare:
    if not false_negative(expdict,key,off): #make sure it is negative 
     dictmat_compare[key] = dict_roc['tn'] #by default everything is negative
     keys_true_negatives_size.append(key)
    else:
     print("Error out: positive sample missing - %s" %key)
     exit()
  print('# %d' %(len(keys_true_positives_size)))
  print('# %d' %(len(keys_true_negatives_size)))
  tp=len(keys_true_positives_size)
  tn=len(keys_true_negatives_size) #2 cases where the data is not entered correctly in the matrix
  fp=len(keys_false_positives_size)
  fn=len(keys_false_negatives_size)
  acc = float(tp+tn)/(tp+tn+fp+fn)
  prec=float(tp)/(tp+fp)
  rec=float(tp)/(tp+fn)
  tpr=float(tp)/(tp+fn)
  tnr=float(tn)/(tn+fp)
  bal_acc=0.5*(tpr+tnr)
  fpr=float(fp)/(tn+fp)
  infd=tpr+tnr-1.0
  outf.write('Acc/Prec/tpr/tnr/fpr/bal_acc/infd %f/%f/%f/%f/%f/%f/%f\n'%(acc,prec,tpr,tnr,fpr,bal_acc,infd))
  outf.write('TP/TN/FP/FN %d/%d/%d/%d\n'%(tp,tn,fp,fn))
  outf.close()
  return dictmat,dictmat_compare , keys_false_positives_size,keys_false_negatives_size , keys_true_positives_size,keys_true_negatives_size

def get_false_positives_and_false_negatives_for_criteria(criteria,heatmap_type,field_name,off=0.10,pick_tpr=0.99,pick_threshold=None,pick_fpr=None,outfile='out.out',file_handle_auc=open('aucscores.txt','w')):

  if heatmap_type==heatmap.value:
    dictmat = get_topN_sequon_unglycosylated( criteria=criteria , heatmap_type=heatmap_type , column_name = labels_df[field_name])
  if heatmap_type==heatmap.size:
     dictmat = get_topN_sequon_unglycosylated( criteria = criteria ,  heatmap_type =  heatmap_type)

  prefix= 'results_aucscores/FullAccuracyMetrics_%s_%s_off%f' %(criteria,field_name,off)
  suffix='tpr%f' %pick_tpr
  if not pick_threshold is None:
       suffix='threshold%f' %pick_threshold
  if not pick_fpr is None:
       suffix='fpr%f' %pick_fpr
  outfile=prefix+'_'+suffix+'.out'
   
  return get_false_positives_and_false_negatives(dictmat,field_name,off=off,pick_tpr=pick_tpr,outfile=outfile,pick_threshold=pick_threshold,pick_fpr=pick_fpr,file_handle_auc=file_handle_auc)

def categorize_exp_data(array):
  temp = array.ravel()
  ret_array = temp
  i=0
  for element in temp:
    if element >= 0.70:
      ret_array[i] = 1
    elif element <=0.10:
      ret_array[i] = 0
    else:
      ret_array[i]=2
    i+=1
  return ret_array.reshape(array.shape)

def feature_importance(listofdictmats,name_suffix=None,show=False,title=None,metadata=None,list_value_names=['interaction_energy']):
  name_suffix = '' 

  print(name_suffix)
  from compare_simulations_to_experiments import getexperimentaldatafrompickle,pfile_exp
  x_glyc, y_expdata,y_expdata_onoff = getexperimentaldatafrompickle(pfile_exp)
  y_expdata_cat = categorize_exp_data(y_expdata)
  data_y_cat = np.ravel(y_expdata_cat)#np.ravel(y_expdata)
  data_y_onoff = np.ravel(y_expdata_onoff)
  filename = "pic.png"

  combined_mat = np.full((361,len(list_value_names)),0)
  ivals = 0
  for dictmat,value_name in zip(listofdictmats,list_value_names):
    default_value = get_default_value(value_name)

    mat = np.full((19,19), default_value ,dtype=float)
    df,mat =  makematrixfordict(dictmat,mat)
    combined_mat[:,ivals] = mat.ravel()
    
    ivals += 1
  import functions_roc_calculate as roccalc
  #roccalc.feature_selection_trees(data_y_onoff,combined_mat)
  #roccalc.feature_selection_rfe(data_y_onoff,combined_mat)
  roccalc.decision_tree_regression(data_y_cat,combined_mat,list_value_names)

def mat_for_dict(dictmat, value_name):    
  default_value = get_default_value(value_name)
  mat = np.full((19,19), default_value ,dtype=float)
  df,mat =  makematrixfordict(dictmat,mat)
  return df, mat

def plot_heatmap_for_dict(dictmat,name_suffix=None,show=False,title=None,metadata=None,value_name='interaction_energy',vmax=None,vmin=None,label_colorbar = None,mask=None,annotate=False,clus2d=False,listres1=None,cbar_kws_fmt=None,fmt='0.1f'):

  default_value = get_default_value(value_name)
  if listres1 is None:
    mat = np.full((19,19), default_value ,dtype=float)
    df,mat =  makematrixfordict(dictmat,mat)
  else:
    mat = np.full((len(listres1),19), default_value ,dtype=float)
    df,mat =  makematrixfordict(dictmat,mat,listres_1=listres1)

  print(df)
  #print(value_name,'\n', mat)

  if mask is None:
    mask = (mat == default_value)

  if vmax is None:
    vmax = round(np.amax(mat[mat != default_value ]))

  if vmin is None:
    vmin = round(np.amin(mat[mat != default_value ]))

  #print(vmax,vmin)

  outfile=None
  basename = ''
  if not name_suffix is None:
    suffix = get_suffix(name_suffix)
    title = '%s\n%s%s\n' %(datestring,suffix,value_name)
    if not clus2d:
      if not annotate:
        basename = 'results/heatmaps_no_annotate_ta%d' %(top_average)
        outfile = '%s/%s_%s.png' %(basename,suffix,value_name)
      else:
        basename = 'results/heatmaps_ta%d' %top_average
        outfile = '%s/%s_%s.png' %(basename,suffix,value_name)
    else:
      if not annotate:
        basename = 'results_clus2d/heatmaps_no_annotate_ta%d' %top_average
        outfile = '%s/%s_%s.png' %(basename,suffix,value_name)
      else:
        basename = 'results_clus2d/heatmaps_ta%d' %top_average
        outfile = '%s/%s_%s.png' %(basename,suffix,value_name) 
  
  os.system('mkdir -p %s' %basename)
  cmap =  cmap_for_field(value_name)

  if value_name == 'compare':
    mask=None

  if listres1 is None:
    plotdataframesingle( df , cbar=True , title=title, cmap=cmap , mask=mask, vlimits=[float(vmin),float(vmax)], outfile = outfile,show=show,label_colorbar =label_colorbar ,annotate=annotate,cbar_kws_fmt=cbar_kws_fmt,fmt=fmt)
  else:
    plotdataframesingle( df , cbar=True , title=title, cmap=cmap , mask=mask, vlimits=[float(vmin),float(vmax)], outfile = outfile,show=show,label_colorbar =label_colorbar ,annotate=annotate,listres1=listres1,cbar_kws_fmt=cbar_kws_fmt,fmt=fmt)


def get_statistic_for_cluster(clus,column_name,labels,df,core_samples_,stat='mean'):
  dfvalues = df[column_name].values
  indexlist = get_df_indexlist_for_clus(df,column_name,labels,clus)
  ie = [dfvalues[j] for j in indexlist]
  if stat=='mean':
    if len(ie)>0:
      return statistics.mean(ie)

  return None

def get_largest_cluster(clusters,sizedict_normalized):
  maxsize = -1
  maxclus = -1
  for clus in clusters:
    if sizedict_normalized[clus]  > maxsize:
      maxsize = sizedict_normalized[clus]
      maxclus = clus
  return maxclus

def get_value_for_cluster(df,column_name,clus,centroids,labels,core_samples_,type_value='centroids',clus2d=False):
    fieldaxis = df.columns.get_loc(column_name) #can be different from criteria

    if type_value == 'centroids':
      if column_name != labels_df['interaction_energy']:
        if not clus2d:
          if fieldaxis < centroids.shape[1]:
            return  centroids[clus,fieldaxis]
        else:
          if fieldaxis < centroids.shape[1]:
            return  centroids[clus,fieldaxis] #NOT TESTED
      else:
        ie = get_statistic_for_cluster(clus,column_name,labels,df,core_samples_)
        if not ie is None:
          return ie
    if type_value == 'sinks' or type_value =='sinks_TopKCal':
        #print('get_sink ',column_name,clus)
        return get_sink(df,column_name,labels,clus,type_value=type_value)

    return None


def get_topN_sequon_glycosylated(residues=['all'],topN=20,heatmap_type=heatmap.value):
  pattern = 'PackandRelax_498*/glycosylated_hr/score.sc'
  files = get_files(pattern,size=0.110 )

  dictmat={}
  for curfile in files:
    print(curfile)

    df = get_filtered_combined_dataframe([curfile],filterTopN=True)
    sinks = get_sinks_for_columns(df,type_value='sinks',columns = [ 'interaction_energy'] )
    #print(sinks)
    dictmat[df['key'][0]] = sinks[0]
  
  return dictmat


def get_topN_sequon_glycosylated(column_name,residues=['all'],topN=20,heatmap_type=heatmap.value):
  pattern = 'PackandRelax_498*/glycosylated_hr_bk/score.sc'
  files = get_files(pattern,size=0.110 )

  dictmat={}
  for curfile in files:
    print(curfile)

    df = get_filtered_combined_dataframe([curfile],filterTopN=True)
    clus = 0
    labels = [ clus for _ in range(0,df.shape[0]) ]
    sink_value =  get_sink(df,column_name,labels,clus,type_value='sinks')
    print('sinks',sink_value)
    dictmat[df['key'][0]] = sink_value

  return dictmat


def get_topN_sequon_unglycosylated(pattern = 'PackandRelax_498*/unglycosylated/score.sc',residues=['all'],cluster=True,clusterOnly=True,topN=200,column_name=labels_df['distance_catalysis_HWUO1B-THR7N'],heatmap_type=heatmap.value,criteria=hb_bond_criteria.largest_cluster,cluster_size_sorted_id=0,usePickled=False,clus2d=False):

  if usePickled:
    pattern = 'PickleFiles/20191124/Un*.p'
  
  files = []
  if type(pattern) is list:
   for f in pattern:
    files += get_files(f)
  else:
    files = get_files(pattern )  

  #print(files)
  dictmat={}
  for curfile in files:
    #print(curfile)

    if not usePickled:
      df, centroids, sizedict, sizedict_normalized,labels,n_clusters,core_samples_ =  get_clusters_from_combined_dataframe([curfile],filterTopN=True,clus2d=clus2d,TopN=topN)
      #sinks must be a dictionary like centroids for similar usage
      type_value = get_type_value(criteria) 
      if clus2d:
        #for ignoring d_mc in this case
        columns= [info['xfield'],info['yfield'],info['zfield'],info['z3field']] 
        sinks = get_sinks_for_columns(df,labels,type_value=type_value)
      else:
        sinks = get_sinks_for_columns(df,labels,type_value=type_value)

    else:
      data_ = pickle.load(open(curfile,'rb'))
      df = data_['dataframe']
      centroids = data_['centroids']
      sizedict = data_['sizedict']
      sizedict_normalized = data_['sizedict_normalized']
      labels = data_['labels']
      n_clusters = data_['n_clusters']
      core_samples_ = data_['core_samples_']
      del data_

    fieldaxis = df.columns.get_loc(labels_df['distance_catalysis_HWUO1B-THR7N'])

    if criteria==hb_bond_criteria.significant_clusters_sinks_cutoff or criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd or criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff_rmsd or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only  or criteria==hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_rmsd or criteria==hb_bond_criteria.significant_clusters_sinks_GTX_TopKCal_cutoff_rmsd or criteria==hb_bond_criteria.significant_clusters_sinks_TTY_cutoff or criteria==hb_bond_criteria.significant_clusters_sinks_TTY_TopKCal_cutoff or criteria == hb_bond_criteria.significant_clusters_sinks_cutoff_sequonrmsd_only or criteria == hb_bond_criteria.significant_clusters_sinks_cutoff_sequonrmsd:
      clusters_criteria = check_hb_criteria(sinks,sizedict, sizedict_normalized, criteria=criteria,axis=fieldaxis)
    else:
      clusters_criteria = check_hb_criteria(centroids,sizedict, sizedict_normalized, criteria=criteria,axis=fieldaxis) #cluster_ids - meet criteria
    
    # WRONG sorted_cluster_keys = sorted(sizedict_normalized,reverse=True)
    sorted_d = dict([(k, sizedict_normalized[k]) for k in sorted(sizedict_normalized, key=sizedict_normalized.get, reverse=True)])
   
    if len(clusters_criteria)>0:

      if criteria==hb_bond_criteria.largest_cluster:

          maxclus = get_largest_cluster(clusters_criteria,sorted_d) 
          fieldaxis = df.columns.get_loc(column_name) #can be different from criteria

          if heatmap_type == heatmap.size:
            if maxclus != -1:
              dictmat[df['key'][0]] = sorted_d[maxclus]

          if heatmap_type == heatmap.value:
            value = get_value_for_cluster(df,column_name,maxclus,centroids,labels,core_samples_,clus2d=clus2d)
            if not value is None:
              dictmat[df['key'][0]] = value

      elif criteria==hb_bond_criteria.significant_clusters or criteria == hb_bond_criteria.significant_clusters_sinks or criteria ==  hb_bond_criteria.significant_clusters_sinks_TopKCal:
            type_value = get_type_value(criteria) 

            list_of_keys = list(sorted_d.keys())
            for curkey in sorted_d:
              if not curkey in clusters_criteria:
                list_of_keys.remove(curkey)

            if cluster_size_sorted_id < len(list_of_keys):

              clus_key = list_of_keys[cluster_size_sorted_id]

              if heatmap_type == heatmap.value:
                value = get_value_for_cluster(df,column_name,clus_key,centroids,labels,core_samples_,type_value=type_value,clus2d=clus2d)
                if not value is None:
                  dictmat[df['key'][0]] = value

              if heatmap_type == heatmap.size:
                dictmat[df['key'][0]] = sorted_d[clus_key]

      else:
      # criteria==hb_bond_criteria.significant_clusters_cutoff or criteria==hb_bond_criteria.significant_clusters_cutoff_rmsd or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd or  criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff_rmsd  or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff or  criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only:

            type_value = get_type_value(criteria) 
            #print("CRITERIA: ",criteria,type_value)
            list_of_keys = list(sorted_d.keys())
            total_size = 0 #sum up all the clusters
            sink_value = get_default_value('interaction_energy')#sinkiest of the sinks
            count = 0.0
            for curkey in sorted_d:
              if not curkey in clusters_criteria:
                continue

              if heatmap_type == heatmap.size:
                total_size  += sorted_d[curkey] #TODO: get number of decoys which satisfy criteria
              
              if heatmap_type == heatmap.size_filtered:

                if criteria == hb_bond_criteria.significant_clusters_cutoff or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff:
                  column_name = labels_df['distance_catalysis_HWUO1B-THR7N']
                  count += filter_clus_by_criteria(df,[column_name,None],labels,curkey,criteria=criteria,type_value='sinks')

                if criteria == hb_bond_criteria.significant_clusters_cutoff_rmsd_only or criteria == hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only:
                  column_name = labels_df['substrate_ca_no_super_rmsd']
                  count += filter_clus_by_criteria(df,[None,column_name],labels,curkey,criteria=criteria,type_value='sinks')

                if criteria == hb_bond_criteria.significant_clusters_cutoff_rmsd or criteria == hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd:
                  count += filter_clus_by_criteria(df,[ labels_df['distance_catalysis_HWUO1B-THR7N'], labels_df['substrate_ca_no_super_rmsd'] ],labels,curkey,criteria=criteria,type_value='sinks')

                if criteria == hb_bond_criteria.significant_clusters_cutoff_sequonrmsd_only or criteria == hb_bond_criteria.significant_clusters_sinks_cutoff_sequonrmsd_only:
                  column_name = labels_df['substrate_sequon_ca_no_super_rmsd']
                  count += filter_clus_by_criteria(df,[None,column_name],labels,curkey,criteria=criteria,type_value='sinks')
 
              if heatmap_type == heatmap.value:
                cluster_sink = get_value_for_cluster(df,labels_df['interaction_energy'],curkey,centroids,labels,core_samples_,type_value=type_value)
                if not cluster_sink is None:
                   if cluster_sink < sink_value:
                    value = get_value_for_cluster(df,column_name,curkey,centroids,labels,core_samples_,type_value=type_value)
                    


            if heatmap_type == heatmap.size:
              dictmat[df['key'][0]] = total_size            

            if heatmap_type == heatmap.size_filtered:
              fraction = count/float(df.shape[0])
              print('fraction',df['key'][0],count,df.shape[0],fraction)
              dictmat[df['key'][0]] = fraction

            if heatmap_type == heatmap.value:
              if not value is None:
                dictmat[df['key'][0]] = value
    else:
      if heatmap_type == heatmap.size:
        dictmat[df['key'][0]] = 0.0

  return dictmat

def heatmap_topN_sequon(dictmat,heatmap_type,criteria=None,column_name='',filterTopN=True,show=True,cluster_size_sorted_id=None,b_label_colorbar = False,annotate=True,roc=True,suffix=None,clus2d=False,off=0.0,listres1=None):

    print(listres1)
    name_suffix = '%s' %heatmap_type
    if not criteria is None:
      name_suffix += '_%s' %(criteria)
    if not suffix is None:
      name_suffix += '_%s' %suffix

    if filterTopN:
      name_suffix += '_TopN'

    if not cluster_size_sorted_id is None:    
      name_suffix += '_clus%d' %cluster_size_sorted_id

    clabel=None
    if b_label_colorbar:
      clabel = get_colorbar_label(criteria)

    print(name_suffix)
    if heatmap_type == heatmap.value:

      if column_name == labels_df['distance_catalysis_HWUO1B-THR7N']: #fieldaxis==1 -> hbond
        plot_heatmap_for_dict(dictmat,name_suffix, value_name = 'distance_catalysis_HWUO1B-THR7N',vmax=5.0,vmin=3.0, show=show,clus2d=clus2d,annotate=annotate,listres1=listres1)
        if roc:
          plot_roc_auc_for_dict(dictmat,name_suffix=name_suffix,value_name='distance_catalysis_HWUO1B-THR7N',clus2d=clus2d,off=off)

      if column_name == labels_df['substrate_ca_no_super_rmsd']:
        if criteria==hb_bond_criteria.significant_clusters_sinks or criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal or criteria==hb_bond_criteria.significant_clusters:
          vmax=2.5
        else:
          vmax=1.5
        plot_heatmap_for_dict(dictmat,name_suffix, value_name = 'substrate_ca_no_super_rmsd',vmax=vmax,vmin=0.5,show=show,clus2d=clus2d,annotate=annotate,listres1=listres1)
        if roc:
          plot_roc_auc_for_dict(dictmat,name_suffix=name_suffix,value_name='substrate_ca_no_super_rmsd',clus2d=clus2d,off=off)

      if column_name == labels_df['substrate_sequon_ca_no_super_rmsd']:
        plot_heatmap_for_dict(dictmat,name_suffix, value_name = 'substrate_sequon_ca_no_super_rmsd',vmax=1.0,vmin=0.5,show=show,clus2d=clus2d,annotate=annotate,listres1=listres1)
        if roc:
          plot_roc_auc_for_dict(dictmat,name_suffix=name_suffix,value_name='substrate_sequon_ca_no_super_rmsd',clus2d=clus2d,off=off)

      if column_name == labels_df['distance_catalysis']:
        plot_heatmap_for_dict(dictmat,name_suffix, value_name = 'distance_catalysis',vmax=3.0,vmin=4.5,show=show,clus2d=clus2d,annotate=annotate,listres1=listres1)
        if roc:
          plot_roc_auc_for_dict(dictmat,name_suffix=name_suffix,value_name='distance_catalysis',clus2d=clus2d,off=off)

      if column_name == labels_df['interaction_energy']:
        plot_heatmap_for_dict(dictmat,name_suffix, value_name = 'interaction_energy',vmax=-32,vmin=-40,show=show,clus2d=clus2d,annotate=annotate,listres1=listres1)
        if roc:
          plot_roc_auc_for_dict(dictmat,name_suffix=name_suffix,value_name='interaction_energy',clus2d=clus2d,off=off)

    if heatmap_type==heatmap.size:
      plot_heatmap_for_dict(dictmat,name_suffix,value_name='cluster_size',show=show,vmax=1.0,vmin=0.0,label_colorbar = clabel,clus2d=clus2d,annotate=annotate,listres1=listres1)
      if roc:
        plot_roc_auc_for_dict(dictmat,name_suffix=name_suffix,value_name='cluster_size',show=show,clus2d=clus2d,off=off)

    if heatmap_type==heatmap.size_filtered:
      plot_heatmap_for_dict(dictmat,name_suffix,value_name='cluster_size_filtered',show=show,vmax=1.0,vmin=0.0,label_colorbar = clabel,clus2d=clus2d,annotate=annotate,listres1=listres1)
      if roc:
        plot_roc_auc_for_dict(dictmat,name_suffix=name_suffix,value_name='cluster_size',show=show,clus2d=clus2d,off=off)
    
    if heatmap_type==heatmap.delta:
      plot_heatmap_for_dict(dictmat,name_suffix,value_name='interaction_energy',show=show,vmax=2.0,vmin=-2.0,label_colorbar = clabel,clus2d=clus2d,annotate=annotate,listres1=listres1)
      if roc:
        plot_roc_auc_for_dict(dictmat,name_suffix=name_suffix,value_name='delta',show=show,off=off)

    if heatmap_type==heatmap.exp:
      print("EXP")
      cbar_kws_fmt = '%.0f%%'
      plot_heatmap_for_dict(dictmat,name_suffix,value_name='experimental',show=show,vmax=100.0,vmin=0.0,label_colorbar = clabel,listres1=listres1,cbar_kws_fmt=cbar_kws_fmt,annotate=annotate,fmt='d')

    if heatmap_type==heatmap.compare:
      plot_heatmap_for_dict(dictmat,name_suffix,value_name='compare',show=show,vmax=1.0,vmin=0.0,label_colorbar = clabel,clus2d=clus2d,annotate=annotate,listres1=listres1)
