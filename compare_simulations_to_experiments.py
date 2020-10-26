import numpy as np
import functions_scores_to_matrix_3 as S2M
import matplotlib.pyplot as plt; plt.style.use('classic')
import matplotlib.ticker as ticker
import pickle 
from datetime import date

large_num = 30.0
pfile_exp = "PickleFiles/experimentaldata_Jewett2018_yashes_dict.p"
pfile_exp_array = "PickleFiles/experimentaldata_19x19array.p"
npfile_exp_array = "PickleFiles/nptxt_experimentaldata_19x19array.txt"
pfile_exp_onoff_array = "PickleFiles/experimentaldata_onoff_19x19array.p"

d = date.today()
datestring = d.isoformat()
colors= {'H1-Glyc':'orangered','H2-Aglyc':'limegreen','H3-Aglyc':'royalblue'}

def get_experimental_data_for_sequon(key1,key2,format_aa = 'threeletter'):
  yexp = pickle.load(open(pfile_exp_array,'rb'))
  if format_aa=='threeletter':
    iM = S2M.listres_3letter.index(key1)
    iP = S2M.listres_3letter.index(key2)
  elif format_aa == 'oneletter':
    iM = S2M.listres.index(key1)
    iP = S2M.listres.index(key2)
  yexp_sequon = yexp[iM,iP]
  return yexp_sequon

def get_experimental_data_for_pos_residue(key1,pos,format_aa = 'oneletter'):
  if key1=='C' or key1=='CYS':
    return 0.0
  yexp = pickle.load(open(pfile_exp_array,'rb'))
  if format_aa=='threeletter':
    iM = S2M.listres_3letter.index(key1)
  elif format_aa == 'oneletter':
    iM = S2M.listres.index(key1)
  if pos==-1:
    yexp_sequon = yexp[iM,:]
  elif pos==1:
    yexp_sequon = yexp[:,iM]
  else:
    exit()
  y_ = yexp_sequon.ravel() 
  return np.mean(y_)

def get_experimental_data_postionwise_mean(aa_list,position_list=[-1,1]):
  pos_array = np.zeros((len(position_list),len(aa_list)))
  for ipos, pos in enumerate(position_list):
    for iaa,aa in enumerate(aa_list):
      pos_array[ipos,iaa] = get_experimental_data_for_pos_residue(aa,pos)
  return pos_array

def getexperimentaldatafrompickle(pfile=pfile_exp):
    superdata = {}
    superdata = pickle.load(open(pfile,'rb'))
    y_expdata_onoff = np.zeros((19,19))
    y_expdata = np.zeros((19,19)) 
    x_glyc = np.zeros((19,19))
    y_expdata_onoff_df, y_expdata_onoff = S2M.makematrixfordict(superdata['experimentalresults_2018']['interaction_energy'],y_expdata_onoff)
    x_glyc_df, x_glyc = S2M.makematrixfordict(superdata['yashes']['interaction_energy'],x_glyc)
    y_expdata = np.zeros((19,19))
    y_expdata_df, y_expdata = S2M.makematrixfordict(superdata['experimentalresults_2018']['true_signal'],y_expdata)
    return x_glyc, y_expdata,y_expdata_onoff

def getexperimentaldatafrompickle_dict(pfile=pfile_exp,percent=False,round_=True):
    superdata = {}
    superdata = pickle.load(open(pfile,'rb'))
    dictmat = superdata['experimentalresults_2018']['true_signal']
    if percent:
      dictmat_ = {}
      for key in dictmat:
        value = dictmat[key]*100.0
        if round_:
          dictmat_[key]=int(value)
          print(key,dictmat_[key])
        else:
          dictmat_[key]=value
      del dictmat
      dictmat = dictmat_
    return dictmat

