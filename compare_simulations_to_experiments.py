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

def plotexperimentaldatafrompickle(lineidentity):
    yexp = pickle.load(open(pfile_exp_array,'rb'))
    iT = S2M.listres.index('T')
    yexp_T = np.ravel(yexp[iT,:])
    print(yexp_T.shape)
    ann = np.chararray.ravel(S2M.getannotations()[iT,:]) #get annotations or threonine
    yexp_wt = pickle.load(open("experimentaldatamutants_"+lineidentity[0]+".p",'rb'))
    print ('new',yexp_wt)
    print ('old',yexp_T)
    #exit()
    fig,ax = plt.subplots(figsize=(8,8))
    fig.suptitle("Experimental Measurements Variation",fontsize=20)
    plt.plot([0,1],[0,1],color='dimgrey',linestyle='dashed')
    plt.scatter(yexp_T, yexp_wt,color='mediumvioletred',s=60,edgecolor='black',alpha=0.9)
    for i in range(0,len(yexp_T)):
      label = ann[i].decode()
      plt.annotate(label,(yexp_T[i],yexp_wt[i]),textcoords='offset points',xytext=(12,0),ha='center',fontsize=8)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.ylim(-0.05,1.05)
    plt.xlim(-0.05,1.05)
    axis_font = { 'size':20}
    plt.xlabel("Wildtype(Published 2018)",**axis_font)
    plt.ylabel("Wildtype(Measurement 2019)",**axis_font)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    plt.savefig("Exp2ExpCompare_PubvsRecent.png",transparent=True,dpi=fig.dpi)
    plt.show()
    for entry in lineidentity[1:]:
      yexp_more = pickle.load(open("experimentaldatamutants_"+entry+".p",'rb'))
      print(yexp_more.shape)
      fig,ax = plt.subplots(figsize=(8,8))
      fig.suptitle(entry,fontsize=25)
      plt.plot([0,1],[0,1],color='dimgrey',linestyle='dashed')
      plt.scatter(yexp_wt, yexp_more,color='purple',s=60,edgecolor='black',alpha=0.9)
      for i in range(0,len(yexp_T)):
        label = ann[i].decode()
        plt.annotate(label,(yexp_wt[i],yexp_more[i]),textcoords='offset points',xytext=(12,0),ha='center',fontsize=8)
      plt.xticks(fontsize=30)
      plt.yticks(fontsize=30)
      plt.ylim(-0.05,1.05)
      plt.xlim(-0.05,1.05)
      axis_font = { 'size':20}
      plt.xlabel("Wildtype",**axis_font)
      plt.ylabel(entry,**axis_font)
      plt.tight_layout()
      plt.subplots_adjust(top=0.9)
      plt.savefig("Exp2ExpCompare_"+ entry + ".png",transparent=True,dpi=fig.dpi)
      plt.show()   
    
  plotexperimentaldatafrompickle(lineidentity)

def getexperimentaldatafromfile(infile,serialize=True):
  S2M.setup()
  superdata,_ = S2M.parsecsvfile(infile)
  y_expdata_onoff = np.zeros((19,19))
  x_glyc = np.zeros((19,19))
  y_expdata_onoff_df, y_expdata_onoff = S2M.makematrixfordict(superdata['experimentalresults_2018']['interaction_energy'],y_expdata_onoff)
  x_glyc_df, x_glyc = S2M.makematrixfordict(superdata['yashes']['interaction_energy'],x_glyc)

  #print("Y:",y_expdata_onoff)#,"\n X:",x_glyc)
  y_expdata = np.zeros((19,19))
  y_expdata_df, y_expdata = S2M.makematrixfordict(superdata['experimentalresults_2018']['true_signal'],y_expdata)
  #serialize y_expdata
  print(y_expdata)
  f1 = open(npfile_exp_array,'wb')
  np.savetxt(f1,y_expdata)
  f1.close()
  if serialize:
    f1 = open(pfile_exp_array,'wb')
    pickle.dump(y_expdata,f1)
    f1.close()
    f2 = open(pfile_exp_onoff_array,'wb')
    pickle.dump(y_expdata_onoff, f2)
    f2.close()
    f3 = open(pfile_exp,'wb')
    pickle.dump(superdata,f3)
    f3.close()   
    f4 = open(pfile_glyc,'wb')
    pickle.dump(x_glyc,f4)
    f4.close()
  return x_glyc,y_expdata, y_expdata_onoff

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

def getexperimentaldatafrompickle_dict(pfile=pfile_exp,percent=False):
    superdata = {}
    superdata = pickle.load(open(pfile,'rb'))
    dictmat = superdata['experimentalresults_2018']['true_signal']
    if percent:
      dictmat_ = {}
      for key in dictmat:
        dictmat_[key] = dictmat[key]*100.0
      del dictmat
      dictmat = dictmat_
    return dictmat

def plotscattergeneric(x,y):
  fig = plt.figure(figsize=(8,8))
  plt.scatter(x, y,color='orchid',s=20)
  plt.xticks(fontsize=15)
  plt.yticks(fontsize=15)
  plt.ylim(-0.1,1.1)
  #plt.xlim(-15,-50)
  plt.show()

def plotscattergeneric3(x1,y1,x2,y2,x3,y3):
  fig = plt.figure(figsize=(8,8))
  plt.scatter(x1, y1,s=20,color='orange',label="Michaelis")
  plt.scatter(x2, y2,s=20,color='orchid',label='S2')
  plt.scatter(x3, y3,s=20,color='mediumseagreen',label='Glyc')
  plt.xticks(fontsize=15)
  plt.yticks(fontsize=15)
  plt.legend()
  #plt.ylim(-0.1,1.1)
  plt.ylim(0,-50)
  plt.show()


def plotcombinedroc(pfilesdict,filename):
  x_glyc, y_expdata,y_expdata_onoff = getexperimentaldatafrompickle(pfile_exp)
  data_y = np.ravel(y_expdata)
  data_y_onoff = np.ravel(y_expdata_onoff)
  import functions_roc_calculate as roccalc
  fig,ax = plt.subplots(figsize=(12,12))
  fig.suptitle(datestring,fontsize=25)
  plt.plot([0.0,1.0],[0.0,1.0],color='black',linestyle='dashed',linewidth=2)
  for key in pfilesdict:
     pfile = pfilesdict[key]
     x_ = getunglycosylatedfrompickle(pfile)
     data_x = -1.0 * np.ravel(x_)
     fpr,tpr,thresholds,score = roccalc.simplistic_roc_curve(data_y_onoff,data_x)
     str_score = "%s AUC = %4.3f" %(key,score)
     plt.plot(fpr,tpr,color=colors[key],linewidth=4,label=str_score) 
  #plt.text(0.6,0.1,str_score,size=40,color='black')
  plt.ylim(0.0,1.0)
  plt.xlim(0.0,1.0)
  axis_font = { 'size':35}
  plt.xlabel("False Positive Rate (FPR)",**axis_font)
  plt.ylabel("True Positive Rate (TPR)",**axis_font)
  plt.tick_params(axis='both',labelsize=35)
  plt.legend()
  plt.tight_layout()
  plt.subplots_adjust(top=0.85)
  plt.savefig(filename,transparent=True,dpi=fig.dpi)
  plt.show()


