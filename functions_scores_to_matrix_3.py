import os
import sys
import numpy as np
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
from datetime import date
large_num = 30.0

#from functions_parsefiles import *

from functions_filterpdbforparams import *


listnotfound = []
'''
inputsfile = "inputs_generate_matrices.py"
print "Executing file ", inputsfile
execfile(inputsfile)
'''
listres_3letter = ["LEU","VAL","ALA","ILE","MET","PRO","GLY","PHE","TRP","TYR","SER","THR","ASN","GLN","HIS","LYS","ARG","ASP","GLU"]
listres = ['L','V','A','I','M','P','G','F','W','Y','S','T','N','Q','H','K','R','D','E']
listres1 = listres
listres2 = listres
assert(len(listres)==19)
N= len(listres)
data = dict()
columns = dict()
columns = {'interaction_score':7,'tag':14,'exp_binary':8,'true_signal':20}
datatype=[ ]
datalabels = ['interaction_energy', 'distance_catalysis', 'substrate_ca_no_super_rmsd', 'substrate_sequon_ca_no_super_rmsd', 'substrate_ca_rmsd', 'kdescore','phi_500','psi_500','rmsddih','rmsddih_xp1','description','HWUO1B-THR7N','HWUO1-THR7N','UDP2OPB-THR7N', 'UDP3OPB-THR7N','yashes','true_value']
reslist1 = ["LEU","VAL","ALA","ILE","MET","PRO","GLY","PHE","TRP","TYR","SER","THR","ASN","GLN","HIS","LYS","ARG","ASP","GLQ"]
reslistnames1 = listres
reslist2 = reslist1
reslistnames2 = reslistnames1

def getfilenames_unglycosylated():

  info ={}
  info["N"] = [10]
  peptidetemplate = "AXTXAPRC_498%s_500"
  info["files"]=[]
  info["files_filtered"]=[]
  info["sortfield"] = ["interaction_energy"]

  info["keys"]=[]
  outputtype = "Unglycosylated_fw_lr_hr"
  filetemplate = "%s%s/score_only_%s.sc"
  filetemplate_filtered = "AdditionalProps_%s_%s_%s_Top%03d.txt"

  info['outtag'] = 'July17_%s_RMSD_Distance_State1' %(outputtype)
  info['statfile'] = "Statistics/Stats_%s" %info["outtag"] + "_%s_Top%03d.txt"

  info['outfile_list'] = []
  tagtemplate = info['outtag'] + "_%s_Top%03d"
  tagtemplate2 = info['outtag'] + "_Top%03d"
  info['filelabels']=[]
  for num in info['N']:
    for sf in info['sortfield']:
      for k,res1 in enumerate(reslist1):
        peptidename = peptidetemplate %res1
        for i,res2 in enumerate(reslist2):
          info["files"].append(filetemplate %(peptidename,res2,outputtype))
          #print filetemplate %(peptidename,res2,outputtype)
          curkey = "A%sT%sAPRC" %(reslistnames1[k],reslistnames2[i])
          info["keys"].append(curkey)
          a = filetemplate_filtered %(curkey,info['outtag'],sf,num)
          b = tagtemplate %(sf,num)
          c = tagtemplate2 %(num)
          d = sf
          info['files_filtered'].append((a,b,c,d))
          info['filelabels'].append(b)
  return info['files_filtered'], info["keys"],info['filelabels']

def getfilenames_unglycosylated_State2():
  info ={}
  info["N"] = [10]
  peptidetemplate = "AXTXAPRC_498%s_500"
  info["files"]=[]
  info["files_filtered"]=[]
  info["sortfield"] = ["interaction_energy"]

  outputtype = "Unglycosylated_fw_lr_hr"
  filetemplate = "%s%s/score_only_%s.sc"
  filetemplate_filtered = "AdditionalProps_%s_%s_%s_Top%03d.txt"

  info['filelabels']=[]
  info['keys']=[]
  info['outtag'] = 'July29_%s_RMSD_Distance_State2' %(outputtype)

  info['statfile'] = "Statistics/Stats_%s" %info["outtag"] + "_%s_Top%03d.txt"
  info['outfile_list'] = []
  tagtemplate = info['outtag'] + "_%s_Top%03d"
  tagtemplate2 = info['outtag'] + "_Top%03d"
  info['filelabels']=[]
  for num in info['N']:
    for sf in info['sortfield']:
      for k,res1 in enumerate(reslist1):
        peptidename = peptidetemplate %res1
        for i,res2 in enumerate(reslist2):
          info["files"].append(filetemplate %(peptidename,res2,outputtype))
          #print filetemplate %(peptidename,res2,outputtype)
          curkey = "A%sT%sAPRC" %(reslistnames1[k],reslistnames2[i])
          info["keys"].append(curkey)
          a = filetemplate_filtered %(curkey,info['outtag'],sf,num)
          b = tagtemplate %(sf,num)
          c = tagtemplate2 %(num)
          d = sf
          info['files_filtered'].append((a,b,c,d))
          info['filelabels'].append(b)
  return info['files_filtered'], info["keys"],info['filelabels']

def parsecsvfile(infile):
  f = open(infile,'r')
  lines = f.readlines()
  f.close()
  for line in lines[1:362]:
   vals = line.split(',')

   isc = float(vals[columns['interaction_score']])
   signal = float(vals[columns['true_signal']])
   binval =  int(vals[columns['exp_binary']])
   tag = vals[columns['tag']][7:15]
   superdata['experimentalresults_2018']['interaction_energy'][tag]=binval
   superdata['yashes']['interaction_energy'][tag]=isc
   superdata['experimentalresults_2018']['true_signal'][tag]=signal
   if not tag in filelabeleddata['interaction_energy']:
                filelabeleddata['interaction_energy'][tag]={}
   filelabeleddata['interaction_energy'][tag]['yashes'] = isc
   filelabeleddata['interaction_energy'][tag]['experimentalresults_2018'] = binval
  return superdata, filelabeleddata

def makematrix(filelabel,allvals,curlabel,bFilelabeleddata=False,listres_1=listres,listres_2=listres,dataframe=True):
  i=0
  for xminusone in listres_1: #y
     j=0
     for xplusone in listres_2: #x
        key = 'A%sT%sAPRC' %(xminusone,xplusone)
        if bFilelabeleddata:
                if key in filelabeleddata[curlabel]:
                        allvals[i,j]=filelabeleddata[curlabel][key][filelabel]
        else:
                if key in superdata[filelabel][curlabel]:
                   allvals[i,j]=superdata[filelabel][curlabel][key]
        j+=1
     i+=1
  df = pd.DataFrame()
  if dataframe:
        rows = pd.Index( listres_1 )
        columns = pd.Index( listres_2  )
        df = pd.DataFrame( data = allvals , index=rows , columns = columns )
        #print df
  return df, allvals


def getkeys(listres_1=listres,listres_2=listres):
  keys=[]
  for xminusone in listres_1: #y
     for xplusone in listres_2: #x
        key = 'A%sT%sAPRC' %(xminusone,xplusone)
        keys.append(key)
  return keys

def makematrixfordict(dictdata,allvals,listres_1=listres,listres_2=listres,dataframe=True):
  i=0
  for xminusone in listres_1: #y
     j=0
     for xplusone in listres_2: #x
        key = 'A%sT%sAPRC' %(xminusone,xplusone)
        if key in dictdata:
          allvals[i,j]=dictdata[key]
        j+=1
     i+=1
  df = pd.DataFrame()
  if dataframe:
        rows = pd.Index( listres_1 )
        columns = pd.Index( listres_2  )
        df = pd.DataFrame( data = allvals , index=rows , columns = columns )
        #print df
  return df, allvals

def makedataframefordict(dictdata,prop,listres_1=listres,listres_2=listres):
  flatarray = []
  minus = []
  plus = []
  i=0
  for xminusone in listres_1: #y
     j=0
     for xplusone in listres_2: #x
        key = 'A%sT%sAPRC' %(xminusone,xplusone)
        #minus.append(xminusone)
        #plus.append(xplusone)
        if key in dictdata:
          minus.append(xminusone)
          plus.append(xplusone)
          flatarray.append(dictdata[key])
        #else:
        #  flatarray.append(float('nan'))
        j+=1
     i+=1
  df = pd.DataFrame()
  df['-1']=minus
  df['+1']=plus
  df[prop]=flatarray
  return df

def getannotations(listres_1=listres,listres_2=listres):
  annotations=np.chararray((len(listres_1),len(listres_2)),itemsize=3)
  for i,xminusone in enumerate(listres_1):
     for j,xplusone in enumerate(listres_2):
        key = '%sT%s' %(xminusone,xplusone)
        annotations[i,j]=key
        print( i,j,key,annotations[i,j])
  return annotations


def addlabeleddata(peptidetag,data,label):

   for curlabel in data:
        if not peptidetag in filelabeleddata[curlabel]:
                filelabeleddata[curlabel][peptidetag]={}
        filelabeleddata[curlabel][peptidetag][label] = float(data[curlabel])

def addlabeleddatalists(peptidetag,data,label):

   for curlabel in data:
        if not peptidetag in filelabeleddata[curlabel]:
                filelabeleddata[curlabel][peptidetag]={}
        filelabeleddata[curlabel][peptidetag][label] = data[curlabel]

def get_max_value(input_matrix, discard_value):
        max_value = -100
        N1 = input_matrix.shape[0]
        N2 = input_matrix.shape[1]
        for i in range(0,N1):
                for j in range(0,N2):
                        current_value = input_matrix[i,j]
                        if current_value != discard_value and current_value >= max_value:
                                max_value = current_value
                                i_max = i
                                j_max = j
        return max_value, i_max, j_max

def get_min_value(input_matrix, discard_value):
        min_value = 100
        N1 = input_matrix.shape[0]
        N2 = input_matrix.shape[1]
        for i in range(0,N1):
                for j in range(0,N2):
                        current_value = input_matrix[i,j]
                        if current_value != discard_value and current_value <= min_value:
                                min_value = current_value
        return min_value

def getweightedaverage(vals,weights):
        if len(vals) != len(weights):
                raise MyError('cannot calculate weighted average - weights and vals dont match')
        else:
                wavg = 0.0
                for v,w in zip(vals,weights):
                        wavg += v * w
                wavg /= float(sum(weights))
                return wavg

def getcutoffaverage(vals,cutoffvals,cutoff,default,verbose=False):
        avg = 0.0
        counts = 0.0
        for v,cv in zip(vals,cutoffvals):
                if cv <= cutoff:
                     if verbose:
                        print( v , cv)
                     avg += v
                     counts += 1.0
        if counts > 0.0:
                avg /= float(counts)
                return avg
        else:
                return default


def getcutoffaverage_multiple(vals,cutoffvalsdict,cutoffdict,default,verbose=False):
        avg = 0.0
        counts = 0.0
        for i , v in enumerate(vals):
                bcutoff = False
                for prop in cutoffdict:
                        bcutoff = False
                        if cutoffvalsdict[prop] <= cutoffdict[prop]:
                                if verbose:
                                        print(i , v , prop, cutoffvalsdict[prop])
                                bcutoff = True
                if bcutoff:
                     avg += v
                     counts += 1.0
        if counts > 0.0:
                avg /= float(counts)
                return avg
        else:
                return default


def makematrixfromlist(filelabel,allvals,curlabel,default,calcparams={},bFilelabeleddata=True,listres_1=listres,listres_2=listres,dataframe=True):
  verbose = False
  for i,xminusone in enumerate(listres_1): #y
     for j,xplusone in enumerate(listres_2): #x
        key = 'A%sT%sAPRC' %(xminusone,xplusone)
        if key in filelabeleddata[curlabel]:
                                calclist  = filelabeleddata[curlabel][key][filelabel]
                                Nvals = len(calclist)
                                if 'N' in calcparams:
                                        Nvals = min(calcparams['N'],Nvals)
                                calctype = 'avg'
                                if 'calctype' in calcparams:
                                        calctype = calcparams['calctype']
                                if Nvals<1:
                                        continue
                                        #print "no vals - do nothing - everything initialized to default"
                                elif calctype == 'avg':
                                        allvals[i,j] = sum (calclist[:Nvals]) / float(Nvals)
                                elif calctype == 'minscore':
                                        allvals[i,j] = calclist[0] #sorted vals by score
                                elif calctype == 'weightedavg': #weighting by interaction score
                                        weightlist = filelabeleddata['interaction_energy'][key][filelabel]
                                        allvals[i,j] = getweightedaverage(calclist[:Nvals],weightlist[:Nvals])
                                elif calctype == "cutoffavg":
                                        if key == 'ASTPAPRC':
                                                verbose = True
                                        else:
                                                verbose = False
                                        cutoffvals = filelabeleddata[calcparams['cutoff']['label']][key][filelabel]
                                        cutoff = calcparams['cutoff']['val']
                                        allvals[i,j] = getcutoffaverage(calclist[:Nvals],cutoffvals[:Nvals],cutoff,default,verbose)
                                elif calctype == 'cutoff_multiple':
                                        cutoffvalsdict = dict()
                                        bConsider = False
                                        for entry in calcparams['cutoff_multiple']:

                                                cutoffvals = filelabeleddata[entry][key][filelabel]
                                                avgcutoff =  sum (cutoffvals) / float(Nvals)
                                                bConsider =  avgcutoff <= calcparams['cutoff_multiple'][entry]
                                        if bConsider:
                                                allvals[i,j] = sum (calclist[:Nvals]) / float(Nvals)
                                        #getcutoffaverage_multiple(calclist,cutoffvalsdict,calcparams['cutoff_multiple'],default)
                                else:
                                        print("Not implemented")
  df = pd.DataFrame()
  if dataframe:
        rows = pd.Index( listres_1 )
        columns = pd.Index( listres_2  )
        df = pd.DataFrame( data = allvals , index=rows , columns = columns )
        #print df
  return df, allvals

def populateprefiltereddata_files(curfile,fields,key,uniquelabel):
                        if not os.path.exists(curfile):
                                listnotfound.append(curfile)
                                print("Not found: ",curfile)
                                return
                        f = open(curfile,'r')
                        lines = f.readlines()
                        f.close()

                        tempdata = dict()
                        dictdata = dict()
                        labeleddata, fields = readfiltereddatafromlines(lines,fields,dictdata) #no min max

                        addlabeleddatalists(key,labeleddata,uniquelabel) #lists instead of single avg values

def unglycosylated_MC_cutoff_prefiltereddata():
  default = large_num
  calctype_cutoff = { 'label':'HWUO1-THR7N', 'val':5.0}
  fields_mc = {"interaction_energy":[-50,0,0], 'HWUO1-THR7N':[4.7,5.7,0]}
  calcparams = {'calctype': 'cutoffavg', 'cutoff':calctype_cutoff}
  files_mc,keys,labels = getfilenames_unglycosylated()
  for i,entry in enumerate(files_mc):
      key = keys[i]
      populateprefiltereddata_files(entry[0],fields_mc,key,entry[1])
      matrix_1_mc = np.full((N,N),default)
      df_1_mc,matrix_1_mc = makematrixfromlist(labels[0],matrix_1_mc,"interaction_energy",default,calcparams)
  return df_1_mc, matrix_1_mc,superdata

def unglycosylated_State2_prefiltereddata():
  default = large_num
  calcparams = {'calctype': 'avg'}
  fields_s2 = {"interaction_energy":[-50,0,0]}
  files_s2,keys,labels = getfilenames_unglycosylated_State2()
  for i,entry in enumerate(files_s2):
      key = keys[i]
      populateprefiltereddata_files(entry[0],fields_s2,key,entry[1])
      matrix_1_s2 = np.full((N,N),default) #fill with large default values - good when some boxes dont have values yet
      df_1_s2,matrix_1_s2 = makematrixfromlist(entry[1],matrix_1_s2,"interaction_energy",default,calcparams) 
  return df_1_s2, matrix_1_s2,superdata

def unglycosylated_MC_cutoff_DeltaS2_prefiltered(matrix_1_mc, matrix_1_s2):
    
    #_,matrix_1_mc, _ = unglycosylated_MC_cutoff_prefiltereddata()
    #_,matrix_1_s2, _ = unglycosylated_State2_prefiltereddata()
    iT = listres.index('T')
    jF = listres.index('F') 

    delta_min = matrix_1_mc[iT,jF] - matrix_1_s2[iT,jF]
    matrix_1_mc_deltacutoff = matrix_1_mc
   
    mc_cutoff,i,j = get_max_value(matrix_1_mc, 1.0)#np.amax(matrix_1_mc) #negative values - min is max
    print( "mc min ", mc_cutoff, listres[i], listres[j])
    s2_cutoff, i, j = get_max_value(matrix_1_s2, 1.0)
    print(  "s2 min ",s2_cutoff, listres[i], listres[j])
    mc_cutoff = min(mc_cutoff,-20)
    s2_cutoff = mc_cutoff

    
    for i in range(0,matrix_1_mc.shape[0]):
      for j in range(0,matrix_1_mc.shape[1]):
        delta = matrix_1_mc[i,j] - matrix_1_s2[i,j]
        if matrix_1_mc[i,j] < mc_cutoff and matrix_1_s2[i,j] < s2_cutoff:
          if delta >= delta_min:
            matrix_1_mc_deltacutoff[i,j]=0.0
        else:
          matrix_1_mc_deltacutoff[i,j]=6.0 #so that binning works with numpy histogram
    return matrix_1_mc_deltacutoff

def setup():
  global fig
  global axes
  global superdata
  global filelabeleddata
  global statematrix
  global figRoc
  global axesRoc

  superdata = dict()
  filelabeleddata = dict() #Use when using 2 types of files for the same matrix - like unglycosylated and glycosylated - safer option always
  for curlabel in datalabels:
                        filelabeleddata[curlabel]={}
  statematrix = {}
  superdata['experimentalresults_2018']={'interaction_energy':{},'true_signal':{}}
  superdata['yashes']={'interaction_energy':{}}

