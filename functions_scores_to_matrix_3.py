import os
import sys
import numpy as np
import pandas as pd

listres_3letter = ["LEU","VAL","ALA","ILE","MET","PRO","GLY","PHE","TRP","TYR","SER","THR","ASN","GLN","HIS","LYS","ARG","ASP","GLU"]
listres = ['L','V','A','I','M','P','G','F','W','Y','S','T','N','Q','H','K','R','D','E']
listres1 = listres
listres2 = listres
assert(len(listres)==19)
N= len(listres)
reslist1 = ["LEU","VAL","ALA","ILE","MET","PRO","GLY","PHE","TRP","TYR","SER","THR","ASN","GLN","HIS","LYS","ARG","ASP","GLQ"]
reslistnames1 = listres
reslist2 = reslist1
reslistnames2 = reslistnames1

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

def getannotations(listres_1=listres,listres_2=listres):
  annotations=np.chararray((len(listres_1),len(listres_2)),itemsize=3)
  for i,xminusone in enumerate(listres_1):
     for j,xplusone in enumerate(listres_2):
        key = '%sT%s' %(xminusone,xplusone)
        annotations[i,j]=key
        print( i,j,key,annotations[i,j])
  return annotations


