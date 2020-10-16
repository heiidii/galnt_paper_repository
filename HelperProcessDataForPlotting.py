import os
import sys
import pandas as pd
import pickle

def get_df_from_dict(dictdata,labels_df,entries=None):
    cleandict = {}
    if entries==None:
      entries = ['xfield','yfield','zfield','cfield','z2field']
    for curfield in entries:
      cleandict[labels_df[info[curfield]]] = dictdata[curfield]['vals']
    df = pd.DataFrame(cleandict)
    return df

def serialized_scoredata_name(inpfile,outpath):
    outname = get_outfile_name(inpfile)
    outfile = '%s/scorefile_dict_%s.p' %(outpath,outname)
    return outfile

def getdatafromfile(curfile,inuse_fields=None):
    f = open(curfile,'r')
    lines = f.readlines()
    f.close()

    print( fields)
    print( info)
    data = dict()
    from functions_filterpdbforparams import getfiltereddatafromlines
    if use_fields is None:
      fields_local = fields
      data, fields_temp = getfiltereddatafromlines(lines,fields_local)
    if use_fields=='all':
      data, fields_temp = getfiltereddatafromlines(lines)
    return data

def serialize_data(inpfile,use_fields=None,outpath = None):

    if outpath is None:
      outpath = 'data/scoredataframes/'
      os.system('mkdir -p %s' %outpath)

    outfile = serialized_scoredict_name(inpfile,outpath)

    if not os.path.exists(outfile):
      dictdata = getdatafromfile(inpfile)
      key = get_key(inpfile)
      dictdata['sequon']=[key for _ in range(0,Ndf)]
      dictdata['source']=[inpfile for _ in range(0,Ndf)]

      f = open(outfile,'wb')
      pickle.dump(df,f)
      f1.close()

def clusters_to_dict(dictdata,n_clusters=None):
    from ClusterPoints import get_clusterpoints
    if n_clusters is not None:
      centroids , sizedict, sizedict_normalized = get_clusterpoints(dictdata,n_clusters=n_clusters)
    else:
      centroids , sizedict, sizedict_normalized = get_clusterpoints(dictdata)

    cdict={}
    for k in sizedict_normalized:
      cdict['centroid_id']=k
      cdict['dimensions']=centroids.shape[1]
      cdict['n_clusters']=centroids.shape[0]
      cdict['center']=centroids[k,:]
      cdict['size']=sizedict[k]
      cdict['size_norm']=sizedict_normalized[k]
      cdict['xfield']=info['xfield']
      cdict['yfield']=info['yfield']
      cdict['zfield']=info['zfield']
    return cdict

def write_clusters_to_file(dictdata,n_clusters=7):
        df = pd.DataFrame()
        cdict = clusters_to_dict(dictdata,n_clusters=n_clusters)
        cdict['sequon']=get_key(curfile)
        cdict['isoenzyme']='T2'
        new_row = pd.DataFrame(cdict)
        df = pd.concat([new_row, df]).reset_index(drop = True)
        df_all = pd.concat([new_row, df_all]).reset_index(drop = True)

        pickle.dump(df,open('results/clusters/pickled/clusters_'+cdict['sequon']+'_'+'_'+cdict['isoenzyme']+info['xfield']+'.p','wb'))
        del df


