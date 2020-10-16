import os
import sys
import pandas as pd
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import numpy as np
import seaborn as sns
import pickle
from HelperPlotting import *
from MultivariablePlottingGlobal import *
from datetime import date
import logging
import json
from DataSortingFiltering import *

sns.set(style="ticks")

#fields = None
#labels = None
#info = None
#labels_df = None
#runtype=None
#dpi = None
#nooutput = []

from MultivariablePlottingGlobal import *

def finishup():
  f=open('nooutput.txt','a')
  d = date.today()
  datestring = d.isoformat()
  #print datestring
  f.write('\n'+datestring+'\n')
  f.write('\n'.join(nooutput))
  f.close()

def get_full_pdbname(curfile,description):
    fullname = curfile.split('/')[0]+'/'+description+'.pdb.gz'
    #print(description,fullname)
    return fullname
    

def serialized_scoredata_name(inpfile,outpath):
    outname = get_outfile_name(inpfile)
    outfile = '%s/scorefile_dict_%s.p' %(outpath,outname)
    return outfile

def write_data_to_json(datadict,outfile,sortfield='interaction_energy'):
  sorteddict = sortdict_by_field_(datadict,sortfield)
  sorteddict['sorted']=True
  sorteddict['sortfield']=sortfield
  print(sorteddict['interaction_energy'][:20])
  json.dump(sorteddict,open(outfile,'w'))

def read_data_from_pickle(filename):
  print("filename")

def write_data_to_pickle(dictdata,outfile,sortfield='interaction_energy'):
  sorteddict = sortdict_by_field_(dictdata,sortfield)
  sorteddict['sorted']=True
  sorteddict['sortfield']=sortfield
  pickle.dump(sorteddict,open(outfile,'wb'))

def serialize_data(inpfile,use_fields=None,outpath = None,type_='serialize_data'):

    if outpath is None:
      outpath = 'data/scoredataframes/'
      os.system('mkdir -p %s' %outpath)

    outfile = serialized_scoredict_name(inpfile,outpath)

    if not os.path.exists(outfile):
      dictdata = getdatafromfile(inpfile)
      key = get_key(inpfile)
      dictdata['sequon']=[key for _ in range(0,Ndf)]
      dictdata['source']=[inpfile for _ in range(0,Ndf)]
      write_data_to_pickle(dictdata,outfile)


def clusters_to_dict(dictdata,n_clusters=None,k_means=False):
    from ClusterPoints import get_clusterpoints
    if n_clusters is not None:
      centroids , sizedict, sizedict_normalized,labels = get_clusterpoints(dictdata,n_clusters=n_clusters,k_means=k_means)
    else:
      centroids , sizedict, sizedict_normalized,labels = get_clusterpoints(dictdata,k_means)

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

def write_clusters_to_file(dictdata,curfile,n_clusters=7,k_means=False,use_topN =True,topN=200):
    df = pd.DataFrame()
    cdict = clusters_to_dict(dictdata,n_clusters=n_clusters,k_means=k_means)
    cdict['sequon']=get_key(curfile)
    cdict['isoenzyme']='T2'
    new_row = pd.DataFrame(cdict)
    df = pd.concat([new_row, df]).reset_index(drop = True)
    
    outfile = get_full_name(df,type_='serialize_clusters_df')
    pickle.dump(df,open(outfile,'wb'))
    del df


def get_axes_ticks():
    limits = {}
    limits_r = {} #reversed limits
    for axfield in info:
      if info[axfield] in fields:
        limits[axfield] = fields[ info[ axfield ] ][:2]

    ticks = {}
    ticks_r = {} #reversed ticks
    for axfield in limits:
      hi = limits[axfield][1]
      lo = limits[axfield][0]
      delta = 1.0
      if info[axfield] in delta_ticks:
        delta = delta_ticks[info[axfield]]
      #print("using ticks",delta_ticks)
      N = round((hi-lo) / delta) + 1
      ticks[axfield] = [ lo + delta*i for i in range(0,N)]

    for axfield in limits:
      ticks_r[axfield] = ticks[ axfield ][::-1] #reverse it
      limits[ axfield ][1] = ticks[ axfield ][-1] #fix max limits value
    for axfield in limits:
      limits_r[axfield] = limits[ axfield ][::-1]
    return limits,ticks,limits_r,ticks_r

def setup_axes(ax,raw=False,key=None,presets=True,ticksize=None,labelsize=18,view=None,view_abs=False,key_fontsize=20,**kwargs):

    if presets:
      axes_presets(ax)
    if not ticksize is None:
      ax.tick_params(labelsize=ticksize, grid_linewidth=1.0, grid_linestyle='-', grid_color='black',grid_alpha=0.6)

    ax.set_xlabel( labels[info['xfield']] ,fontsize=labelsize)
    ax.set_ylabel( labels[info['yfield']] ,fontsize=labelsize)
    ax.set_zlabel( labels[info['zfield']] ,fontsize=labelsize)

    limits , ticks, limits_r, ticks_r = get_axes_ticks()
    ax.set_yticks(ticks_r["yfield"])
    ax.set_xticks(ticks_r['xfield'])
    ax.set_zticks(ticks['zfield'])
    if not raw:
      ax.set(xlim = limits_r['xfield'] , ylim =  limits_r['yfield'] ,zlim = limits['zfield'])

    if not key is None:
      ax.text(fields[info['xfield']][0]+2.0,fields[info['yfield']][0]+2.0,fields[info['zfield']][1]-1.0,key,fontsize=key_fontsize)
    if view is None:
      ax.view_init(ax.elev  - 15 , ax.azim + 28)
    elif view_abs:
      ax.view_init(view['elev'], view['azim'])
    else:
      ax.view_init(ax.elev + view['elev'],ax.azim + view['azim'])
    

def plot3d(dictdata,key=None,raw=False,outfile=None, show=False,**kwargs):
    fig = plt.figure(figsize=(12,12))
    ax  = fig.add_subplot(111,projection='3d')
    im = ax.scatter( dictdata['xfield']['vals'] , dictdata['yfield']['vals'], dictdata['zfield']['vals'],c=dictdata['cfield']['vals'] , cmap=plt.cm.RdPu_r,vmin=-35, vmax=-5, alpha=0.8, edgecolors='w',lw=0,s=35)

    setup_axes(ax,raw,key)
    cbar = plt.colorbar(im ,ticks=[-5,-15,-25,-35],shrink = 0.4)
    cbar.ax.tick_params(labelsize=20)
    
    plt.tight_layout()
    if not outfile is None:
        plt.savefig(outfile,transparent=True,dpi=dpi[runtype])
    else:
        plt.show()
    if show:
        plt.show()
    plt.close()

def plot2dhist(collection,keys=None):
    limits , ticks, limits_r, ticks_r = get_axes_ticks()
    fig = plt.figure(figsize=(12,6))
    for i,dictdata in enumerate(collection):
      ax = fig.add_subplot(1,len(collection),i+1)
      im = ax.hist2d(dictdata['yfield']['vals'] , dictdata['zfield']['vals'],bins=250,cmap=plt.cm.GnBu,cmax=4)
      #ax = plt.gca()
      ax.set_yticks(ticks["zfield"])
      ax.set_xticks(ticks['yfield'])
      ax.set(xlim = limits['yfield'] , ylim =  limits['zfield'] )
      if not keys is None:
        ax.text(limits['yfield'][1]-2.0,limits['zfield'][1]-1.0,keys[i],fontsize=20)
      #plt.colorbar(im)
    #plt.show()
    plt.close()
    fig = plt.figure()
    for i,dictdata in enumerate(collection):
      ax = fig.add_subplot(1,len(collection),i+1)
      im = sns.kdeplot(dictdata['yfield']['vals'] , dictdata['zfield']['vals'],cmap=plt.cm.GnBu,ax=ax)
      ax.set_yticks(ticks["zfield"][::2])
      ax.set_xticks(ticks['yfield'][::2])
      ax.set(xlim = limits['yfield'] , ylim =  limits['zfield'] )
    plt.show()


def plot3d_compare(dictdata1, dictdata2):
    fig = plt.figure(figsize=(12,12))
    ax  = fig.add_subplot(111,projection='3d')
    markers = ['o','D']
    cmaps = [plt.cm.RdPu_r , plt.cm.GnBu_r]
    i=0
    for dictdata in [dictdata1 , dictdata2]:
      im = ax.scatter( dictdata['xfield']['vals'] , dictdata['yfield']['vals'], dictdata['zfield']['vals'],c=dictdata['cfield']['vals'] , cmap=cmaps[i],vmin=-35, vmax=-5, alpha=1.0, edgecolors='w',lw=0,s=60,marker=markers[i])
      i+=1
    setup_axes(ax)
    plt.colorbar(im)
    plt.show()

def k_means_distance(data,centroid,i_centroid,cluster_labels):
  if len(centroid)==3:
          [cx,cy,cz] = centroid
          distances = [np.sqrt((x-cx)**2+(y-cy)**2+(z-cz)**2) for (x, y, z) in data[cluster_labels == i_centroid]]
          return distances
  else:
          [cx,cy] = centroid
          distances = [np.sqrt((x-cx)**2+(y-cy)**2) for (x, y) in data[cluster_labels == i_centroid]]
          return distances



def add_cluster_centroids_to_axes_from_list(listnd,ax,annotate=True,n_clusters=None,nodiamonds=False,write_to_file=None,write_pdbs=True,**kwargs):
    from ClusterPoints import get_clusterpoints_from_list
    if n_clusters is not None:
      centroids , sizedict, sizedict_normalized, labels = get_clusterpoints_from_list(listnd,n_clusters=n_clusters)
    else:
      centroids , sizedict, sizedict_normalized,labels = get_clusterpoints_from_list(listnd)
   
    
    if not nodiamonds:
      if centroids.shape[1]==3:
        im2 = ax.scatter(centroids[:, 0], centroids[:, 1],centroids[:,2],marker='o', s=110, linewidths=3,
            color='black',alpha=1.0)
      else:
        im2 = ax.scatter(centroids[:, 0], centroids[:, 1],marker='o', s=110, linewidths=3,
            color='black',alpha=1.0)

    if annotate:
          i=0
          for k in sizedict_normalized:
            #TODO: for 3d as well
            #print(sizedict_normalized[k])
            label = "{:0.2f}".format(sizedict_normalized[k])
            #label2 = "{:.1f},\n{:.1f},\n{:.1f}".format(centroids[k, 0], centroids[k, 1],centroids[k,2])
            label3 = sizedict[k]
            if centroids.shape[1]==3:
              ax.text(centroids[k, 0]+0.2+0.01*i, centroids[k, 1]+0.2,centroids[k,2]+0.2,label,fontsize=11)
            else:
              #ax.text(centroids[k, 0]+0.2+0.01*i, centroids[k, 1]+0.2,label,fontsize=11)
              ax.text(centroids[k, 0]+0.2, centroids[k, 1]+0.3,label3,fontsize=14)
            i+=1

def add_cluster_centroids_to_axes(dictdata,ax,annotate=True,n_clusters=None,nodiamonds=False,write_to_file=None,k_means=False,fontsize_clusterlabel=14,**kwargs):
    from ClusterPoints import get_clusterpoints
    print(n_clusters)
    if n_clusters is not None:
      print(n_clusters)
      centroids , sizedict, sizedict_normalized,labels,n_clusters = get_clusterpoints(dictdata,n_clusters=n_clusters,k_means=k_means)
    else:
      print(n_clusters)
      centroids , sizedict, sizedict_normalized,labels,n_clusters = get_clusterpoints(dictdata,k_means=k_means)

    if not nodiamonds:
      im2 = ax.scatter(centroids[:, 0], centroids[:, 1],centroids[:,2],marker='o', s=110, linewidths=3,
            color='black',alpha=1.0)

    if annotate:
          i=0
          for k in sizedict_normalized:
            print("k ",k,sizedict_normalized[k])
            #print("k centroids ",k,centroids[k,:])
            if k<0: continue
            #print(sizedict_normalized[k])
            label = "{:0.2f}".format(sizedict_normalized[k])
            label2 = "{:.1f},\n{:.1f},\n{:.1f}".format(centroids[k, 0], centroids[k, 1],centroids[k,2])
            label3 = sizedict[k]
            #ax.text(centroids[k, 0]+0.2+0.01*i, centroids[k, 1]+0.2,centroids[k,2]+0.2,label,fontsize=11)
            #ax.text(centroids[k, 0]+0.23, centroids[k, 1]+0.23,centroids[k,2]+0.45,label2,fontsize=8)
            ax.text(centroids[k, 0]+0.23, centroids[k, 1]+0.23,centroids[k,2]+0.4,label3,fontsize=fontsize_clusterlabel)
            i+=1



def plot3d_compare_n(collection,keys=None,n=3,cluster=True,clusterOnly = False,annotate=True,outfile=None,**kwargs):

    fig = plt.figure(figsize=(22,10))
    markers = ['o', 'o','o']
    print("keys",keys)

    for i,dictdata in enumerate(collection):
      ax  = fig.add_subplot(1,n,i+1,projection='3d')
     
      if not clusterOnly: 
        alpha=0.8
        if cluster: alpha=0.2
        im1 = ax.scatter( dictdata['xfield']['vals'] , dictdata['yfield']['vals'], dictdata['zfield']['vals'],c=dictdata['cfield']['vals'] , cmap=plt.cm.RdPu_r,vmin=-35, vmax=-5, alpha=alpha, edgecolors='w',lw=0,s=20,marker=markers[i])
      if cluster:
        add_cluster_centroids_to_axes(dictdata,ax,n_clusters=7)
      if not keys is None:
        print("inside ",keys[i])
        setup_axes(ax,key=keys[i],ticksize=16,labelsize=18,view={'elev':-15, 'azim':55})
      else:
        setup_keys(ax)
      limits , ticks, limits_r, ticks_r = get_axes_ticks()
      ax.set_yticks(ticks_r["yfield"])
      ax.set_xticks(ticks_r['xfield'])
      ax.set_zticks(ticks['zfield'])
      ax.xaxis.set_rotate_label(True)
      ax.yaxis.set_rotate_label(True)
    plt.subplots_adjust(wspace=0.1)
    plt.subplots_adjust(hspace=0.1)
    plt.tight_layout()
    if not outfile is None:
      plt.savefig(outfile, transparent = True, dpi=dpi[runtype])
    plt.show()

def sortdict_by_field_(newdata,sortfield):
  if not sortfield in newdata:
    return None

  sorteddict ={}
  sortedindices = np.argsort(newdata[sortfield])
  for field in newdata:
    sorteddict[field] = []

  for ind in sortedindices:
    for field in newdata:
      sorteddict[field].append(newdata[field][ind])
  return sorteddict


def sortdict_by_field(org_dictdata,sortfield='cfield'):
  newdata = {'xfield':org_dictdata['xfield']['vals'],'yfield':org_dictdata['yfield']['vals'],'zfield':org_dictdata['zfield']['vals'],'cfield':org_dictdata['cfield']['vals'],'z2field':org_dictdata['z2field']['vals'],'description':org_dictdata['description']['vals']}
  return sortdict_by_field_(newdata,sortfield=sortfield)

def get_slice(sorteddict,threshold,field='cfield',presorted=True):
  return len(list(filter(lambda x: x <= threshold,sorteddict[field])))


def sort_cluster_by_field(sorteddict,pointindices_cluster,labels,distances=None,core_samples=None,sortorder = ['distance'],topN=20,**kwargs):
      indices_sorted = []
      for sfield in sortorder:

        if sfield=='distance':

          point_distances_cluster = [distances[val] for val in pointindices_cluster]
          
          indices_distance_sorted = np.argsort(point_distances_cluster)
          pointindices_cluster_sorted  = [pointindices_cluster[i] for i in indices_distance_sorted ]
          #indices_sorted = pointindices_cluster #WORKS ALWAYS
          
          indices_sorted = pointindices_cluster_sorted[:min(len(pointindices_cluster_sorted),topN)]
          for k in indices_sorted:
            print(k,distances[k],labels[k])

        if sfield=='interaction_energy':

          point_ie_cluster = [sorteddict['cfield'][i] for i in pointindices_cluster]
          indices_ie_sorted = np.argsort(point_ie_cluster)
          pointindices_cluster_sorted  = [pointindices_cluster[i] for i in indices_ie_sorted ]
          #indices_sorted = pointindices_cluster #WORKS ALWAYS
          indices_sorted = pointindices_cluster_sorted[:min(len(pointindices_cluster_sorted),topN)]
          for k in indices_sorted:
            print(k,sorteddict['cfield'][k],labels[k])

      tempdict={}
      for field in sorteddict:
        tempdict[field]=[sorteddict[field][k] for k in indices_sorted]
      return tempdict

def sort_clusters_by_field(sorteddict,centroids , sizedict, sizedict_normalized, labels, distances=None,core_samples=None,sortorder = ['interaction_energy'],topN=20):

  if distances is None:
    assert(core_samples is not None)
    assert(not 'distances' in sortorder)
  
  if core_samples is None:
    assert(distances is not None)

  def is_core_sample(isample):
    if core_samples is None:
      return True
    else:
      return (isample in core_samples)
  newclusterdict = {}    
  sorted_d = sorted(sizedict.items(), key=lambda x: x[1])
  print("SORTED dict",sorted_d)
  for cluster_id in sizedict:
      if cluster_id <0: continue #dbscan return negative clusters
      print("cluster ",cluster_id)
      pointindices_cluster=[]
      for i,v in enumerate(list(labels)):
        if v==cluster_id and is_core_sample(i):
          pointindices_cluster.append(i)
      print( pointindices_cluster)
      print("cluster size", len(pointindices_cluster),sizedict[cluster_id])
      if distances is None: 
        tempdict = sort_cluster_by_field(sorteddict,pointindices_cluster,labels,sortorder = sortorder,topN=topN)
      else:
        tempdict = sort_cluster_by_field(sorteddict,pointindices_cluster,labels,distances=distances[:,cluster_id],sortorder = sortorder,topN=topN)
      newclusterdict[cluster_id]={}
      newclusterdict[cluster_id]['cluster_size']=sizedict[cluster_id]
      newclusterdict[cluster_id]['dict']=tempdict

  return newclusterdict

def slice_and_cluster(sorteddict,thr,k_means=False,clus2d=False):
    
      topN = get_slice(sorteddict,thr)
      print(thr,topN)
      sorteddict_mod={}
      
      for ckey in sorteddict:
        print(ckey)
        sorteddict_mod[ckey]={}
        sorteddict_mod[ckey]['vals'] = sorteddict[ckey][:topN] #must do because of how i hav it setup :( NH
      if topN>15:
        from ClusterPoints import get_clusterpoints
        if k_means:
          #not tested for 2d option
          centroids , sizedict, sizedict_normalized, labels, distances= get_clusterpoints(sorteddict_mod,n_clusters=3,distances=True,k_means=k_means)
          return sort_clusters_by_field(sorteddict,centroids , sizedict, sizedict_normalized, labels, distances=distances)
        else:
          #listfromdf = np.array([ combined_df[labels_df[info['xfield']]] , combined_df[labels_df[info['yfield']]], combined_df[labels_df[info['zfield']]] ])
          #print(listfromdf.shape)
          centroids , sizedict, sizedict_normalized,labels,n_clusters,core_samples_ = get_clusterpoints(sorteddict_mod,k_means=k_means,core_samples=True)
          return sort_clusters_by_field(sorteddict_mod,centroids , sizedict, sizedict_normalized, labels, core_samples=core_samples_,sortorder=['interaction_energy'])
      else:
        return None
      
def write_sliced_cluster_pdbs(dictdata,thresholds=None,outfile=None,clus2d=False,k_means=False):

  if thresholds is None:
      thresholds=[-35]
  
  if outfile is None:
    if clus2d:
      outfile = 'results/clusters_2d/pdbfiles/test'
      os.system('mkdir -p results/clusters_2d/pdbfiles/')
    else:
      outfile = 'results/clusters/pdbfiles/test'
      os.system('mkdir -p results/clusters/pdbfiles/')
  thresholds.sort(reverse=True) #least negative to most negative 

  sorteddict = sortdict_by_field(dictdata)
  for thr in thresholds:
      newclusterdict = slice_and_cluster(sorteddict,thr,clus2d=clus2d,k_means=k_means)
      if newclusterdict is None:
        filename_p = outfile + '_thresh%03d' %thr + '_clustersize0000_pdbs.txt'
        outf=open(filename_p,'w')
        outf.write('\n')
        outf.close()
        return

      for cluster_id in newclusterdict:
        filename_p = outfile + '_thresh%03d' %thr + '_clustersize%04d_pdbs.txt'%newclusterdict[cluster_id]['cluster_size']
        #NEXT write out vals
        filename_v = outfile + '_thresh%03d' %thr + '_clustersize%04d_dict.p'%newclusterdict[cluster_id]['cluster_size']
        outf = open(filename_p,'w')
        for i,name in enumerate( newclusterdict[cluster_id]['dict']['description'] ):
            outstr=[]
            for axfield in ['description','xfield','yfield','zfield','cfield']:
              outstr.append( str(newclusterdict[cluster_id]['dict'][axfield][i]) )
            outstr.append( str(cluster_id) )
            outstr.append( str(newclusterdict[cluster_id]['cluster_size']) )
            outf.write('\t'.join(outstr)+'\n')
        outf.close()
        pickle.dump(newclusterdict,open(filename_v,'wb'))


def plot3d_slice_topN(collection,keys=None,outfile=None,show=False,cluster=False,color_field='cfield',cmap_field='cfield',dim3d=True,dim2d=False,topNs=[200],use_topN=True,**kwargs):


  figsize = (5*len(topNs),4*len(collection))
  fig = plt.figure(figsize=figsize)
  markers = ['o', 'o','o']
  print("keys",keys)
  if topNs is None:
    topNs = [200]#10%
  
  count=1
  cols = len(topNs)
  rows = len(collection)
  vmin = -40
  vmax = -25
  for i,org_dictdata in enumerate(collection):
    for topN in topNs:
      sorteddict = sortdict_by_field(org_dictdata)
      if not dim2d:
          ax  = fig.add_subplot(rows,cols,count,projection='3d')
          im1 = ax.scatter( sorteddict['xfield'][:topN] , sorteddict['yfield'][:topN], sorteddict['zfield'][:topN],c=sorteddict[cmap_field][:topN] , cmap=plt.cm.RdPu_r,vmin=vmin, vmax=vmax, alpha=1.0, edgecolors='w',lw=0,s=25,marker=markers[i])
          if cluster:
            sorteddict_mod={}
            for ckey in sorteddict:
              sorteddict_mod[ckey]={}
              sorteddict_mod[ckey]['vals'] = sorteddict[ckey][:topN]
            if topN>15:
              from ClusterPoints import get_clusterpoints
              centroid_clusters, counts_cluster,_,labels,dbscan_clusters = get_clusterpoints(sorteddict_mod,k_means=False)
              add_cluster_centroids_to_axes(sorteddict_mod,ax,annotate=True,n_clusters=3,nodiamonds=True)
          if not keys is None:
            setup_axes(ax,key=keys[i],ticksize=12,labelsize=10,view={'elev':17, 'azim':-37},view_abs=True,key_fontsize=15)
          else:
            setup_keys(ax,ticksize=15)
          label = '<{}'.format(topN)
          delta_ticks['substrate_ca_no_super_rmsd']=0.5
          limits , ticks, limits_r, ticks_r = get_axes_ticks()
          ax.set_yticks(ticks_r["yfield"])
          ax.set_xticks(ticks_r['xfield'])
          ax.set_zticks(ticks['zfield'][::2])
          if info['xfield']=='substrate_ca_no_super_rmsd':
            print('setting new limits')
            ax.set_xticks(ticks['xfield'])
            ax.set(xlim = [2.0,0.0] , ylim =  [6.5,2.5] ,zlim = limits['zfield'])
            ax.text(limits['xfield'][0]+0.5,6.5-0.3,limits['zfield'][1]-0.5,label,fontsize=14)
          ax.xaxis.set_rotate_label(True)
          ax.yaxis.set_rotate_label(True)
      count+=1
  plt.tight_layout()
  if not outfile is None:
      plt.savefig(outfile,transparent=True,dpi=dpi[runtype])
  if show:
    plt.show()
  plt.close()

def plot3d_slice_and_compare_n(collection,keys=None,n=3,thresholds=None,outfile=None,show=False,cluster=False,size=None,color_field='cfield',cmap_field='cfield',dim3d=True,dim2d=False,**kwargs):
    figsize = (12,4*len(collection))
    if not size is None:
      figsize=size
    fig = plt.figure(figsize=figsize)
    markers = ['o', 'o','o']
    print("keys",keys)
    
    if thresholds is None:
      thresholds=[-32,-30,-27]
   
    #sort threshold
    thresholds.sort(reverse=True) #least negative to most negative 
    count=1
    cols = len(thresholds)
    rows = len(collection)
    for i,org_dictdata in enumerate(collection):
      
      sorteddict = sortdict_by_field(org_dictdata)
      for j,thr in enumerate(thresholds):
      
        topN = get_slice(sorteddict,thr)
        print(thr,topN)
        vmax=-10
        vmin=-35
        if cmap_field=='z2field':
          vmax=3.0
          vmin=0.4
        if not dim2d:
          ax  = fig.add_subplot(rows,cols,count,projection='3d')
          im1 = ax.scatter( sorteddict['xfield'][:topN] , sorteddict['yfield'][:topN], sorteddict['zfield'][:topN],c=sorteddict[cmap_field][:topN] , cmap=plt.cm.RdPu_r,vmin=vmin, vmax=vmax, alpha=1.0, edgecolors='w',lw=0,s=25,marker=markers[i])
          if cluster:
            sorteddict_mod={}
            for ckey in sorteddict:
              sorteddict_mod[ckey]={}
              sorteddict_mod[ckey]['vals'] = sorteddict[ckey][:topN] #must do because of how i hav it setup :( NH
            if topN>15:
              from ClusterPoints import get_clusterpoints
              centroid_clusters, counts_cluster,_,labels,dbscan_clusters = get_clusterpoints(sorteddict_mod,k_means=False)
              print('DBSCAN', dbscan_clusters, counts_cluster,centroid_clusters)
              add_cluster_centroids_to_axes(sorteddict_mod,ax,annotate=True,n_clusters=3,nodiamonds=True) #max(dbscan_clusters-1,1)
              set_colors = plt.cm.rainbow(np.linspace(0,1,dbscan_clusters))
              clusters = list(set(labels))
              point_colors={}
              for ik,key in enumerate(clusters):               
                  point_colors[key] = set_colors[ik]
              colors = [point_colors[t] for t in labels]
              for ic,col in enumerate(colors):
                if labels[ic]==-1:
                  colors[ic]='black'
              print(len(sorteddict['xfield'][:topN]),len(colors))
              #im1 = ax.scatter( sorteddict['xfield'][:topN] , sorteddict['yfield'][:topN], sorteddict['zfield'][:topN],c=colors, alpha=1.0, edgecolors='w',lw=0,s=30,marker='o')
          if not keys is None:
            setup_axes(ax,key=keys[i],ticksize=12,labelsize=10,view={'elev':17, 'azim':-37},view_abs=True,key_fontsize=15)
          else:
            setup_keys(ax,ticksize=15)
          label = '<{}'.format(thr)
          delta_ticks['substrate_ca_no_super_rmsd']=0.5
          limits , ticks, limits_r, ticks_r = get_axes_ticks()
          ax.set_yticks(ticks_r["yfield"])
          ax.set_xticks(ticks_r['xfield'])
          ax.set_zticks(ticks['zfield'][::2])
          if info['xfield']=='substrate_ca_no_super_rmsd':
            print('setting new limits')
            ax.set_xticks(ticks['xfield'])
            ax.set(xlim = [2.0,0.0] , ylim =  [6.5,2.5] ,zlim = limits['zfield']) 
            ax.text(limits['xfield'][0]+0.5,6.5-0.3,limits['zfield'][1]-0.5,label,fontsize=14)
          ax.xaxis.set_rotate_label(True)
          ax.yaxis.set_rotate_label(True)
        if dim2d:
          ax  = fig.add_subplot(rows,cols,count)
          vmin=-40
          vmax=-20
          im1 = ax.scatter( sorteddict['yfield'][:topN], sorteddict['zfield'][:topN],c=sorteddict[cmap_field][:topN] , cmap=plt.cm.RdPu_r,vmin=vmin, vmax=vmax, alpha=1.0, edgecolors='w',lw=0,s=25,marker=markers[i])
          if cluster:
            sorteddict_mod={}
            for ckey in sorteddict:
              sorteddict_mod[ckey]={}
              sorteddict_mod[ckey] = sorteddict[ckey][:topN] #must do because of how i hav it setup :( NH
            if topN >15:
              add_cluster_centroids_to_axes_from_list([sorteddict_mod['yfield'],sorteddict['zfield']],ax,annotate=True,n_clusters=3,nodiamonds=True)
          label = '<{}'.format(thr)
          limits , ticks, limits_r, ticks_r = get_axes_ticks()
          ax.set_yticks(ticks["zfield"][::2])
          ax.set_xticks(ticks['yfield'][::2])
          ax.set(xlim =  limits['yfield'] ,ylim = limits['zfield'])
          ax.text(limits['yfield'][1]-1.5,limits['zfield'][1]-0.5,label,fontsize=20)
          ax.text(limits['yfield'][0]+0.1,limits['zfield'][1]-0.5,keys[i],fontsize=20)
        count+=1
        
    plt.tight_layout()
    if not outfile is None:
      plt.savefig(outfile,transparent=True,dpi=dpi[runtype])
      if show:
        plt.show()
    else:
      plt.show()

def write_sliced_cluster_pdbs_for_files(files,thresholds=None,k_means = False):
  basedir = 'results/clusters/pdbfiles/'
  os.system('mkdir -p %s' %basedir)

  if thresholds is None:
    thresholds = [-34,-35,-36,-37]

  for curfile in files:
    dictdata = getdata_xyzc(curfile)
    basedir = 'results/clusters/pdbfiles/'
    os.system('mkdir -p %s' %basedir)
    outfile = basedir +'/clusters_' + info['xfield']+'_'+ get_outfile_name(curfile)
    write_sliced_cluster_pdbs(dictdata,thresholds=thresholds,outfile=outfile,k_means=k_means)


def plot_seaborn_pairgrid_n(listofdicts,keys,outfile,show=True):
   df = get_combined_df(listofdicts,keys)
   sns.set_palette("husl")
   g = sns.PairGrid(df,vars=[labels_df[info['xfield']],labels_df[info['yfield']],labels_df[info['cfield']]],hue="sequon",palette="Set2")
   g.map_diag(sns.distplot, hist_kws={"histtype":'step','linewidth':2,"alpha":1},kde=False)
   g.map_offdiag(plt.scatter,s=12)
   g.add_legend()
   #for ax in g.axes.flat:
   # ax.plot((0, 50), (0, .2 * 50), c=".2", ls="--")
   plt.savefig(outfile,transparent=True,dpi=dpi[runtype])
   if show:
     plt.show()

def setup_pairgrid_axes(g):
   '''
    grid is drawn as follows
    
    x 00 01 02 
    y 10
    z 20    22
      x  y  z
   '''
   limits , ticks, limits_r, ticks_r = get_axes_ticks()
   for i in range(0,3):
    for j,curaxes in enumerate(['xfield','yfield','cfield']):
      g.axes[j,i].set_ylim(limits[curaxes][0],limits[curaxes][1])
      g.axes[i,j].set_xlim(limits[curaxes][0],limits[curaxes][1])
      g.axes[i,j].grid(b=True, which='major', color='grey', linewidth=0.3)
      if curaxes == 'cfield':
        g.axes[j,i].set_yticks(ticks[curaxes][::4])
        g.axes[i,j].set_xticks(ticks[curaxes][::4])
        g.axes[j,i].set_ylim(-42,-22)
        g.axes[i,j].set_xlim(-42,-22)
        #g.axes[j,i].xaxis.set_rotate_label(True)
      elif curaxes == 'xfield':
        g.axes[j,i].set_yticks([0.0,0.5,1.0,1.5,2.0,2.5])
        g.axes[i,j].set_xticks([0.0,0.5,1.0,1.5,2.0,2.5])
        g.axes[j,i].set_ylim(0.0,2.5)
        g.axes[i,j].set_xlim(0.0,2.5)
        #g.axes[i,j].grid(b=True, which='major', color='k', linewidth=0.6)
      else:
        g.axes[j,i].set_ylim(limits[curaxes][0],limits[curaxes][1])
        g.axes[i,j].set_xlim(limits[curaxes][0],limits[curaxes][1])
        g.axes[j,i].set_yticks(ticks[curaxes][::2])
        g.axes[i,j].set_xticks(ticks[curaxes][::2])

def plot_seaborn_pairgrid(dictdata,outfile=None,key=None,violin=False,show=True):
   mpl.rcParams.update(mpl.rcParamsDefault)
   sns.set(style="ticks")
   cleandict = {}
   listoffields = ['xfield','yfield','zfield','cfield','z2field']
   for curfield in listoffields:
     cleandict[labels_df[info[curfield]]] = dictdata[curfield]['vals']
   df = pd.DataFrame(cleandict)
   sns.set_palette("husl")
   g = sns.PairGrid(df,vars=[labels_df[info['xfield']],labels_df[info['yfield']],labels_df[info['cfield']]])
   if violin:
    g.map_diag(sns.violinplot)#, hist_kws={ "histtype":"step","linewidth":1.5,"alpha":1},kde=False,rug=False)
   else:
    g.map_diag(sns.distplot, hist_kws={ "histtype":"step","linewidth":1.5,"alpha":1},kde=False,rug=False)
   g.map_offdiag(plt.scatter,color='pink',s=7)

   
   setup_pairgrid_axes(g)
   if not key is None:
    g.axes[0,2].set_title(key)
   if not outfile is None:
    plt.savefig(outfile)
   if show:
     plt.show()

def plot_seaborn_distplot(dictdata,outfile,field,key=None,violin=False,show=True,raw=False):
   mpl.rcParams.update(mpl.rcParamsDefault)
   sns.set(style="ticks")
   cleandict = {}
   listoffields = ['xfield','yfield','zfield','cfield','z2field']
   for curfield in listoffields:
     if len(dictdata[curfield]['vals'])<1: continue
     cleandict[labels_df[info[curfield]]] = dictdata[curfield]['vals']
   df = pd.DataFrame(cleandict)
   sns.set_palette("husl")
   
   if violin:
    ax = sns.violinplot(df[labels_df[info[field]]],color='lightskyblue')#, hist_kws={ "histtype":"step","linewidth":1.5,"alpha":1},kde=False,rug=False)
   else:
    ax = sns.distplot(df[labels_df[info[field]]],color='lightskyblue')

   limits , ticks,limits_r, ticks_r = get_axes_ticks()

   ax.tick_params(labelsize=30, grid_linewidth=1.0, grid_linestyle='-', grid_color='black',grid_alpha=0.6)
   ax.grid(b=True, which='major', color='grey', linewidth=0.3)
   if not raw:
    ax.set_xlim(limits[field][0],limits[field][1])
    ax.set_xticks(ticks[field])

   if not key is None:
    ax.set_title(key)
   plt.savefig(outfile,transparent=True,dpi=600)
   if show:
     plt.show()

def plot_seaborn_jointplot(dictdata,outfile,field,key=None,violin=False,show=True,raw=False):
   mpl.rcParams.update(mpl.rcParamsDefault)
   sns.set(style="ticks")
   cleandict = {}
   listoffields = ['xfield','yfield','zfield','cfield','z2field']
   for curfield in listoffields:
     if len(dictdata[curfield]['vals'])<1: continue
     cleandict[labels_df[info[curfield]]] = dictdata[curfield]['vals']
   df = pd.DataFrame(cleandict)
   sns.set_palette("husl")

   if violin:
    ax = sns.violinplot(df[labels_df[info[field]]],color='pink')#, hist_kws={ "histtype":"step","linewidth":1.5,"alpha":1},kde=False,rug=False)
   else:
    ax = sns.distplot(df[labels_df[info[field]]],color='pink')

   limits , ticks,limits_r, ticks_r = get_axes_ticks()

   ax.tick_params(labelsize=30, grid_linewidth=1.0, grid_linestyle='-', grid_color='black',grid_alpha=0.6)
   ax.grid(b=True, which='major', color='grey', linewidth=0.3)
   if not raw:
    ax.set_xlim(limits[field][0],limits[field][1])
    ax.set_xticks(ticks[field])

   if not key is None:
    ax.set_title(key)
   plt.savefig(outfile)
   if show:
     plt.show()

def plot3d_from_combined_dataframe(files,outfile=None,key_all ='',show=True,view = [16,-26],category='all',expdata='all',filterbyfield=None):

  combined_df = get_filtered_combined_dataframe(files,category='all',expdata='all',filterbyfield=None)
  if combined_df is None: return
  if combined_df.shape[0]<5: return  

  print("Indices: ",combined_df.index[:])
  for col in combined_df.columns:
    print("col ",col)

  cbar_ticks = [-5,-15,-25,-35]
  vmax=-5
  vmin=-35
  print(outfile)

  if not filterbyfield is None:
    if 'interaction_energy' in filterbyfield:
      if filterbyfield['interaction_energy']=='high':
        vmin=-38
        vmax=-30
        cbar_ticks = [-32,-35,-38]

  fig = plt.figure(figsize=(12,12))
  ax  = fig.add_subplot(111,projection='3d')
  im = ax.scatter( combined_df[labels_df[info['xfield']]] , combined_df[labels_df[info['yfield']]], combined_df[labels_df[info['zfield']]],c=combined_df[labels_df[info['cfield']]] , cmap=plt.cm.RdPu_r,vmin=vmin, vmax=vmax, alpha=0.4, edgecolors='w',lw=0,s=15)

  setup_axes(ax,raw=False,key=key_all,view={'elev':view[0],'azim':view[1]},view_abs=True)
  ax.xaxis.set_rotate_label(True)
  ax.yaxis.set_rotate_label(True)
  if view[1]<-60:
    ax.zaxis.labelpad = 22
  if not filterbyfield is None:
    for key in filterbyfield:
      if key== 'interaction_energy' and filterbyfield[key]=='high':
       ax.set(xlim = [4.0,0.0] )
 
  cbar = plt.colorbar(im ,ticks=cbar_ticks,shrink = 0.4)
  cbar.ax.tick_params(labelsize=20)

  plt.tight_layout()
  if not outfile is None:
        plt.savefig(outfile,transparent=True,dpi=dpi[runtype])
  if show:
        plt.show()
  plt.close()



def plot_clusters_from_combined_dataframe(files,outfile=None,key_all ='',show=True,view = [16,-26],category='all',expdata='all',filterbyfield=None,k_means=False,nodiamonds=False,annotate=True,plot_core_samples=True):

  combined_df = get_filtered_combined_dataframe(files,category=category,expdata=expdata,filterbyfield=filterbyfield)
  if combined_df is None: return
  if combined_df.shape[0]<5: return
  print("Indices: ",combined_df.index[:])
  for col in combined_df.columns:
    print("col ",col)
  print("Shape ",combined_df.shape)

  listfromdf = np.array([ combined_df[labels_df[info['xfield']]] , combined_df[labels_df[info['yfield']]], combined_df[labels_df[info['zfield']]] ])
  print(listfromdf.shape)
  from ClusterPoints import get_clusterpoints_from_list
  centroids , sizedict, sizedict_normalized,labels,n_clusters,core_samples_ = get_clusterpoints_from_list(listfromdf,k_means=k_means,core_samples=True)  
  cbar_ticks = [-5,-15,-25,-35]
  vmax=-5
  vmin=-35
  print(outfile)

  if not filterbyfield is None:
    if 'interaction_energy' in filterbyfield:
      if filterbyfield['interaction_energy']=='high':
        vmin=-38
        vmax=-30
        cbar_ticks = [-32,-35,-38]

  fig = plt.figure(figsize=(12,12))
  ax  = fig.add_subplot(111,projection='3d')
  sizeMax = max(sizedict.values())
  fivePercent = int(sizeMax*0.01)
  print("Max, 1\% ",sizeMax,fivePercent)
  deleterows =[]
  for j in range(0,centroids.shape[0]):
    if sizedict[j] <=fivePercent:
      deleterows.append(j)
      del sizedict[j]

  clean_centroids = np.delete(centroids , deleterows,0)
  print(len(core_samples_))
  im2 = None
  if plot_core_samples:
    n_components = 3
    from sklearn.decomposition import IncrementalPCA
    ipca = IncrementalPCA(n_components=n_components, batch_size=10)
    X = np.array(list(map(list, zip(*listfromdf))))
    #X_ipca = ipca.fit_transform(X)
    del X
    #im2 = ax.scatter(X_ipca[:,0] , X_ipca[:,2], X_ipca[:,2],marker='o', s=10, linewidths=3,color='pink',alpha=0.5)

  im = ax.scatter(clean_centroids[:, 0], clean_centroids[:, 1],clean_centroids[:,2],marker='o', s=110, linewidths=3,
            color='black',alpha=1.0)

  if annotate:
          i=0
          sizeMax = max(sizedict.values())
          for k in sizedict:
            print("k ",k,sizedict_normalized[k])
            print("k centroids ",k,centroids[k,:])
            if k<0: continue
            label = "{:0.2f}".format(sizedict_normalized[k])
            label2 = "{:.1f}, {:.1f}, {:.1f}".format(centroids[k, 0], centroids[k, 1],centroids[k,2])
            label3 = sizedict[k]
            text_params = {'ha': 'center', 'va': 'center', 'family': 'sans-serif',
                   'fontweight': 'bold'}
            ax.text(centroids[k, 0]+0.1, centroids[k, 1]+0.1,centroids[k,2]+0.1,label3,fontsize=14,**text_params)
            ax.text(centroids[k, 0]+0.1, centroids[k, 1]+0.1,centroids[k,2]-0.1,label,fontsize=14,**text_params)
            ax.text(centroids[k, 0]+0.1, centroids[k, 1]+0.1,centroids[k,2]+0.25,label2,fontsize=11,color='red',**text_params)
            i+=1

  setup_axes(ax,raw=False,key=key_all,view={'elev':view[0],'azim':view[1]},view_abs=True)
  ax.xaxis.set_rotate_label(True)
  ax.yaxis.set_rotate_label(True)
  if view[1]<-60:
    ax.zaxis.labelpad = 22
  if not filterbyfield is None:
    for key in filterbyfield:
      if key== 'interaction_energy' and filterbyfield[key]=='high':
       ax.set(xlim = [4.0,0.0] )

  if not im2 is None:
    cbar = plt.colorbar(im ,ticks=cbar_ticks,shrink = 0.4)
    cbar.ax.tick_params(labelsize=20)

  plt.tight_layout()
  if not outfile is None:
    plt.savefig(outfile,transparent=True,dpi=dpi[runtype])
  if show:
    plt.show()
  plt.close()

def plot_clusters(dictdata,curfile,show=True,n_clusters=None,**kwargs):
      bData = True
      for axfield in info:
        if not 'vals' in dictdata[axfield]:
          bData=False
          break
        else:
          if len(dictdata[axfield]['vals'])<2:
            bData=False
            break
      if not bData:
        #print( dictdata['xfield'])
        print(dictdata['zfield']['vals'][:10])
        return bData

      print('here')
      fig = plt.figure(figsize=(8,8))
      ax = fig.add_subplot(111,projection='3d')

      add_cluster_centroids_to_axes(dictdata,ax,annotate=True,n_clusters=n_clusters,fontsize_clusterlabel=18)

      for axfield in info: delta_ticks[info[axfield]]=1.0
      setup_axes(ax,key=get_key(curfile),ticksize=16,labelsize=20,view={'elev':18,'azim':-45},view_abs=True)
      plt.tight_layout()

      outfile = get_full_name(curfile,'cluster')
      plt.savefig(outfile,transparent=True,dpi=dpi[runtype])
      if show:
          plt.show()
      plt.close()
      return bData

def plot_clusters_from_scorefiles(files,n_clusters=7,show=True,plot=True,write=False,**kwargs):
    os.system('mkdir -p results/clusters/pickled/')
    for curfile in files:
      fields = default_fields
      dictdata={}
      dictdata = getdata_xyzc(curfile)
      if dictdata is None:
        print("No output for file",curfile)
        nooutput.append(curfile)
        continue
      if plot:
        
        bData =  plot_clusters(dictdata,curfile,show=show,n_clusters=n_clusters)
        if not bData:
          print("No output for file",curfile)
          nooutput.append(curfile)
          continue 
      if write:
        write_clusters_to_file(dictdata,curfile,n_clusters=n_clusters)


 
def plot_pairgrid_from_scorefiles_n(files,add_xfield=True):
    mpl.rcParams.update(mpl.rcParamsDefault)
    basedir_pp = 'results/pairplots/'
    os.system('mkdir -p %s' %basedir_pp)
    sns.set(style="ticks")
    keys = [get_key(fcur) for fcur in files[:3]] #3 at a time
    if not add_xfield:
      outfile =  basedir_pp + 'pairgrid_'+ get_outfile_name_n(files)+'.png'
    else:
      outfile =  basedir_pp + 'pairgrid_'+info['xfield']+'_' + get_outfile_name_n(files)+'.png'
    listofdicts = [getdata_xyzc(fcur) for fcur in files]
    plot_seaborn_pairgrid_n(listofdicts,keys,outfile,show=False)

def plot_pairgrid_from_scorefiles(files,violin=False,add_xfield=False):
    mpl.rcParams.update(mpl.rcParamsDefault)
    sns.set(style="ticks")
    basedir_pp = 'results/pairplots/'
    os.system('mkdir -p %s' %basedir_pp)

    for curfile in files:
      dictdata = getdata_xyzc(curfile)
      if not add_xfield:
        if violin:
          outfile = basedir_pp + 'pairgridviolin_'+ get_outfile_name(curfile)+'.png'
        else:
          outfile = basedir_pp + 'pairgridviolin_'+ get_outfile_name(curfile)+'.png'
          
      else:
        outfile = basedir_pp + 'pairgrid_'+ info['xfield']+'_' + get_outfile_name(curfile)+'.png'
      key=get_key(curfile)
      
      plot_seaborn_pairgrid(dictdata,outfile,key=key,violin=violin)
      #except:
      #  print('could not generate plot for ',curfile)

def plot_distribution_from_scorefiles(files,field,violin=False,show=True,tag='',raw=False,threshold=None):
    mpl.rcParams.update(mpl.rcParamsDefault)
    sns.set(style="ticks")
    basedir_pp = 'results%s/distribution/' %tag
    os.system('mkdir -p %s' %basedir_pp)

    for curfile in files:
      dictdata = getdata_xyzc(curfile)
      if not threshold is None:
        dictdata_ = sortdict_by_field(dictdata)
        dictdata = get_slice(dictdata_,threshold)
      if violin:
          outfile = basedir_pp + '_' + info[field] + '_violin_'+ get_outfile_name(curfile)+'.png'
      else:
          outfile = basedir_pp + '_' + info[field] + '_dist_'+ get_outfile_name(curfile)+'.png'

      key=get_key(curfile)

      plot_seaborn_distplot(dictdata,outfile,field,key=key,violin=violin,show=show,raw=raw)

def plot3d_from_scorefiles(files,show=False):

    basedir_3d = 'results/scatter3d/'
    os.system('mkdir -p %s' %basedir_3d)

    for curfile in files:
      dictdata = getdata_xyzc(curfile)
      outfile = get_full_name(curfile,type_='3d')
      #try:
      plot3d(dictdata,key=get_key(curfile),outfile = outfile,show=show)
      #except:
      #  print('could not generate plot for ',curfile)

def plot_graphs_from_scorefiles(files):

    basedir_3d = 'results/scatter3d/'
    os.system('mkdir -p %s' %basedir_3d)
    basedir_pp = 'results/pairplots/'
    os.system('mkdir -p %s' %basedir_pp)

    for curfile in files:

      dictdata = getdata_xyzc(curfile)
      print( curfile)
      outfile = basedir_3d + 'scatter3d_'+ get_outfile_name(curfile)+'.png'
      if not os.path.exists(outfile):
        plot3d(dictdata,key=get_key(curfile),outfile = outfile)

      outfile = basedir_pp + 'pairgrid_'+ get_outfile_name(curfile)+'.png'
      if not os.path.exists(outfile):
        plot_seaborn_pairgrid(dictdata,outfile)

