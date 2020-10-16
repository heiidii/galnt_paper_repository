import os
import sys
import pandas as pd
import numpy as np
import logging

def really_safe_normalise_in_place(d):
      import math
      import operator
      factor=1.0/math.fsum(d.values())
      d_n = {}
      for k in d:
        d_n[k] = d[k]*factor
      key_for_max = max(d_n.items(), key=operator.itemgetter(1))[0]
      diff = 1.0 - math.fsum(d_n.values())
      d_n[key_for_max] += diff
      return d_n

def get_clusterpoints(dictdata,n_clusters=7,k_means=True,cluster_fields=['xfield','yfield','zfield'],distances=False,core_samples=False):

    l = [dictdata['xfield']['vals'],dictdata['yfield']['vals'],dictdata['zfield']['vals']]
    #l = [ [ dictdata[field]['vals'] for field in cluster_fields ] ]#generalized for nay number of dims
    
    return get_clusterpoints_from_list(l,n_clusters=n_clusters,k_means=k_means,distances=distances,core_samples=core_samples)

def get_clusterpoints_from_list(l,n_clusters=7,k_means=True,distances=False,core_samples=False):
    if k_means:
      return  k_means_cluster(l,n_clusters,distances=distances)
    else:
      return dbscan_cluster(l,core_samples=core_samples)

def k_means_cluster(l,n_clusters,distances=False,**kwargs):
    from sklearn import cluster
    #print('k_means; using clusters ',n_clusters)
    k_means = cluster.KMeans(init='k-means++', n_init=10,n_clusters=n_clusters)
    X = np.array(list(map(list, zip(*l))))
    k_means.fit(X)
    centroids = k_means.cluster_centers_
    from collections import Counter
    counts_cluster  = Counter(k_means.labels_)
    counts_cluster_normalized = really_safe_normalise_in_place(counts_cluster)
    if distances:
      X_dist = k_means.transform(X)
      return centroids,counts_cluster,counts_cluster_normalized,list(k_means.labels_),X_dist
    else:
      return centroids,counts_cluster,counts_cluster_normalized,list(k_means.labels_)

def get_centroids(labels,X):
    unique_labels = set(labels)
    for elem in unique_labels:
      if elem<0:
        unique_labels.remove(elem)
        break
    centroids = np.zeros((len(unique_labels),X.shape[1])) #removing -ve 
    for ic,clus in enumerate(unique_labels):
      points_indices = []
      if clus>=0: #Non-negative - zero included. Neg clusters are non-clusters
        for ip,il in enumerate(labels):
          if il==clus:
            points_indices.append(ip)
        if len(points_indices)>0:
          Y = np.array([X[ip,:] for ip in points_indices])
          Ymean = np.mean(Y,axis=0)
          centroids[clus,:]=Ymean
    
    return centroids
    

def dbscan_cluster(l,core_samples=True,**kwargs):
    from sklearn.cluster import DBSCAN
    X = np.array(list(map(list, zip(*l))))
    db_ = DBSCAN(eps=0.3, min_samples=10)
    db = db_.fit(X) 
    from collections import Counter
    counts_cluster  = Counter(db.labels_)
    labels = list(db.labels_)
    n_clusters = len(set(labels))
    counts_cluster_normalized = really_safe_normalise_in_place(counts_cluster) 
    centroids = get_centroids(labels,X) #for labeling plots

    if core_samples:
      return centroids, counts_cluster,counts_cluster_normalized,labels,n_clusters,list(db.core_sample_indices_)
    else:
      return centroids, counts_cluster,counts_cluster_normalized,labels,n_clusters
    
