from enum import Enum, unique , auto

hb_cutoff = 4.0
rmsd_cutoff = 1.0
sequonrmsd_cutoff = 0.9

class hb_bond_criteria(Enum):
  significant_clusters_cutoff = auto()
  significant_clusters = auto()
  significant_clusters_largest_cluster = auto()
  largest_cluster = auto()

  significant_clusters_cutoff_rmsd = auto()
  significant_clusters_cutoff_rmsd_only = auto()
  significant_clusters_cutoff_sequonrmsd = auto()
  
  significant_clusters_cutoff_sequonrmsd_only = auto()

  significant_clusters_sinks = auto()
  significant_clusters_sinks_TopKCal = auto()
  significant_clusters_sinks_cutoff = auto()
  significant_clusters_sinks_TopKCal_cutoff = auto()
  significant_clusters_sinks_cutoff_rmsd = auto()
  significant_clusters_sinks_cutoff_rmsd_only = auto()
  significant_clusters_sinks_cutoff_dmc_only = auto()

  significant_clusters_sinks_TopKCal_cutoff_rmsd = auto()
  significant_clusters_sinks_TopKCal_cutoff_rmsd_only = auto()

  significant_clusters_sinks_GTX_cutoff_rmsd = auto()
  significant_clusters_sinks_GTX_TopKCal_cutoff_rmsd = auto()

  significant_clusters_sinks_cutoff_sequonrmsd = auto()
  significant_clusters_sinks_cutoff_sequonrmsd_only = auto()
  significant_clusters_sinks_TopKCal_cutoff_sequonrmsd = auto()
  significant_clusters_sinks_TopKCal_cutoff_sequonrmsd_only = auto()

  significant_clusters_sinks_GTX_cutoff_sequonrmsd = auto()
  significant_clusters_sinks_GTX_TopKCal_cutoff_sequonrmsd = auto()

  significant_clusters_sinks_TTY_cutoff = auto()
  significant_clusters_sinks_TTY_TopKCal_cutoff = auto()

  significant_clusters_sinks_cutoff_rmsd_ie = auto()

def get_type_value(criteria):
  type_value = 'centroids'
  if criteria==hb_bond_criteria.significant_clusters_sinks or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only or criteria == hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_rmsd  or criteria == hb_bond_criteria.significant_clusters_sinks_TTY_cutoff or criteria == hb_bond_criteria.significant_clusters_sinks_cutoff_sequonrmsd or criteria == hb_bond_criteria.significant_clusters_sinks_cutoff_sequonrmsd_only or criteria == hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_sequonrmsd or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_ie:
    type_value='sinks'
  if criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal or criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff or criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff_rmsd or criteria == hb_bond_criteria.significant_clusters_sinks_TTY_TopKCal_cutoff  or criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff_sequonrmsd or criteria == hb_bond_criteria.significant_clusters_sinks_GTX_TopKCal_cutoff_sequonrmsd:
    type_value='sinks_TopKCal'

  return type_value


def check_criteria_decoy(values,criteria=hb_bond_criteria.significant_clusters_sinks_cutoff):

  [value1,value2] = values

  print('values',value1)
  if criteria==hb_bond_criteria.significant_clusters_cutoff or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff or  criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff:
      return value1 < hb_cutoff

  if criteria==hb_bond_criteria.significant_clusters_sinks_TTY_cutoff or  criteria==hb_bond_criteria.significant_clusters_sinks_TTY_TopKCal_cutoff:
      return (value1 > hb_cutoff and value1 < 5.5 and value2 < 1.6)

  if criteria==hb_bond_criteria.significant_clusters_cutoff_rmsd_only or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only or  criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff_rmsd_only:
      return value2 < rmsd_cutoff

  if criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_sequonrmsd_only or criteria==hb_bond_criteria.significant_clusters_cutoff_sequonrmsd_only:
      return value2 < sequonrmsd_cutoff

  if criteria==hb_bond_criteria.significant_clusters_cutoff_rmsd or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd or  criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff_rmsd:
      return ( value1 < hb_cutoff and value2 < rmsd_cutoff )

  if criteria==hb_bond_criteria.significant_clusters_cutoff_sequonrmsd or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_sequonrmsd or  criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff_sequonrmsd:
      return ( value1 < hb_cutoff and value2 < sequonrmsd_cutoff )

  if criteria==hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_rmsd or  criteria==hb_bond_criteria.significant_clusters_sinks_GTX_TopKCal_cutoff_rmsd:
      return ( value1 <= hb_cutoff and value2 >= rmsd_cutoff )

def check_hb_criteria(centroids,sizedict,sizedict_normalized,criteria=hb_bond_criteria.significant_clusters,axis=1):
  fieldname = 'distance_catalysis_HWUO1B-THR7N'

  cluster_ids=[]
  axis_rmsd = 0
  axis_sequonrmsd = 3 # 0-x,1-y,2-z,3-z3
  if criteria==hb_bond_criteria.significant_clusters_cutoff or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff or  criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff:
    for k in sizedict:
      if sizedict_normalized[k]<0.05: continue
      if k == -1: continue
      if centroids[k,axis]<hb_cutoff:
        cluster_ids.append(k)


  if criteria==hb_bond_criteria.significant_clusters_sinks_TTY_cutoff or  criteria==hb_bond_criteria.significant_clusters_sinks_TTY_TopKCal_cutoff:
    for k in sizedict:
      if sizedict_normalized[k]<0.05: continue
      if k == -1: continue
      if (centroids[k,axis]>hb_cutoff and centroids[k,axis]<5.5 and centroids[k,axis_rmsd] < 1.6):
        cluster_ids.append(k)

  if criteria==hb_bond_criteria.significant_clusters_cutoff_rmsd_only or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd_only or  criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff_rmsd_only:
    for k in sizedict:
      if sizedict_normalized[k]<0.05: continue
      if k == -1: continue
      if centroids[k,axis_rmsd]<rmsd_cutoff:
        cluster_ids.append(k)

  if criteria==hb_bond_criteria.significant_clusters_cutoff_sequonrmsd_only or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_sequonrmsd_only or  criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff_sequonrmsd_only:
    for k in sizedict:
      if sizedict_normalized[k]<0.05: continue
      if k == -1: continue
      if centroids[k,axis_sequonrmsd]<sequonrmsd_cutoff and centroids[k,axis_rmsd]<rmsd_cutoff:
        cluster_ids.append(k)

  if criteria==hb_bond_criteria.significant_clusters_cutoff_rmsd or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_rmsd or  criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff_rmsd:
    for k in sizedict:
      if sizedict_normalized[k]<0.05: continue
      if k == -1: continue
      #print('centroids ', centroids[k,:])
      if centroids[k,axis]<=hb_cutoff and centroids[k,axis_rmsd]<rmsd_cutoff:
        cluster_ids.append(k)
    #print("CRITERIA CUTOFF ",cluster_ids)

  if criteria==hb_bond_criteria.significant_clusters_cutoff_sequonrmsd or criteria==hb_bond_criteria.significant_clusters_sinks_cutoff_sequonrmsd or  criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal_cutoff_sequonrmsd:
    for k in sizedict:
      if sizedict_normalized[k]<0.05: continue
      if k == -1: continue
      #print('centroids ', centroids[k,:])
      if centroids[k,axis]<=hb_cutoff and centroids[k,axis_sequonrmsd]<sequonrmsd_cutoff:
        cluster_ids.append(k)
    #print("CRITERIA CUTOFF ",cluster_ids)

  if criteria==criteria==hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_rmsd or  criteria==hb_bond_criteria.significant_clusters_sinks_GTX_TopKCal_cutoff_rmsd:
    for k in sizedict:
      if sizedict_normalized[k]<0.05: continue
      if k == -1: continue
      #print('centroids ', centroids[k,:])
      if centroids[k,axis]<=hb_cutoff and centroids[k,axis_rmsd] >= rmsd_cutoff:
        cluster_ids.append(k)
    #print("CRITERIA CUTOFF ",cluster_ids)

  if criteria==criteria==hb_bond_criteria.significant_clusters_sinks_GTX_cutoff_sequonrmsd or  criteria==hb_bond_criteria.significant_clusters_sinks_GTX_TopKCal_cutoff_sequonrmsd:
    for k in sizedict:
      if sizedict_normalized[k]<0.05: continue
      if k == -1: continue
      #print('centroids ', centroids[k,:])
      if centroids[k,axis]<=hb_cutoff and centroids[k,axis_sequonrmsd] >= sequonrmsd_cutoff:
        cluster_ids.append(k)
    #print("CRITERIA CUTOFF ",cluster_ids)

  if criteria==hb_bond_criteria.significant_clusters or criteria==hb_bond_criteria.significant_clusters_sinks or criteria==hb_bond_criteria.significant_clusters_sinks_TopKCal:
    for k in sizedict:
      #if sizedict_normalized[k]<0.05: continue
      if k == -1: continue
      cluster_ids.append(k)

  if criteria==hb_bond_criteria.significant_clusters_largest_cluster:
      maxclus = max(sizedict_normalized, key = lambda x: sizedict_normalized.get(x) )
      cluster_ids.append(maxclus)

  if criteria==hb_bond_criteria.largest_cluster:
      maxclus = max(sizedict_normalized, key = lambda x: sizedict_normalized.get(x) )
      cluster_ids.append(maxclus)

  return cluster_ids
