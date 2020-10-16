import os
import sys
from pyrosetta import *
from rosetta import *
pyrosetta.init('-include_sugars','-mute core -mute basic.io.database')
from rosetta.core.scoring import *
from rosetta.protocols.docking import *

import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
#Setting up mover for viewing structures

import argparse
from rosetta.core.pose import *
import math
from rosetta.core.simple_metrics import metrics
sf = ScoreFunction()
sf = get_fa_scorefxn()

def get_score_string_from_score_type(st):
  strst = str(st)
  return strst.split('.')[1]

def PairwiseEnergyMapForPdbFile(curfile,tag='notag',residues=[500],filter_pep_pep_energies=True):
  print(tag)
  curpose = pose_from_pdb(curfile)
  sf_pose = sf(curpose)
  dict_ = {}
  columns=['residue_i','residue_j','tag','file','score_type','score_value']
  for key in columns:
    dict_[key] = []
  for resi in residues:
    score_types = sf.get_nonzero_weighted_scoretypes()
    mylist = ['rama','fa_rep','fa_atr','hbond_bb_sc','hbond_sc','fa_sol','fa_intra_rep']
    #mylist = ['fa_elec']
    for st in score_types:
      strst_clean = get_score_string_from_score_type(st)
      if strst_clean in mylist:
        score_tuple = pyrosetta.toolbox.atom_pair_energy._reisude_pair_energies(resi,curpose,sf,st,0.25)
        #print(strst_clean,score_tuple[0],score_tuple[1])
        for entry in score_tuple:
          if filter_pep_pep_energies:
            if entry[0] in range(497,505):
              continue
          #print(strst_clean,entry[0],entry[1])
          dict_['score_type'].append(strst_clean)
          dict_['score_value'].append(entry[1])
          dict_['residue_j'].append(entry[0])
          dict_['residue_i'].append(resi)
          dict_['file'].append(curfile)
          dict_['tag'].append(tag)
  return dict_


def EnergyMapForPdbFile(curfile,tag='notag',residues=[498,499,500],columns=['residue','tag','file','score_type','score_value']):

  columns_total = ['tag','file','score_type','score_value']
  dict_ ={}
  dict_total = {}
  for key in columns_total:
    dict_total[key]=[]
  for key in columns:
    dict_[key]=[]
  curpose = pose_from_pdb(curfile)
  sf_pose = sf(curpose)
  pose_energies = curpose.energies()
  total_e = pose_energies.total_energies()
  total_e = pose_energies.active_total_energies()
  score_types = sf.get_nonzero_weighted_scoretypes()
  mylist = ['rama','fa_rep','fa_atr','hbond_bb_sc','hbond_sc','fa_sol','fa_intra_rep']
  for st in score_types:
      for resi in residues:
        strst = str(st)
        strst_clean = strst.split('.')[1]
        if strst_clean in mylist:
          value = pose_energies.residue_total_energies(resi)[st] 
          
          dict_['score_value'].append(value)
          dict_['score_type'].append(strst_clean)
          dict_['residue'].append(resi)
          dict_['file'].append(curfile)
          dict_['tag'].append(tag)
  return dict_



def EnergyMapForPdbFiles(files,tags=None,plotresenergies=True,plottotenergies=False,suffix=None,residues=[498,499,500],columns=['residue','tag','file','score_type','score_value'],dump_single_df=False,plot_single_df=False):
  dfs = []
  for ifile,curfile in enumerate(files):
    dict_ = EnergyMapForPdbFile(curfile,tag=tags[ifile],plotresenergies=plotresenergies,plottotenergies=plottotenergies,residues=residues,columns=columns)
    tidy_df = pd.DataFrame()
    if dump_single_df:
      print(tidy_df)
      pickle.dump(tidy_df,open('PickleFiles/Dataframetidy_residues.p','wb'))
    if plot_single_df:
      sns.barplot(x='score_type',y='score_value',data=tidy_df,hue='tag')
      plt.show()
      plt.close()
    for column in columns:
      tidy_df[column] = dict_[column]
    dfs.append(tidy_df)
  df_combined = pd.concat(dfs,ignore_index=True)

  print(df_combined)
  pickle.dump(df_combined,open('PickleFiles/Dataframetidy_residues.p','wb'))

  for residue in residues:
    if suffix is None:
      outfilename = 'results/residue_energies/residuewiseenergies_%s.png' %(residue)
    else:
      outfilename = 'results/residue_energies/residuewiseenergies_%s_%s.png' %(residue,suffix)
    fig = plt.figure()
    df_combined_filter  = df_combined[ df_combined['residue'] == residue ]
    ax = sns.barplot(x='score_type',y='score_value',data=df_combined_filter,hue='tag')
    xlabelsize=12
    ax.tick_params(axis='x',labelsize=xlabelsize, grid_linewidth=0.5, grid_linestyle='-', grid_color='black',grid_alpha=0.6,labelrotation=45)
    ax.tick_params(axis='y', labelcolor='gray',labelsize=15, grid_linewidth=0.5, grid_linestyle='-', grid_color='black',grid_alpha=0.6)
    ax.text(4.5,1.5,'res %s'%residue)
    ax.set(ylim=[-6.5,6.5])
    plt.tight_layout()
    plt.savefig(outfilename,transparent=True,dpi=1200)
    del df_combined_filter
    #plt.show()
    plt.close()



