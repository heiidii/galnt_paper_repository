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

def EnergyMapForPdbFilesOld(files,tags=None,plotresenergies=True,plottotenergies=False):

  energy_terms = [fa_atr,fa_rep]
  residues=range(498,500+1)
  dict_array={'tag':[], 'file':[]}
  pmm = PyMOLMover()
  dict_array_total = {'tag':[], 'file':[]}
  tidy_df = pd.DataFrame(columns=['residue','tag','file','score_type','score_value'])
  
  for ifile,curfile in enumerate(files):
    curpose = pose_from_pdb(curfile)
    sf_pose = sf(curpose)
    pmm.apply(curpose)
    #sf.show(curpose)
    pose_energies = curpose.energies()
    total_e = pose_energies.total_energies()
    total_e = pose_energies.active_total_energies()
    #sf_info = core.scoring.ScoreFunctionInfo()
    #sf_info.initialize_from(sf)
    #print(sf_info.scores_present())
    #for resi in residues:
    #  pyrosetta.toolbox.atom_pair_energy.print_residue_pair_energies(resi,curpose,sf)
    
    score_types = sf.get_nonzero_weighted_scoretypes()
    dict_total = {}
    for st in score_types:
      #print(st)
      #print(str(st))
      dict_total[str(st)] = pose_energies.total_energies()[st]
      #for resi in [498]:
      #  print('\n%s %s\n'%(resi,str(st)))
      #  pyrosetta.toolbox.atom_pair_energy.print_residue_pair_energies(resi,curpose,sf,st,0.25)
    #exit()
    for key in dict_total:
        if not key in dict_array_total:
          dict_array_total[key]=[]
        dict_array_total[key].append(dict_total[key])
    dict_array_total['tag'].append(tags[ifile])
    dict_array_total['file'].append(files[ifile])
    
    dict_={}
    for resi in residues:
      
      erama = pose_energies.residue_total_energies(resi)[rama]
      e_fa_rep = pose_energies.residue_total_energies(resi)[fa_rep]
      e_fa_atr = pose_energies.residue_total_energies(resi)[fa_atr]
      e_hb_bb_sc = pose_energies.residue_total_energies(resi)[hbond_bb_sc]
      e_hb_sc = pose_energies.residue_total_energies(resi)[hbond_sc]
      e_sol = pose_energies.residue_total_energies(resi)[fa_sol]
      e_intra_rep = pose_energies.residue_total_energies(resi)[fa_intra_rep]
      dict_['fa_sol']=e_sol
      dict_['hb_bb_sc']=e_hb_bb_sc
      dict_['hb_sc']=e_hb_sc
      dict_['fa_intra_rep'] = e_intra_rep
      
      #res_ene_array_2b = pose_energies.residue_pair_energies_array([resi])
      #sf.eval_cd_2b_sc_sc(curpose.residue(94),curpose.residue(98),curpose,emap)
      #e1 = pose_energies.residue_total_energies(resi)
      dict_['total_fa_rep']=e_fa_rep
      dict_['total_fa_atr']=e_fa_atr
      dict_['rama']=erama
      dict_['residue']=resi
      for key in dict_:
        if not key in dict_array:
          dict_array[key]=[]
        dict_array[key].append(dict_[key])
      dict_array['tag'].append(tags[ifile])
      dict_array['file'].append(files[ifile])
  df = pd.DataFrame.from_dict(dict_array)
  pickle.dump(df,open('Dataframe_H3files_residues.p','wb'))
  df_total = pd.DataFrame.from_dict(dict_array_total)
  #df_temp = df_total[df_total['tag']=='native']
  #df_temp.reset_index(drop=True)
  #print(df_temp)

  pickle.dump(df_total,open('Dataframe_H3files_total.p','wb'))
  import seaborn as sns
  import matplotlib.pyplot as plt
  if plotresenergies:
    print('plotting res energies')
    for thirdrow in ['fa_sol','hb_bb_sc','hb_sc','fa_intra_rep']:
      fig, ax = plt.subplots(nrows=3,ncols=1,figsize=(10,18))
      sns.barplot(x='residue',y='total_fa_rep',hue='tag',data=df,palette='Set2',ax=ax[0])
      sns.barplot(x='residue',y='total_fa_atr',hue='tag',data=df,palette='Set2',ax=ax[1])
      sns.barplot(x='residue',y=thirdrow,hue='tag',data=df,palette='Set2',ax=ax[2])
      #plt.savefig('residuewiseenergies_%s.png' %thirdrow,transparent=True,dpi=1200)
      #plt.show()
      plt.close()
    

  '''
  print(df_total)


  if plottotenergies:
    score_types_filter=['fa_rep','fa_atr']
    for col in score_types_filter:
      
      fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(12,8))
      sns.pointplot(x='tag',y=str(col),data=df_total,palette='Set2',ax=ax,s=65,order=['native','best','topscore_1'])
      #plt.savefig('totalenergies_important_AMAMut3_%s.png' %(col),transparent=True,dpi=1200)
      plt.show()
      plt.close()
      #exit()
  
  #return df
  '''


if __name__ == '__main__':
  pose_mod = Pose()
  pose_ref = Pose()
  #pose_mod = pose_from_pdb(sys.argv[1])
  #pose_ref = pose_from_pdb(sys.argv[2])
  #rmsd =CalculateRMSDCA(pose_mod,pose_ref)
  #print(rmsd)
  #files = ['session_JeffRuffolo_AMAmutant3_bestsampled.pdb',
  #        'session_JeffRuffolo_AMAmutant3_native.pdb',
  #        'session_JeffRuffolo_AMAmutant3_topscore_1.pdb']
  files = ['PackandRelax_498ALA_500HIS/unglycosylated/T2_with_UDPGalNAc_withS_from_4d0z_ATTAAPRC_498ALA_500HIS.pk_1251.pdb.gz',
  'PackandRelax_498THR_500PRO/unglycosylated/T2_with_UDPGalNAc_withS_from_4d0z_ATTAAPRC_498THR_500PRO.pk_1935.pdb.gz',
  'PackandRelax_498PRO_500PRO/unglycosylated/T2_with_UDPGalNAc_withS_from_4d0z_ATTAAPRC_498PRO_500PRO.pk_1099.pdb.gz',
  'PackandRelax_498THR_500GLN/unglycosylated/T2_with_UDPGalNAc_withS_from_4d0z_ATTAAPRC_498THR_500GLN.pk_1705.pdb.gz']
  tags=['ALA,HIS','THR,PRO','PRO,PRO','THR,GLN']
  #EnergyMapForPdbFiles(files,tags,suffix='alahis',residues=[501,502,503,504])
  for ifile,cfile in enumerate(files[:1]):
    dict_ = PairwiseEnergyMapForPdbFile(cfile,tag=tags[ifile],residues=[498,500])
  tidy_df = pd.DataFrame()
  for key in dict_:
    tidy_df[key] = dict_[key]
