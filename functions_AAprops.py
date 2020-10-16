import os
import sys
from pyrosetta import *
from rosetta import *
pyrosetta.init('-include_sugars -mute core')
from rosetta.core.scoring import *
from rosetta.protocols.docking import *

#Setting up mover for viewing structures
import argparse
from rosetta.core.pose import *
import math

PI = math.pi
calculatermsd=False
scriptmode=False
verbose=False
dihedrals = ["phi","psi","omega"]

atomdict = {'HWUO1B-THR7N' : ['O1B' , 'N'] , 'HWUO1-THR7N' : ['O1\'' , 'N'], 'UDP2OPB-THR7N' : ['2OPB','N'] , 'UDP3OPB-THR7N' : ['3OPB', 'N'] }



distdict = dict()


inputfile=''
if scriptmode:
	inputfile = sys.argv[1]
	execfile(inputfile)


def getdistanceforpose_and_atoms(pose,res1,res2,a1,a2):
  a1 = pose.residue(res1).xyz(a1)
  a2 = pose.residue(res2).xyz(a2)
  dispvec = a1 - a2
  return dispvec.norm()

def get_resids_for_pose(pose):
  for i in range(pose.total_residue()-1,1,-1):
    print(pose.residue(i).name())
    if name.find('HWU'):
      udp = i
    
def getdistanceforfile_OHsidechain_Obackbone_bond(fname,pep_r=498,enz_r=288,enz_a='O'):
    pose = pose_from_pdb(fname)
    resname = pose.residue(pep_r).name()
    if resname == 'THR':
      pep_a = 'OG1'
    elif resname == 'SER':
      pep_a = 'OG'
    else:
      return None
    d = getdistanceforpose_and_atoms(pose,pep_r,enz_r,pep_a,enz_a)
    return {'dist_%dOG-%d%s' %(pep_r,enz_r,enz_a) :d}

def getdistanceforfile_plusone(fname,pep_r=500,enz_r=287,enz_a='CB',pep_a='CA'):
    pose = pose_from_pdb(fname)
    d = getdistanceforpose_and_atoms(pose,pep_r,enz_r,pep_a,enz_a)
    return {'dist_%d%s-%d%s' %(pep_r,pep_a,enz_r,enz_a) :d}  

def getdistancemetricforfile(fname,peptide_residues=None,udpres=None):
  #print(fname,peptide_residues,udpres)
  pose = pose_from_pdb(fname)
  udp_atoms = ['C1\'','O1\'','O1B']
  udp_names = ['C1','O1','O1B']
  backbone_atoms = ['N','C','CA']
  if peptide_residues is None:
    peptide_residues = [498,499,500]
    #peptide_residues, udpres = get_resiids_for_pose(pose)
  if udpres is None:
    udpres = 495
  distdict ={}
  tupledict = {}
  for iudpa,udpa in enumerate(udp_atoms):
    distsum = 0.0
    distsum_CA = 0.0
    for ip,pepres in enumerate(peptide_residues):
      for bka in backbone_atoms:
        tempkey =  'd_pep%d%s_udp_%s' %(ip-1 , bka , udp_names[ iudpa ])
        temptuple = (ip-1 , bka , udp_names[ iudpa ])

        tupledict[tempkey]=temptuple
        distdict[tempkey]= getdistanceforpose_and_atoms(pose,udpres,pepres,udpa,bka)

        distsum +=  distdict[tempkey]
      distsum_CA += getdistanceforpose_and_atoms(pose,udpres,pepres,udpa,'CA') 

    distdict[ 'distance_bk_udp_' + udp_names[ iudpa ] ] = distsum
    temptuple = ('9','bk',udp_names[ iudpa ])
    tupledict[ 'distance_bk_udp_' + udp_names[ iudpa ] ]= temptuple
    
    temptuple = ('3','CA',udp_names[ iudpa ])
    distdict[ 'distance_ca_udp_' + udp_names[ iudpa ] ] = distsum_CA
    tupledict[ 'distance_ca_udp_' + udp_names[ iudpa ] ]= temptuple

  return distdict,tupledict

def getdistanceforpose(pose):
        for entry in resdict:
                atom1 = pose.residue(resdict[entry][0]).xyz(atomdict[entry][0])
                atom2 = pose.residue(resdict[entry][1]).xyz(atomdict[entry][1])
                dispvec = atom1 - atom2
                dist = dispvec.norm()
                distdict[entry] = dist
        return distdict

def getdistances(files):
        propertiesdict = dict()
        for entry in resdict:
          propertiesdict[entry]=[]
        for f in files:
                pose = pose_from_pdb(f)
                dist = getdistanceforpose(pose)
                for entry in dist:
			              propertiesdict[entry].append(dist[entry])
        return propertiesdict

def getdihedralforpose(pose,residue):
        phi = pose.phi(residue)
        psi = pose.psi(residue)
        omega = pose.omega(residue)
        return phi,psi,omega


			
def CalculateRMSDCA(pose_mod,pose_ref):
  rmsd = CA_rmsd(pose_mod,pose_ref)
  return rmsd

def CalculateSASA(pose):
  sasaC = core.scoring.sasa.SasaCalc()
  totalsasa = sasaC.calculate(pose) 
  v_ressasa = sasaC.get_residue_sasa()
  l_ressasa = list(v_ressasa)
  print(len(l_ressasa))
  for i,entry in enumerate(l_ressasa):
    print(i+1,pose.residue(i+1).name(),entry)
  

def get_sc_for_pose(pose,molecule_1='A',molecule_2='P'):
  scC = core.scoring.sc.ShapeComplementarityCalculator()
  for ires in range(1,pose.size()+1):
    chain =  pose.pdb_info().chain(ires)
    if chain==molecule_2:
      scC.AddResidue(1, pose.residue(ires))
    elif chain==molecule_1:# or chain=='B':
      scC.AddResidue(0, pose.residue(ires))
    else:
      continue
  ran = scC.Calc()
  results = scC.GetResults()
  return results

def get_sc_for_file(infile,molecule_1='A',molecule_2='P',removeidentity=False,residue_remove=None):
  pose = Pose()
  try:
    pose = pose_from_pdb(infile)
  except RuntimeError:
    return {}
  #if removeidentity:
    
  sc = get_sc_for_pose(pose,molecule_1=molecule_1,molecule_2=molecule_2)
  dict_ = {}
  dict_['sc_shapecomplementarity']=sc.sc
  dict_['sc_area'] = sc.area
  dict_['sc_distance']=sc.distance
  return dict_

def get_fnat_and_fnotnat(file1,file2='T2_with_UDPGalNAc_withS_from_4d0z_ATTAAPRC_498THR_500PRO.pk.pdb'):
  pose_ref = pose_from_pdb(file2)
  pose = pose_from_pdb(file1)
  sf = ScoreFunction()
  sf = get_fa_scorefxn()
  setup_foldtree(pose_ref,'AB_P', Vector1([1]))
  setup_foldtree(pose,'AB_P', Vector1([1]))
  fnat = pyrosetta.rosetta.protocols.docking.calc_Fnat(pose,pose_ref,sf,Vector1([1]))
  fnonnat = pyrosetta.rosetta.protocols.docking.calc_Fnat(pose,pose_ref,sf,Vector1([1]))
  dict_ = {}
  dict_['fnat']=fnat
  dict_['fnonnat']=fnonnat
  return dict_

  
