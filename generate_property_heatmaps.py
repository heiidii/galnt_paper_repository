import os
import sys
import logging
from MultivariablePlotting import *
from DataSortingFiltering import getdata_xyzc
logging.basicConfig(filename='scatterplot.log',level=logging.DEBUG)


def get_files(pattern,size=0.520):
  import glob
  raw_files = glob.glob(pattern)
  print(len(raw_files))
  files=[]
  for rfile in raw_files:
      if os.path.exists(rfile):
        s = os.path.getsize(rfile) #s in bytes
        print(rfile,s/(1000*1000.0),"mb")
        if s/(1000.0*1000.0) >size:
          files.append(rfile)
  print( "Files:",'\n'.join(files))
  print("TOTALFILES: ", len(files))
  return files


def get_combined_plots_all_categories(all_categories=['all'],all_expdata=['on','medium','high','off'],filterbyfield={'interaction_energy':'high'},show=True,residues=['all'],cluster=True,clusterOnly=True):
  pattern = 'PackandRelax_498*/unglycosylated/score.sc'
  files = get_files(pattern )
  for residue in residues:
    
    if residue != 'all':
      pattern = 'PackandRelax_498%s*/unglycosylated/score.sc' %residue
      files = get_files(pattern )

    for category in all_categories:
      for expdata in all_expdata:
        for key in filterbyfield:
          key_all=''
          suffix=''

          if residue != 'all':
            key_all = 'X-1 %s\n'%residue
            suffix='X-1_%s_' %residue

          if category != 'all':
            key_all += category
            suffix += '%s' %category

          if expdata != 'all':
            key_all += 'exp '+expdata
            suffix += 'exp_%s_'%expdata

          if key != 'all':
            key_all += '\n%s %s'%(key,filterbyfield[key])
            suffix += '%s_%s'%(key,filterbyfield[key])
          
          if not clusterOnly:

            #plot3d_from_combined_dataframe(files,outfile=get_full_name(files[0],type_='3d_all',suffix=suffix),key_all=key_all,category=category,expdata=expdata,filterbyfield=filterbyfield,show=show)
            suffix += '_view3'
            plot3d_from_combined_dataframe(files,outfile=get_full_name(files[0],type_='3d_all',suffix=suffix),key_all=key_all,view=[15,-61],category=category,expdata=expdata,filterbyfield=filterbyfield,show=show)

          if cluster:
            suffix += '_clusters_view3'
            plot_clusters_from_combined_dataframe(files,outfile=get_full_name(files[0],type_='3d_all',suffix=suffix),key_all=key_all,view=[15,-61],category=category,expdata=expdata,filterbyfield=filterbyfield,show=show)

def run_plot3d_slice_topN(pattern,show=False):
  files = get_files(pattern)
  basedir_3d_compare = 'results/scatter3dTopN/'
  os.system('mkdir -p %s' %basedir_3d_compare)
  cluster=True
  for curfile in files:
    #for topN in topN_list:
    outfile = basedir_3d_compare + 'scatter3dTopN_'+info['xfield']+'_'+ get_outfile_name(curfile) +'.png'
    collection=[getdata_xyzc(curfile)]
    keys = [get_key(curfile)]
    plot3d_slice_topN(collection,keys,outfile=outfile,show=show,cluster=cluster,topNs=[100,200,500]) 


def get_for_files_from_commandline():
  collection =[]
  keys = []
  files = []
  if len(sys.argv)>1:
    file1 = sys.argv[1]
    files.append(file1)
    collection.append(getdata_xyzc(file1))
    keys.append(get_key(file1))
  if len(sys.argv)>2:
    file2 = sys.argv[2]
    files.append(file2)
    collection.append(getdata_xyzc(file2))
    keys.append(get_key(file2))
  if len(sys.argv)>3:
    file3 = sys.argv[3]
    files.append(file3)
    collection.append(getdata_xyzc(file3))
    keys.append(get_key(file3))
  plot_clusters_from_scorefiles(files,show=True,write=False)
  for i,dictdata in enumerate(collection):
    basedir = 'results/clusters/pdbfiles/'
    os.system('mkdir -p %s' %basedir)
    outfile = basedir + info['xfield']+'_'+ get_outfile_name(files[i])
    #write_sliced_cluster_pdbs(dictdata,thresholds=thresholds,outfile=outfile)
  if 'description' in collection[0]:
    print( collection[0]['description']['vals'][:10])

  return collection,keys,files


def plot3d_slice_topN_from_commandline(thresholds=[-20,-30,-35]):

  collection,keys,files = get_for_files_from_commandline()
  dim2d=False
  show=True
  cluster=False
  if dim2d:
    basedir_2d_compare = 'results/scatter2dcompare/'
    os.system('mkdir -p %s' %basedir_2d_compare)
    outfile = basedir_2d_compare + 'scatter2dslices_'+info['xfield']+'_'+ get_outfile_name_n(files) +'.png'
    plot3d_slice_and_compare_n(collection,keys,n=2,outfile=outfile,thresholds=[-32,-35,-37,-39],show=show,cluster=cluster,dim2d=dim2d,dim3d=(not dim2d))
  else:
    basedir_3d_compare = 'results/scatter3dTopN/'
    os.system('mkdir -p %s' %basedir_3d_compare)
    #outfile = basedir_3d_compare + 'scatter3dslices_'+info['xfield']+'_'+ get_outfile_name_n(files) +'.png'
    outfile = basedir_3d_compare + 'scatter3dTopN_'+info['xfield']+'_'+ get_outfile_name_n(files) +'.png'
    #plot3d_slice_topN(collection,keys,outfile=outfile,show=show,cluster=cluster,topNs=[100,200,500])
    print(thresholds)
    plot3d_slice_and_compare_n(collection,keys,n=2,outfile=outfile,thresholds=thresholds,show=show,cluster=cluster,dim2d=dim2d,dim3d=(not dim2d))


def plot3d_slices_from_files(pattern,show=False,thresholds=[-20],suffix='',cluster=True):
  files = get_files(pattern)
  for file1 in files:
    collection = [getdata_xyzc(file1)]
    keys = [get_key(file1)]
    basedir_3d_compare = 'results/scatter3dcompare/'
    os.system('mkdir -p %s' %basedir_3d_compare)
    outfile = basedir_3d_compare + 'scatter3dslices_'+info['xfield']+'_'+ get_outfile_name(file1) +'_%s.png' %suffix
    plot3d_slice_and_compare_n(collection,keys,n=2,outfile=outfile,thresholds=thresholds,show=show,cluster=cluster,dim2d=False,dim3d=True) 

if __name__ == '__main__':

  '''
  usage:
  files = get_files(pattern)
  run_plot3d_slice_topN(pattern)
  get_combined_plots_all_categories(all_expdata=['off','medium','on','high'],show=False,residues = ['ALA','THR','PRO','GLY','SER'],filterbyfield={'interaction_energy':'medium'})
  get_combined_plots_all_categories(all_expdata=['off'],show=True,filterbyfield={'interaction_energy':'high'})
  plot3d_from_scorefiles(files)
  plot_clusters_from_scorefiles(files,show=False,write=True)
  thresholds = [-20,-23,-25]
  suffix = 'low'
  show=False
  cluster = True
  write_sliced_cluster_pdbs_for_files(files,thresholds=thresholds)
  
  show=True
  files = []
  for res in ['ILE','GLU','THR']:
    pattern = 'PackandRelax_T12_YY_correct_531%s*/unglycosylated/score.sc' %res
    files += get_files(pattern)
  print( files)
  '''

  show = True
  pattern = 'PackandRelax_498GLU_500PRO/unglycosylated/score.sc'
  #files = get_files(pattern,size=0.480)
  thresholds = [-20]
  #write_sliced_cluster_pdbs_for_files(files,thresholds=thresholds,k_means=True)
  #plot_distribution_from_scorefiles(files,'yfield',violin=True,show=show,tag='',raw=False)
  #plot3d_slices_from_files(pattern,show=True,thresholds=thresholds,suffix='high')
  from MultivariablePlotting_Paper import plot_sinks_from_pickled_dataframe,plot_ramas
  for aa in ['Y']:#['P','G','S','T','V','A']:
    for Ncur in [2000]:
      print(aa)
      #plot_sinks_from_pickled_dataframe(clusterid=0,pickledfiles='PickleFiles/20200120/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_A%sTA*.p' %aa,basename='results_collated/scatterplots/',N=Ncur,x_max=2.5)
      #pickledfiles = 'PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_AAT*.p'
      #plot_heatmaps_from_pickled_dataframe(pickledfiles=pickledfiles,N=20,prop='sc_shapecomplementarity')
  #pickledfiles_aa = 'PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_A%sT*.p'
  #plot_ramas(pickledfiles_aa=pickledfiles_aa,aa=['S','A'],scatter=False,show=False,individual=True)
  #pickledfiles = 'PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_AAT*.p'
  #plot_sinks_from_pickled_dataframe(pickledfiles=pickledfiles)

  from functions_additional_heatmaps_paper import plot_heatmaps_from_pickled_dataframe, plot_heatmaps_from_pickled_dataframe_energies
  pickledfiles = 'PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_A*.p'
  pickledfiles_pointwise = 'PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_A[ASGT]T*.p'
  pickledfiles_energies ='PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_Energies_A[AGSTP]T*.p' 
  rmsd_cutoff = 1.0
  reverse=False #First filter then get top N
  applyfraction=True
  #maxmin = [0.75,0.80]
  for Ncur in [1,10]:#[1,10]:
    for off in [0.10]:
      for criteria in ['rmsd_less_than']:
        for prop in ['sc_shapecomplementarity']:
          print(Ncur)
          #plot_heatmaps_from_pickled_dataframe(pickledfiles = pickledfiles_pointwise,rmsd_cutoff=rmsd_cutoff,criteria=criteria,prop=prop,cmap='Greys',N=Ncur,off=off,reverse=reverse,annotate=False,applyfraction=applyfraction,listres1=['G','A','S','T'],show=False,positionwise = True,pointplot=False,heatmaps=False,boxplot=True,by='median')
          #plot_heatmaps_from_pickled_dataframe(pickledfiles = pickledfiles,rmsd_cutoff=rmsd_cutoff,criteria=criteria,prop=prop,cmap='Greys',N=Ncur,off=off,reverse=reverse,annotate=False,applyfraction=applyfraction,listres1=['G','A','S','T','P'],show=False,positionwise = False,pointplot=False,heatmaps=True,boxplot=False,by='median')
          plot_heatmaps_from_pickled_dataframe(pickledfiles = pickledfiles,rmsd_cutoff=rmsd_cutoff,criteria=criteria,prop=prop,cmap='Greys',N=Ncur,off=off,reverse=reverse,annotate=False,applyfraction=applyfraction,show=False,positionwise = False,pointplot=False,heatmaps=True,boxplot=False,by='median',rocs=True,off_list=[0.10])
          #plot_heatmaps_from_pickled_dataframe_energies(pickledfiles=pickledfiles_energies,basename='results/heatmaps_energies_pairwise_filtered/',N=Ncur,rmsd_cutoff=rmsd_cutoff,criteria=criteria,cmap=None,off=off,plotdist=False,fitstuff=False,reverse=False,maxmin=None,annotate=False,suffix=None,applyfraction=False,show=False,donotapplytopNtoparent=True,pairwise=True,listres1=['G','A','S','T','P'])
