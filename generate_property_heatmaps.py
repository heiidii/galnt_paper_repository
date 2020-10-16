import os
import sys
import logging
logging.basicConfig(filename='scatterplot.log',level=logging.DEBUG)
from functions_additional_heatmaps_paper import plot_heatmaps_from_pickled_dataframe, plot_heatmaps_from_pickled_dataframe_energies


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


if __name__ == '__main__':
  pickledfiles = 'PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_A*.p'
  pickledfiles_pointwise = 'PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_A[ASGT]T*.p'
  pickledfiles_energies ='PickleFiles/20200207/UnglycosylatedPeptideData_Filtered_ClusteredDBScan_Energies_A[AGSTP]T*.p' 
    

  rmsd_cutoff = 1.0
  reverse=False #First filter then get top N
  applyfraction=True
  for Ncur in [1,10]:
    for off in [0.10]:
      for criteria in ['rmsd_less_than']:
        for prop in ['sc_shapecomplementarity']:
          print(Ncur)
          plot_heatmaps_from_pickled_dataframe(pickledfiles = pickledfiles_pointwise,rmsd_cutoff=rmsd_cutoff,criteria=criteria,prop=prop,cmap='Greys',N=Ncur,off=off,reverse=reverse,annotate=False,applyfraction=applyfraction,listres1=['G','A','S','T'],show=False,positionwise = True,pointplot=False,heatmaps=False,boxplot=True,by='median')
          plot_heatmaps_from_pickled_dataframe(pickledfiles = pickledfiles,rmsd_cutoff=rmsd_cutoff,criteria=criteria,prop=prop,cmap='Greys',N=Ncur,off=off,reverse=reverse,annotate=False,applyfraction=applyfraction,show=False,positionwise = False,pointplot=False,heatmaps=True,boxplot=False,by='median',rocs=True,off_list=[0.10])
          