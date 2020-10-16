import os
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import glob
import mpl_toolkits
#from matplotlib.mlab import griddata
plt.style.use('classic')

from functions_parsefiles import *
import functions_AAprops

listnotfound = []

i0=1

def setupproperties(info):
   properties = []
   properties_residues=[]
   if "additionalprops" in info:
                                properties = info["additionalprops"]
   else:
                           properties.append("backbone_dihedrals")
   if "additionalprops_residues" in info:
                           properties_residues = info["additional_props_residues"]
   else:
                           properties_residues.append(498)
                           properties_residues.append(499)
                           properties_residues.append(500)
   return properties_residues, properties
      
def getadditionalproperties(info,top_files):
  properties_residues, properties = setupproperties(info)
  for prop in properties:
    if prop=='backbone_dihedrals':
                  additionalproperties = functions_AAprops.getdihedrals(top_files,properties_residues)
    elif prop=='distances':
      additionalproperties = functions_AAprops.getdistances(top_files)
    elif prop=='none':
      return None
    else:
      return None
  return additionalproperties


def filterpdb(sfile,basename,createlink,fields,outputdir):
  f = open(sfile,'r')
  lines = f.readlines()
  f.close()
  outf = open("FilteredFile_%s.txt" %basename,'w')
  outfdata = open("FilteredFileData_%s.txt" %basename,'w')
  mapheader = dict()
  histdict = dict()
  mapheader = getmappingfromlines(lines)
  for field in fields:
                minfield = fields[field][0]
                maxfield = fields[field][1]
                outf.write("#%s\t%f\t%f\n" %(field, minfield, maxfield))
                if field in mapheader:
                  fields[field].append(mapheader[field])
                  histdict[field]=[]
  for line in lines[i0+1:]:
    data = line.split()
    if line.find("I_sc") != -1: continue
    setTrue = False
    addstring = []
    tempval = dict()
    for field in fields:
      minfield = fields[field][0]
      maxfield = fields[field][1]
      index = fields[field][2]
      val = float(data[index])
      if val<maxfield and val>=minfield:
        setTrue = True
        tempval[field]=val
        addstring.append(str(val))
      else:
        setTrue = False
        break
    if setTrue:
      pdbfilename=data[mapheader["description"]]
      for field in fields:
         histdict[field].append(tempval[field])
      
      outf.write(pdbfilename+'\t'+'\t'.join(addstring)+"\n")
      outfdata.write(line)
      if createlink:
        fullpath = outputdir + "/"+pdbfilename
        #print pdbfilename,filtereddir
        parentfilename = pdbfilename[:-6]
        #print parentfilename
        if os.path.exists(outputdir + "/"+pdbfilename+'.pdb.gz'):
          cmd="gunzip %s.pdb.gz" %fullpath
          os.system(cmd)
        if os.path.exists(outputdir + "/"+pdbfilename+".pdb"):
          cmd='cp %s/%s.pdb %s/.' %(outputdir,pdbfilename,filtereddir)
          os.system(cmd)
        if os.path.exists(pdbfilename+".pdb"):
          cmd = 'cp %s.pdb %s/.' %(pdbfilename,filtereddir)
          os.system(cmd)
        if os.path.exists(parentfilename+".pdb"):
                                        cmd = 'cp %s.pdb %s/.' %(parentfilename,filtereddir)
                                        os.system(cmd)
        
  outf.close()
  outfdata.close()

def mergefiles(scorefile1,scorefile2,type="parentchild"):
  mapping1 = dict()
  mapping2 = dict()
  mapping1 = getmapping(scorefile1)
  mapping2 = getmapping(scorefile2)
  sf1 = open(scorefile1,'r')
  if type=="parentchild":
    outf=open("Mergedfile.sc",'w')
    newfields=""
    for field in mergedata:
          newfields += field + '\t'
          sf11=open(scorefile1,'r')
          headerlines1 = sf11.readlines()[:i0+1]
          sf11.close()
          header = headerlines1[i0].split('\n')[0] + '\t' +newfields + '\n'
          outf.write(headerlines1[0]+header)
    i=0
    for line1 in sf1:
      if line1.find("description")!=-1 or line1.find("SEQUENCE")!=-1:
        i+=1
        continue
      data1 = line1.split()
      filename1 = data1[mapping1["description"]]
      i+=1
      sf2 = open(scorefile2,'r')
      for line2 in sf2:
        if line2.find("description")!=-1 or line2.find("SEQUENCE")!=-1:
                                  continue
        data2 = line2.split()
        filename2 = data2[mapping2["description"]]
        if filename1.find(filename2):
          outline = line1.split('\n')[0]
          for field in mergedata:
            outline += '\t' + data2[mapping2[field]]
          outf.write(outline+'\n')
          break
      sf2.close()
    outf.close()
      
  sf1.close()

def dictfromlines(lines,fields=None):
  mapheader = dict()
  mapheader = getmappingfromlines(lines,"quiet")
  dictofarrays = dict()
  delete_fields =[]
  max_index = 1
  if fields is None:
    for key in mapheader:
        fields[key]=[-20000,20000,0]
  del fields['description'] #max min does not apply
  if not fields is None:
    for field in fields:
                dictofarrays[field]=[]
                if field in mapheader:
                        fields[field][2] = mapheader[field]
                        if fields[field][2] > max_index:
                                max_index = fields[field][2]
                else:
                        delete_fields.append(field)

    for field in delete_fields:
                del fields[ field ]

  dictofarrays["description"]=[]

  for line in lines[i0+1:]:
                data = line.split()
                if line.find("total_score") != -1: continue
                setTrue = False
                addstring = []
                if len(data) < max_index: continue
                tempofarrays = dict()
                if "description" in mapheader:
                   tempofarrays["description"]=data[mapheader["description"]]
                for field in fields:
                        minfield = fields[field][0]
                        maxfield = fields[field][1]
                        index = fields[field][2]
                        val = float(data[index])
                        if val<maxfield and val>=minfield:
                                setTrue = True
                                addstring.append(str(val))
                                tempofarrays[field] = val

                        else:
                                setTrue = False
                                break
                if setTrue:
                        for field in fields:
                                dictofarrays[field].append(tempofarrays[field])
                        if "description" in tempofarrays:
                                dictofarrays["description"].append(tempofarrays["description"])
  return dictofarrays, fields

def serialize_scorefile(infile,outfile,fields=None,serialize=True):
  f = open(infile,'r')
  lines = f.readlines()
  f.close()
  dictofarrays , fields = dictfromlines(lines,fields=fields)
  if serialize:
    import pickle
    pickle.dump(dictofarrays,open(f,'wb')) 
  return dictofarrays, fields

def getfiltereddatafromlines(lines,fields_local=None,dictofarrays=None,nodelete=True):
  mapheader = dict()
  mapheader = getmappingfromlines(lines,"quiet")
  dictofarrays = dict()
  fields_avail = {}
  delete_fields =[]
  max_index = 1
  for field in fields_local:
                dictofarrays[field]=[]
                minfield = fields_local[field][0]
                maxfield = fields_local[field][1]
                if field in mapheader:
                        fields_local[field][2] = mapheader[field]
                        if fields_local[field][2] > max_index:
                                max_index = fields_local[field][2]
                        fields_avail[field] = fields_local[field]
                else:
                        delete_fields.append(field)
  
  dictofarrays["description"]=[]  
  for line in lines[i0+1:]:
                data = line.split()
                if line.find("total_score") != -1: continue
                setTrue = False
                addstring = []
                if len(data) < max_index: continue
                tempofarrays = dict()
                if "description" in mapheader:
                   tempofarrays["description"]=data[mapheader["description"]]
                for field in fields_avail:
                        minfield = fields_avail[field][0]
                        maxfield = fields_avail[field][1]
                        index = fields_avail[field][2]
                        val = float(data[index])
                        if val<maxfield and val>=minfield:
                                setTrue = True
                                addstring.append(str(val))
                                tempofarrays[field] = val
        
                        else:
                                setTrue = False
                                break
                if setTrue:
                        for field in fields_avail:
                                dictofarrays[field].append(tempofarrays[field])
                        if "description" in tempofarrays:
                                dictofarrays["description"].append(tempofarrays["description"])
  return dictofarrays, fields_avail

def readfiltereddatafromlines(lines,fields,dictofarrays):
        mapheader = dict()
        mapheader = getmappingfromlines(lines,"quiet")
        dictofarrays = dict()
        delete_fields =[]
        max_index = 1
        for field in fields:
                #if field=='description': continue
                dictofarrays[field]=[]
                if field in mapheader:
                        fields[field][2] = mapheader[field]
                        if fields[field][2] > max_index:
                                max_index = fields[field][2]
                else:
                        delete_fields.append(field)

        for field in delete_fields:
                del fields[ field ]

        dictofarrays["description"]=[]

        for line in lines[i0+1:]:
                data = line.split()
                if line.find("total_score") != -1: continue
                setTrue = False
                addstring = []
                if len(data) < max_index: continue
                tempofarrays = dict()
                if "description" in mapheader:
                        tempofarrays["description"]=data[mapheader["description"]]
                for field in fields:
                        index = fields[field][2]
                        val = float(data[index])
                        dictofarrays[field].append(val)
                        if "description" in tempofarrays:
                                dictofarrays["description"].append(tempofarrays["description"])
        return dictofarrays, fields

def getcleannames(filebasepath, files):
  newfiles=[]
  for fname in files:
    newname = filebasepath + '/'+fname + '.pdb.gz'
    if os.path.exists(newname):
      newfiles.append(newname)
    else:
      temp = fname.split('.pdb')[0][:-5]
      tempname = filebasepath + '/'+ temp + '.pdb'
      if os.path.exists(tempname):
        newfiles.append(tempname)
        continue
      elif os.path.exists(tempname+'.gz'):
        newfiles.append(tempname+'.gz')
      else:
        newfiles.append('')
        exit()
  assert(len(newfiles)==len(files))
  return newfiles


def getsorteddatafromdict(fields, dictofarrays, sortfield):
 sorteddict = dict()
 sortedindices = np.argsort(dictofarrays[sortfield])
 #fields["description"]=[]
 for field in fields:
  sorteddict[field] = []
 
 for ind in sortedindices:
  for field in fields:
    sorteddict[field].append(dictofarrays[field][ind])

 return sorteddict

def getTopN(info):
 fields_TopN = info["fields"]
 for sortfield in info["sortfield"]:
  for num_N in info["N"]:
   outfvals = open("Stats_%s_%s_Top%03d.txt" %(info["outtag"],sortfield,num_N),'w')
   fieldstring = "\t"
   for field in fields_TopN:
     fieldstring += "%s\t%f - %f ;" %(field,fields_TopN[field][0],fields_TopN[field][1])
     outfvals.write("#%s\n" %fieldstring)
     for ifile,sfile in enumerate(info["files"]): #Scorefiles - one for peptide
        if not os.path.exists(sfile):
           listnotfound.append(sfile)
           continue
        f = open(sfile,'r')
        lines = f.readlines()
        f.close()
        dictofarrays = dict()
        unsorteddictofarrays , fields_TopN = getfiltereddatafromlines(lines,fields_TopN)
        dictofarrays = getsorteddatafromdict(fields_TopN,unsorteddictofarrays, sortfield)
        additionalproperties = dict()
        if "calculate_additionalprop" in info:
           sortedkeys = []
           if info["calculate_additionalprop"]:
              #print dictofarrays["description"][:num_N]  
              #Due to rescoring - sometimes same file repeats
              top_files = getcleannames(info['files_pdb'][ifile],dictofarrays["description"][:num_N]) #pdbfiles ~ Ndecoys per scorefile ~ 2000
      
              additionalproperties = getadditionalproperties(info,top_files)
              if not additionalproperties is None:
                   sortedkeys= sorted(additionalproperties)
              if "save_additionalprop" in info:
                if info["save_additionalprop"]:
                  outfap = open("AdditionalProps_%s_%s_%s_Top%03d.txt" %(info['keys'][ifile],info["outtag"],sortfield,num_N),'w') 
      
                  curind = 0
                  #header = "SEQUENCE:\nSCORE: " + sortfield
                  header = ''
                  if 'distance_catalysis' in fields_TopN:
                    header  = "SEQUENCE:\nSCORE: interaction_energy distance_catalysis substrate_ca_rmsd"
                  else:
                    header  = "SEQUENCE:\nSCORE: interaction_energy substrate_ca_rmsd"  
                  for okey in sortedkeys:
                    header += " "+okey
                  header += " description"
                  outfap.write(header + '\n' )
                  for tfile in top_files:
                    outstring =["SCORE: "] #pdbfilename
                    outstring.append(str(dictofarrays["interaction_energy"][curind]))
                    if 'distance_catalysis' in fields_TopN:
                      outstring.append(str(dictofarrays["distance_catalysis"][curind]))
                    outstring.append(str(dictofarrays["substrate_ca_rmsd"][curind]))
                    for okey in sortedkeys:
                      outstring.append(str(additionalproperties[okey][curind]))
                      outstring.append(tfile)
                      #print outstring
                    outfap.write('\t'.join(outstring)+'\n')
                    curind+=1
                  outfap.close()
        
        if "save_pdbfilenames" in info:
          if info["save_pdbfilenames"]:
            outfpdb = open("PDBFiles_%s_%s_%s_Top%03d.txt" %(info['keys'][ifile],info["outtag"],sortfield,num_N),'w')
            for tfile in top_files:
              outfpdb.write(tfile+'\n')
            outfpdb.close()       
        for apkey in additionalproperties:
          #print apkey
          fields_TopN[apkey]=[]
          dictofarrays[apkey]=additionalproperties[apkey]
          #print apkey, additionalproperties[apkey]
          outstring = ""
          if 'description' in fields_TopN:
            del fields_TopN[ 'description' ]
          for field in fields_TopN:
            #if field=='description': continue
            avg = 0.0
            sum_ = 0.0
            sum_sq = 0.0
            sd = 0.0
            icount = 0
            for val in dictofarrays[field][:num_N]:
              sum_ += val
              sum_sq += val*val
              icount += 1
      
            if icount > 0:
                  avg = sum_/float(icount)
                  var = (sum_sq / float(icount) - avg * avg)
                  sd = var ** 0.5
            else:
              if field in info['defaults_avg']:
                avg = info["defaults_avg"][field]
              else:
                avg = -999
            outstring +=  "%s\t%f\t%f\t%s\n" %(field,avg,sd,info['keys'][ifile])
            outfvals.write(outstring)
        for apkey in additionalproperties:
            del fields_TopN[ apkey ] #before starting next file
   outfvals.close()

