import os
import sys
#from pymol import stored

colors=['lightblue','salmon','limon','pink','lightorange']
color_id_base = 'hydrogen'
sele_Thr_glyc = "/%s//C/THR`14"


def openfromfile(filelist):
  argname_base ='T12_*'
  for i,tempfilename in enumerate(filelist):
   print(tempfilename)
   f = open(tempfilename,'r')
   lines=f.readlines()
   f.close()
   color_id = colors[i]
   for line in lines[:(min(10,len(lines)))]:
      curfile = line.split()[0]
      print(curfile)
      cmd.load(curfile)
      argname_base = curfile.split('/')[-1].split(".pdb.gz")[0]
      argname = "/%s//*" %argname_base
      print("ARGNAME",argname)
      #cmd.select(argname)
      #cmd.color(color_id,"sele")
      argname_P = "/%s//C" %(argname_base)
      cmd.select(argname_P)
      cmd.show_as('sticks',"sele")
      cmd.hide('cartoon',"sele")
      cmd.hide('sticks','hydrogens')
      cmd.color(color_id,"sele")
      #cmd.align("/%s//A" %argname,"/1ACB_b//E")
      argname_AB = "/%s//A or /%s//B" %(argname_base,argname_base)
      print("A:",argname_AB)
      cmd.select(argname_AB)
      cmd.color(color_id_base,"sele")
      argname_HWU = "resname HWU"
      cmd.select(argname_HWU)
      cmd.color('orange',"sele")
      argname_e = argname_HWU+" and element o "
      cmd.select(argname_e)
      cmd.color('red',"sele")
      argname_e = argname_HWU+" and element n "
      cmd.select(argname_e)
      cmd.color('blue',"sele")
      argname_e = argname_HWU+" and element s "
      cmd.select(argname_e)
      cmd.color('yellow',"sele")
      argname_thr_glyc_e = sele_Thr_glyc %argname_base + " and element o "
      cmd.select(argname_thr_glyc_e)
      cmd.color('red',"sele")
      argname_thr_glyc_e = sele_Thr_glyc %argname_base + " and element n "
      cmd.select(argname_thr_glyc_e)
      cmd.color('blue',"sele")
      
  cmd.bg_color(color="white")
  
  
  #argname_HWU

openfromfile(sys.argv[2:])
