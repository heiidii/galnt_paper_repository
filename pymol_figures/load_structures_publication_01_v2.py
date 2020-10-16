import os
import sys
#from pymol import stored
#scandium
colors=['aquamarine','paleyellow','lightblue','salmon','limon','pink','lightorange']
#color_surface = 'praseodymium'
color_surface = 'white'
color_cartoon = 'scandium'
color_HWU = 'boron'
sele_Thr_glyc = "/%s//P/THR`3"
surface_transparency=0.0
maxload=1
pocket_residues = ['257', '289','290','291','288']
pocket_residues_string = ' and resi ' + ' or resi '.join(pocket_residues)
pocket_residues_string = ' ( '+pocket_residues_string+' ) '
savefig=False

def openfromfile(filelist):
  argname_base ='T2_*'
  for i,tempfilename in enumerate(filelist):
   print(tempfilename)
   f = open(tempfilename,'r')
   lines=f.readlines()
   f.close()
   color_id = colors[i]
   for line in lines[:(min(maxload,len(lines)))]:
      curfile = line.split()[0]
      print(curfile)
      cmd.load(curfile)
      if curfile.find('.gz') != -1:
      	argname_base = curfile.split('/')[-1].split(".pdb.gz")[0]
      else:
	      argname_base = curfile.split('/')[-1].split(".pdb")[0]
      argname = "/%s//*" %argname_base
      print("ARGNAME",argname)
      argname_AB = "/%s//A or /%s//B" %(argname_base,argname_base)
      print("A:",argname_AB)
      cmd.hide('cartoon',argname_AB)
      argname_P = "/%s//P" %(argname_base)
      cmd.select(argname_P)
      cmd.show_as('sticks',"sele")
      cmd.hide('cartoon',"sele")
      cmd.hide('sticks','hydrogens')
      cmd.color(color_id,"sele")
      argname_prot_noudp = "%s and chain A and (not resname HWU)" %(argname_base)
      curprot = 'prot_%d' %i
      cmd.create(curprot,argname_prot_noudp)
      cmd.color(color_surface,curprot)
      cmd.show_as('surface',curprot)
      cmd.set('transparency',surface_transparency,curprot)
      argname_HWU = "resname HWU"
      cmd.select(argname_HWU)
      cmd.color('orange',"sele")
      argname_P_subset = argname_P + ' and ( resi 7 or resi 8)'
      cmd.hide('sticks',argname_P_subset)
      #cmd.color('red',argname_HWU + ' and name O1B ')
      #cmd.show('spheres',argname_HWU + ' and name O1B ')
      #cmd.color('blue',argname_P + ' and (resi 3) and name N ')
      #cmd.show('spheres',argname_P + ' and (resi 3) and name N ')

  cmd.bg_color(color="white")
  cmd.set('sphere_scale',0.36)
  return curprot
  
curprot = openfromfile([sys.argv[1]])
cmd.hide('spheres','chain B')
cmd.set('two_sided_lighting')
#cmd.set('ambient_occlusion_scale',11)
#cmd.util.ray_shadows('occlusion2')
cmd.set('stick_radius',0.31, 'resname HWU')
cmd.set('stick_radius',0.3, 'chain P')
if savefig:
  cmd.set('surface_quality',2)
cmd.set('spec_reflect',0.2)
cmd.set_view([0.990826607,   -0.128174603,    0.042645268,
     0.107263848,    0.938434541,    0.328366935,
    -0.082108788,   -0.320781440,    0.943583190,
     0.000177413,   -0.000515442,  -49.219161987,
   -17.985958099,   36.096130371,  -20.925792694,
    -0.336019993,   98.884994507,  -20.000000000])
cmd.scene('001','store')
cmd.set_view([0.981157362,   -0.190879568,    0.029653054,
     0.158735752,    0.884184420,    0.439332157,
    -0.110079117,   -0.426348925,    0.897831857,
     0.000190074,   -0.000424286,  -45.424575806,
   -17.730279922,   35.478080750,  -20.630357742,
    -4.129285336,   95.091735840,  -20.000000000])
cmd.scene('002','store')
cmd.set('ray_trace_mode',1)
cmd.set('ray_trace_gain',0.1)
cmd.set('ray_trace_disco_factor',1)
cmd.set('ray_shadows','off')
if savefig:
  cmd.ray(800,600)
  cmd.png(sys.argv[-2])

cmd.set_view([0.504837930,    0.002740630,   -0.863198340,
     0.243603423,    0.958888113,    0.145517126,
     0.828116775,   -0.283741266,    0.483424038,
     0.000439802,   -0.000729280,  -51.568401337,
   -19.725727081,   38.694969177,  -21.041263580,
    -3.494707108,  106.578361511,  -20.000000000])

argname_pocket  = curprot + pocket_residues_string
cmd.color('lightpink',argname_pocket)
cmd.show('sticks',argname_pocket)
cmd.set('stick_radius',0.35, argname_pocket)
#cmd.select(argname_pocket)
cmd.hide('sticks','hydrogens')
#cmd.show('cartoon',curprot)
#cmd.set('cartoon_transparency',0.,curprot)
cmd.set('stick_radius',0.2, 'resname HWU')
cmd.util.cnc(argname_pocket)
cmd.color('lightpink',argname_pocket+' and hydrogens')
#cmd.color('blue',argname_pocket+' and name N')
#cmd.color('red',argname_pocket+' and name C')
cmd.set('transparency',0.1,argname_pocket)
cmd.scene('003','store')
if savefig:
  cmd.ray(800,600)
  cmd.png(sys.argv[-1])
#cmd.save(sys.argv[-1]+'.pse',state=0)
