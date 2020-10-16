import os
import sys
colors=['aquamarine','paleyellow','lightblue','salmon','limon','pink','lightorange']
#color_surface = 'praseodymium'
color_surface = 'white'
color_cartoon = 'scandium'
color_HWU = 'boron'
sele_Thr_glyc = "/%s//P/THR`3"
surface_transparency=0.1
maxload=1
savefig=True

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
      argname_k207 = curprot + " and resi 207 "
      cmd.show("sticks",argname_k207)
      cmd.color("lightpink",argname_k207)
      cmd.util.cnc(argname_k207)
      cmd.color('lightpink',argname_k207+' and hydrogens')
      cmd.hide("sticks","hydrogens")
      cmd.util.cnc(argname_P + ' and resi 4 ')
      #cmd.set("transparency",0.5,argname_k207)
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
cmd.set('ambient_occlusion_scale',11)
cmd.util.ray_shadows('occlusion2')
cmd.set('stick_radius',0.31, 'resname HWU')
cmd.set('stick_radius',0.25, 'chain P')
if savefig:
  cmd.set('surface_quality',2)
cmd.set('spec_reflect',0.2)
cmd.set_view([0.990826607,   -0.128174603,    0.042645268,
     0.107263848,    0.938434541,    0.328366935,
    -0.082108788,   -0.320781440,    0.943583190,
     0.000177413,   -0.000515442,  -49.219161987,
   -17.985958099,   36.096130371,  -10.925792694,
    -0.336019993,   98.884994507,  -20.000000000])
cmd.scene('001','store')
#cmd.set('transparency_mode',1)
#cmd.set('ray_transparency_oblique')
#cmd.set('ray_transparency_oblique_power',5)
#cmd.set('ray_transparency_contrast',5)
#cmd.set('ray_trace_mode',1)
#cmd.set_view([0.981157362,   -0.190879568,    0.029653054,
#     0.158735752,    0.884184420,    0.439332157,
#    -0.110079117,   -0.426348925,    0.897831857,
#     0.000190074,   -0.000424286,  -45.424575806,
#   -17.730279922,   35.478080750,  -10.630357742,
#    -4.129285336,   95.091735840,  -20.000000000])
cmd.set('ray_trace_mode',1)
cmd.set('ray_trace_gain',0.1)
cmd.set('ray_trace_disco_factor',1)
cmd.set('ray_shadows','off')
if savefig:
  cmd.ray(800,600)
  cmd.png(sys.argv[-1])
