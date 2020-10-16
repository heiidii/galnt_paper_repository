import os
import sys
#from pymol import stored
#scandium
colors=['scandium','aquamarine','paleyellow','lightblue','salmon','limon','pink','lightorange']
#color_surface = 'praseodymium'
color_surface = 'white'
color_cartoon = 'scandium'
color_HWU = 'boron'
sele_Thr_glyc = "/%s//P/THR`3"
surface_transparency=0.0
maxload=1
pocket_residues = ['257', '289','290','291']
pocket_residues_string = ' and resi ' + ' or resi '.join(pocket_residues)
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
      if i==0:
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
  
curprot = openfromfile(sys.argv[1:3])
cmd.hide('spheres','chain B')
cmd.set('two_sided_lighting')
cmd.set('ambient_occlusion_scale',11)
cmd.util.ray_shadows('occlusion2')
cmd.set('stick_radius',0.31, 'resname HWU')
cmd.set('stick_radius',0.25, 'chain P')
if savefig:
  cmd.set('surface_quality',2)
cmd.set('spec_reflect',0.2)
cmd.set_view([0.981157362,   -0.190879568,    0.029653054,
     0.158735752,    0.884184420,    0.439332157,
    -0.110079117,   -0.426348925,    0.897831857,
     0.000190074,   -0.000424286,  -45.424575806,
   -17.730279922,   35.478080750,  -20.630357742,
    -4.129285336,   95.091735840,  -20.000000000])
cmd.set_view([0.628397942,   -0.065572485,   -0.775111496,
     0.270620376,    0.952619731,    0.138809636,
     0.729291201,   -0.296990186,    0.616377652,
     0.000196302,   -0.000449082,  -51.915107727,
   -18.331748962,   37.791656494,  -21.735721588,
     2.367432594,  101.588432312,  -20.000000000])
cmd.scene('002','store')
#cmd.set('transparency_mode',1)
#cmd.set('ray_transparency_oblique')
#cmd.set('ray_transparency_oblique_power',5)
#cmd.set('ray_transparency_contrast',5)
#cmd.set('ray_trace_mode',1)
cmd.set('ray_trace_mode',1)
cmd.set('ray_trace_gain',0.1)
cmd.set('ray_trace_disco_factor',1)
cmd.set('ray_shadows','off')
if savefig:
  cmd.ray(800,600)
  cmd.png(sys.argv[-3])
cmd.set_view([0.406665862,   -0.098616935,   -0.908228874,
     0.356439292,    0.932487488,    0.058348548,
     0.841167390,   -0.347459048,    0.414368480,
     0.000196302,   -0.000449082,  -50.182567596,
   -18.331748962,   37.791656494,  -21.735721588,
     0.634902120,   99.855918884,  -20.000000000])
argname_pocket  = 'prot_0 ' + pocket_residues_string
cmd.color('lightpink',argname_pocket)
cmd.scene('003','store')
if savefig:
  cmd.ray(800,600)
  cmd.png(sys.argv[-2])
cmd.set_view([-0.096462026,   -0.586937845,    0.803848684,
    -0.055954624,    0.809540153,    0.584380448,
    -0.993753493,    0.011392793,   -0.110933736,
    -0.000344856,   -0.000557811,  -54.301692963,
   -17.228042603,   34.733673096,  -23.525236130,
    -0.208567888,  108.934570312,  -20.000000000])
cmd.hide('surface',curprot)
cmd.util.cnc("resname HWU")
cmd.util.cnc("chain P and (resi 2 or resi 3)")
cmd.set('stick_radius',0.23, 'resname HWU')
cmd.set('ray_trace_gain',0.3)
cmd.scene('004','store')
if savefig:
  cmd.ray(800,600)
  cmd.png(sys.argv[-1])
