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
surface_transparency=0.10
surface_transparency_pep=0.50
maxload=1
pocket_residues = ['257', '289','290','291']
pocket_residues_string = ' and resi ' + ' or resi '.join(pocket_residues)
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
      #cmd.show_as('sticks',"sele")
      cmd.hide('cartoon',"sele")
      #cmd.hide('sticks','hydrogens')
      #cmd.color(color_id,"sele")
      curpep = 'pep_%d' %i
      cmd.create(curpep,argname_P)
      cmd.show_as('surface',curpep)
      argname_prot_noudp = "%s and chain A and (not resname HWU)" %(argname_base)
      curprot = 'prot_%d' %i
      cmd.create(curprot,argname_prot_noudp)
      cmd.color(color_surface,curprot)
      cmd.show_as('surface',curprot)
      cmd.set('transparency',surface_transparency,curprot)
      cmd.set('transparency',surface_transparency_pep,curpep)
      argname_HWU = "resname HWU"
      cmd.select(argname_HWU)
      cmd.color('orange',"sele")
      cmd.show_as('lines',argname_HWU)
      argname_P_subset_backbone = curpep + ' and ( name ca or name n or name c )'
      cmd.show('sticks',argname_P_subset_backbone)
      cmd.color(color_id,curpep)
      argname_P_imp_res = curpep + ' and (resi 2)'
      cmd.show('sticks',argname_P_imp_res)
      cmd.color(color_id,argname_P_imp_res)
      argname_his = 'chain A and resi 291'
      cmd.show('sticks',argname_his)
      cmd.hide('sticks','hydrogens')
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
#cmd.util.ray_shadows('occlusion2')
#cmd.set('stick_radius',0.31, 'resname HWU')
cmd.set('stick_radius',0.30, 'pep_0')
cmd.set('stick_radius',0.30, 'prot_0')
if savefig:
  cmd.set('surface_quality',2)
cmd.set('spec_reflect',0.2)
cmd.set('ray_trace_mode',1)
cmd.set('ray_trace_gain',0.1)
cmd.set('ray_trace_disco_factor',1)
cmd.set('ray_shadows','off')
cmd.set_view([0.972537398,   -0.124563880,   -0.196503386,
     0.151300848,    0.980212092,    0.127463728,
     0.176738292,   -0.153691411,    0.972172260,
     0.000861183,   -0.002124823,  -53.522701263,
   -21.330295563,   38.240524292,  -22.969825745,
    10.045980453,   97.567794800,  -20.000000000])

argname_pocket  = 'prot_0 ' + pocket_residues_string
cmd.color('lightpink',argname_pocket)
cmd.scene('003','store')
if savefig:
  cmd.ray(800,600)
  cmd.png(sys.argv[-1])

cmd.set_view ([0.767031491,   -0.012962798,   -0.641445220,
     0.254339725,    0.923999310,    0.285464734,
     0.589011252,   -0.382110775,    0.712062716,
     0.000560009,   -0.002080690,  -74.491249084,
   -14.611545563,   37.630973816,  -20.076419830,
    30.872917175,  118.394783020,  -20.000000000])

cmd.set('transparency',0,'pep_0')
cmd.set('transparency',0,'prot_0')
cmd.scene('004','store')
if savefig:
  cmd.ray(800,600)
  cmd.png(sys.argv[-1] + '_solidsurface')
