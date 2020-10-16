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
pocket_residues = ['257', '289','290','291']
pocket_residues_string = ' and resi ' + ' or resi '.join(pocket_residues)
distance=6.0
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
      cmd.hide('sticks','hydrogens')
      cmd.color(color_id,"sele")
      argname_prot_noudp = "%s and chain A and (not resname HWU)" %(argname_base)
      curprot = 'prot_%d' %i
      cmd.create(curprot,argname_prot_noudp)
      cmd.color(color_surface,curprot)
      cmd.hide('everything',curprot)
      #cmd.set('transparency',surface_transparency,curprot)
      argname_HWU = "resname HWU"
      cmd.hide('sticks',argname_HWU)
      cmd.color('orange',argname_HWU)
      argname_P_subset = argname_P + ' and ( resi 7 or resi 8)'
      cmd.hide('sticks',argname_P)
      #cmd.show('cartoon',curprot)
      #cmd.color('scandium',curprot)
      
      cmd.util.cnc(argname_P + ' and resi 4 ')
      selection_pepres = argname_P + ' and resi 4 '
      plusoneres = 'plusoneres'
      cmd.create('plusoneres',selection_pepres)
      cmd.show('sticks',plusoneres)
      #cmd.show('surface',plusoneres)
      #cmd.set('transparency',0.3,plusoneres)
      cmd.util.cnc(plusoneres)
      cmd.color('aquamarine',plusoneres+ ' and hydrogens')

      contactres = 'contact_residues'
      selection_string = '(%s) within %f of %s ' %(curprot,distance,selection_pepres)
      cmd.select('curselatoms',selection_string)
      cmd.select('contact_residues',' ( br. curselatoms )' )
      cmd.show('sticks','sc. contact_residues')
      #cmd.show('surface','contact_residues')
      #cmd.set('transparency',0.3,'contact_residues')
      cmd.color("hotpink",'contact_residues')
      cmd.util.cnc(contactres)
      cmd.color('hotpink',contactres+' and hydrogens')
      cmd.hide("sticks","hydrogens")

  cmd.bg_color(color="white")
  #cmd.set('sphere_scale',0.36)
  return curprot
  
curprot = openfromfile([sys.argv[1]])
cmd.hide('spheres','chain B')
cmd.set('two_sided_lighting')
cmd.set('ambient_occlusion_scale',11)
cmd.util.ray_shadows('occlusion2')
cmd.set('stick_radius',0.35)
cmd.set('stick_radius',0.20,'resname HWU')
if savefig:
  cmd.set('surface_quality',2)
cmd.set('spec_reflect',0.2)
#cmd.set_view([0.990826607,   -0.128174603,    0.042645268,
#     0.107263848,    0.938434541,    0.328366935,
#    -0.082108788,   -0.320781440,    0.943583190,
#     0.000177413,   -0.000515442,  -49.219161987,
#   -17.985958099,   36.096130371,  -10.925792694,
#    -0.336019993,   98.884994507,  -20.000000000])

cmd.set_view ([
    -0.754826546,   -0.081104644,   -0.650870502,
    -0.176844001,    0.980721891,    0.082886577,
     0.631624460,    0.177689806,   -0.754619777,
     0.000887860,   -0.001544237,  -60.999008179,
   -15.107797623,   33.454196930,  -15.351015091,
     5.642638683,  115.904205322,  -20.000000000 ])

cmd.scene('001','store')
cmd.set('ray_trace_mode',1)
cmd.set('ray_trace_gain',0.1)
cmd.set('ray_trace_disco_factor',1)
cmd.set('ray_shadows','off')
if savefig:
  #cmd.ray(800,600)
  cmd.png(sys.argv[-1]+ '_sticks')
