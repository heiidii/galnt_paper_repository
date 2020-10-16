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
      #cmd.show_as('surface',curprot)
      cmd.set('transparency',surface_transparency,curprot)
      argname_HWU = "resname HWU"
      #cmd.show('sticks',argname_HWU)
      #cmd.color('orange',argname_HWU)
      argname_P_subset = argname_P + ' and ( resi 7 or resi 8)'
      cmd.hide('sticks',argname_P_subset)
      #cmd.show('cartoon',curprot)
      #cmd.color('scandium',curprot)
      
      cmd.util.cnc(argname_P + ' and resi 2 ')
      selection_pepres = argname_P + ' and resi 2 '
      plusoneres = 'plusoneres'
      cmd.create('plusoneres',selection_pepres)
      cmd.show('sticks',plusoneres)
      cmd.show('surface',plusoneres)
      cmd.set('transparency',0.3,plusoneres)
      cmd.util.cnc(plusoneres)
      cmd.color('aquamarine',plusoneres+ ' and hydrogens')

      contactres = 'contact_residues'
      selection_string = ' (%s) within %f of %s ' %(curprot,distance,selection_pepres)
      contact_residues =  'br. ( %s )' %selection_string
      #cmd.select('curselatoms',selection_string)
      #cmd.select('contact_residues','br. curselatoms' )
      cmd.show('sticks',contact_residues)
      cmd.show('surface',contact_residues)
      cmd.set('transparency',0.3,contact_residues)
      cmd.color("lightpink",contact_residues)
      cmd.util.cnc(contact_residues)
      cmd.color('lightpink',' (%s) and hydrogens' %contact_residues)
      cmd.hide("sticks","hydrogens")

  cmd.bg_color(color="white")
  #cmd.set('sphere_scale',0.36)
  cmd.select('%s' %contactres,contact_residues )
  myspace = {'reslist = []'}
  #cmd.iterate( '(%s)' %contactres, 'reslist.append((resi,resn))',space=myspace)
  return curprot
  
curprot = openfromfile([sys.argv[1]])
cmd.hide('spheres','chain B')
cmd.set('two_sided_lighting')
cmd.set('ambient_occlusion_scale',11)
cmd.util.ray_shadows('occlusion2')
cmd.set('stick_radius',0.35)
#cmd.set('stick_radius',0.20,'resname HWU')
cmd.hide('sticks','resname HWU')
if savefig:
  cmd.set('surface_quality',2)
cmd.set('spec_reflect',0.2)
cmd.set_view([0.504837930,    0.002740630,   -0.863198340,
     0.243603423,    0.958888113,    0.145517126,
     0.828116775,   -0.283741266,    0.483424038,
     0.000439802,   -0.000729280,  -51.568401337,
   -19.725727081,   38.694969177,  -21.041263580,
    -3.494707108,  106.578361511,  -20.000000000])
cmd.scene('001','store')
cmd.set('ray_trace_mode',1)
cmd.set('ray_trace_gain',0.1)
cmd.set('ray_trace_disco_factor',1)
cmd.set('ray_shadows','off')
if savefig:
  #cmd.ray(800,600)
  cmd.png(sys.argv[-1])
#print(reslist)
