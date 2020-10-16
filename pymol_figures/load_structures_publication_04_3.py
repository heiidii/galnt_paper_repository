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
      #cmd.show_as('surface',curprot)
      cmd.set('transparency',surface_transparency,curprot)
      argname_HWU = "resname HWU"
      cmd.select(argname_HWU)
      cmd.color('orange',"sele")
      argname_P_subset = argname_P + ' and ( resi 7 or resi 8)'
      cmd.hide('sticks',argname_P_subset)
      argname_k288 = curprot + " and resi 288 "
      cmd.show("sticks",argname_k288 )
      cmd.color("lightpink",argname_k288)
      cmd.set('transparency',0.4,curprot +' and (resi 288) ')
      cmd.set('transparency',1.0,curprot +' and (resi 289) ')
      cmd.util.cnc(argname_P + ' and ( resi 2 )')
      cmd.color('lightpink',argname_k288+' and hydrogens')
      cmd.util.cnc(argname_k288)
      cmd.show('cartoon', curprot+ ' and (resi 286-292)')#284-294)')
      cmd.distance('dhbond',argname_P + ' and ( resi 2 ) ',argname_k288 + '' ,mode=2)
      cmd.hide('labels','dhbond')
      cmd.set('dash_gap',0.5)
      cmd.set('dash_radius',0.18)
      cmd.show("sticks","hydrogens and " +  argname_P + ' and ( resi 2 )')
      
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
cmd.set('stick_radius',0.25, 'resname HWU')
cmd.set('stick_radius',0.25, 'chain P')
if savefig:
  cmd.set('surface_quality',2)
cmd.set('spec_reflect',0.15)
cmd.set_view([0.559028506,    0.038134605,   -0.828269780,
     0.229416862,    0.952828348,    0.198710650,
     0.796770990,   -0.301107675,    0.523911476,
     0.000145853,    0.000030812,  -53.878547668,
   -19.991203308,   38.031669617,  -19.044155121,
  -7436.402832031, 7544.086914062,  -20.000000000 ])

cmd.set_view([ 0.560823977,   -0.057994414,   -0.825901330,
     0.280622661,    0.951809168,    0.123718791,
     0.778919101,   -0.301154196,    0.550075233,
     0.000162907,   -0.000101577,  -37.700313568,
   -20.097080231,   37.180461884,  -19.248620987,
  -7452.602539062, 7527.887207031,  -20.000000000])
cmd.set_view([0.560823977,   -0.057994414,   -0.825901330,
     0.280622661,    0.951809168,    0.123718791,
     0.778919101,   -0.301154196,    0.550075233,
     0.000162907,   -0.000101577,  -50.542243958,
   -20.097080231,   37.180461884,  -19.248620987,
  -7439.758789062, 7540.730957031,  -20.000000000])
cmd.set_view ([ 0.562684298,    0.041209981,   -0.825643718,
     0.146542013,    0.977968335,    0.148682088,
     0.813580871,   -0.204654470,    0.544249058,
    -0.000095510,    0.000182927,  -60.292354584,
   -11.669906616,   36.101730347,  -27.599609375,
     5.149775505,  115.411437988,  -20.000000000] )
cmd.set('ray_trace_mode',1)
cmd.set('ray_trace_gain',0.1)
cmd.set('ray_trace_disco_factor',1)
cmd.set('ray_shadows','off')
#cmd.log_open('log_04_2.pml')
#cmd.get_view()
if savefig:
  cmd.ray(800,600)
  cmd.png(sys.argv[-1])
