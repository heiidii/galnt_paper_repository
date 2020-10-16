delta_ticks = {'substrate_ca_no_super_rmsd':1.0}

fields = {"interaction_energy":[-45,-5,0], 
          "distance_catalysis":[2.5,6.5,0],"substrate_ca_rmsd":[0,6.5,0], "substrate_ca_no_super_rmsd":[0,6.5,1],"substrate_sequon_ca_no_super_rmsd":[0,10.0,1],'distance_catalysis_HWUO1-THR7N':[4.0,8.0,0],'distance_catalysis_HWUO1B-THR7N':[2.5,6.5,0] , 'distance_catalysis_UDP2OPB-THR7N':[1.5,6.5,0] , 'distance_catalysis_UDP3OPB-THR7N': [3.0,8.0,0]}

default_fields=fields

labels = {'interaction_energy':'IE (REU)','distance_catalysis': 'd_MC', 'distance_catalysis_HWUO1-THR7N':'d_PO-Thr', 'distance_catalysis_HWUO1B-THR7N':'d_PO2-Thr',"substrate_ca_no_super_rmsd":'rmsd', "substrate_ca_no_super_rmsd":'rmsd_sequon'}

labels_df = {'substrate_ca_no_super_rmsd':'rmsd','interaction_energy':'IE','distance_catalysis': 'd_MC', 'distance_catalysis_HWUO1-THR7N':'d_UDPO-Thr', 'distance_catalysis_HWUO1B-THR7N':'d_UDPO2-Thr_HB','distance_catalysis_UDP2OPB-THR7N':'d_UDPO2-Thr_HB','distance_catalysis_UDP3OPB-THR7N':'d_UDPO2-Thr','substrate_sequon_ca_no_super_rmsd':'rmsd_sequon','description':'description'}

info ={}
info['zfield'] = 'distance_catalysis'
info['yfield'] = 'distance_catalysis_HWUO1B-THR7N'
info['z2field'] = 'distance_catalysis_HWUO1-THR7N'
info['xfield'] = 'substrate_ca_no_super_rmsd'
info['cfield'] = 'interaction_energy'
info['z3field'] = 'substrate_sequon_ca_no_super_rmsd'
info['description'] = 'description'

categories ={}
categories['non-polar']=['LEU','VAL','ALA','ILE','MET','PRO','GLY']
#categories['gly']=['GLY']
categories['aromatics']=['PHE','TYR','TRP']
#categories['pro']=['PRO']
categories['polar']=['SER','THR','ASN','GLN']
categories['acidic']=['ASP','GLU']
categories['basic']=['ARG','LYS','HIS']

aalist = ['L','V','A','I','M','P','G','F','W','Y','S','T','N','Q','H','K','R','D','E']