limits = { 'substrate_ca_rmsd':[0.4,1.0],'interaction_energy':[-36,-26],'distance_catalysis':[2.8,3.3], 'interaction_energy_binary':[0,1],'normalized_2_combined':[0,2] , 'binary':[0,1], "phi_500":[-100,-50,0],"psi_500":[130,170,0],"rmsddih_xp1":[0.015,0.09] ,'HWUO1B-THR7N':[3.2,4.2,0] ,'HWUO1-THR7N':[4.7,5.7,0]}

limitsbins = { 'interaction_energy':[-42,0],'distance_catalysis':[2.5,4.0],'substrate_ca_no_super_rmsd':[0.15,4.0], 'substrate_sequon_ca_no_super_rmsd':[0.15,4.0], 'substrate_ca_rmsd':[0.15,1.50], 'binary':[0.0,1.0], 'normalized_2_combined':[0.0,2.0], 'HWUO1B-THR7N':[3.2,4.2,0] ,'HWUO1-THR7N':[4.7,5.7,0]}


info['limits']=limits

colors = { 'interaction_energy':'RdPu_r','distance_catalysis':'Greens_r','substrate_ca_no_super_rmsd':'RdPu_r', 'substrate_sequon_ca_no_super_rmsd':'Blues_r', 'substrate_ca_rmsd': 'Oranges_r' , 'normalized_2_combined':'Purples', 'interaction_energy_binary':'Purples'}
bins = {}
for entry in limitsbins:
  dbin = (limitsbins[entry][1] - limitsbins[entry][0])/18.0
  bins[entry]=[ limitsbins[entry][0] + dbin*i for i in range(0,int((limitsbins[entry][1] - limitsbins[entry][0])/dbin)+1)]
  print bins[entry]

info['rocbins'] = bins
info['limitsbins'] = limitsbins
