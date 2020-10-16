
def get_outfile_name(inpfile):
    name = inpfile.split('/')[0]
    return name

def get_keys(inpfile):
    name_base = get_outfile_name(inpfile)
    if name_base.find('T12')==-1:
      key1  = name_base[-10:-7]
      key2 = name_base[-3:]
      return key1,key2
    else:
      basename = name_base[:-13]
      print(name_base,basename)
      key1  = basename.split('PackandRelax_')[1]+'-'+name_base[-10:-7]
      key2 = name_base[-3:]
      print(name_base,key1,key2)
      return key1, key2

def get_key(inpfile):
    key1,key2 = get_keys(inpfile)
    return "%s,%s" %(key1,key2)

def get_keys_from_key(key):
    key1 = key.split(',')[0]
    key2 = key.split(',')[1]
    return key1,key2

def get_matrix_key(key):
    key1,key2 = get_keys_from_key(key)
    listres_3letter = ["LEU","VAL","ALA","ILE","MET","PRO","GLY","PHE","TRP","TYR","SER","THR","ASN","GLN","HIS","LYS","ARG","ASP","GLU"]
    listres = ['L','V','A','I','M','P','G','F','W','Y','S','T','N','Q','H','K','R','D','E']
    matrix_key = 'A%sT%sAPRC' %(listres[listres_3letter.index(key1)],listres[listres_3letter.index(key2)])
    return matrix_key

def get_matrix_key_from_file(inpfile):
  key1,key2 = get_keys(inpfile)
  key =  "%s,%s" %(key1,key2)
  return get_matrix_key(key)

def get_outfile_name_n(listoffiles):
    keys=[]
    for inpfile in listoffiles:
      key1, key2 = get_keys(inpfile)
      keys.append(key1+'-'+key2)
    return '_'.join(keys)

def get_full_name(curfile,type_='3d',ext='.png',info=None,suffix=''):
    import os
    if type_=='3d':
      basedir='results/scatter3d/'
      basename = 'scatter3d'
    if type_=='parigrid':
      basedir = 'results/pairplots/'
      basename = 'pairgrid'
    if type_ =='slice3d':
      basedir = 'results/scatter3dcompare/'
      basename = 'scatter3dslice'
    if type_ == 'cluster':
      basedir = 'results/clusters/'
      basename = 'clusters'
    if type_ == 'serialize_data':
      basedir = 'data/pickled/'
      basename = 'sortedscoredict'
      ext='.p'
    if type_ == 'json_data':
      basedir = 'data/json/'
      basename = 'sortedscoredict'
      ext='.json'
    if type_ == 'serialize_data_clusters':
      basedir = 'data/clusters/'
      basename = 'sortedscoredict_clusters'
      ext='.p'
    if type_== '3d_all':
      basedir = 'results/scatter3dcombine/'
      basename = 'scatter3dcombineAll'
      if not os.path.exists(basedir):
        os.system('mkdir -p %s' %basedir)
      if info is None:
        name = basedir + basename  + '_' + suffix  +ext
      else:
        name = basedir + basename + '_'+ info['xfield'] +'_' + suffix +ext
      return name
    if type_== '3d_X-1':
      basedir = 'results/scatter3dcombine/'
      basename = 'scatter3dcombineX-1'
      if not os.path.exists(basedir):
        os.system('mkdir -p %s' %basedir)
      key1,key2 = get_keys(curfile)
      if info is None:
        name = basedir + basename  + '_' + key1 + '_' + suffix +'_' +ext
      else:
        name = basedir + basename + '_'+ info['xfield'] + '_' + key1 +'_' + suffix +'_' +ext
      return name
    if type_== '3d_X+1':
      basedir = 'results/scatter3dcombine/'
      basename = 'scatter3dcombineX+1'
    if not os.path.exists(basedir):
      os.system('mkdir -p %s' %basedir)
    if info is None:
      name = basedir + basename + '_' + get_outfile_name(curfile)+ext
    else:
      name = basedir + basename + '_'+ info['xfield'] + '_' + get_outfile_name(curfile)+ext
    return name  

def _test_3d_view(ax):
    base_azim =ax.azim
    base_elev = ax.elev
    for angle in range(20,40,10):
      ax.view_init(base_elev  -10 , base_azim + angle)
      plt.draw()
      plt.pause(1)

def axes_presets(ax):

    from matplotlib import rc
    rc('font',size=28)
    rc('font',family='serif')
    rc('axes',labelsize=32)
    #ax.grid(False)
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis._axinfo['tick']['inward_factor'] = 0
    ax.xaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.yaxis._axinfo['tick']['inward_factor'] = 0
    ax.yaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.zaxis._axinfo['tick']['inward_factor'] = 0
    ax.zaxis._axinfo['tick']['outward_factor'] = 0.1
    ax.xaxis.labelpad = 20
    ax.yaxis.labelpad = 20
    ax.zaxis.labelpad = 15
    ax.tick_params(labelsize=16, grid_linewidth=1.0, grid_linestyle='-', grid_color='black',grid_alpha=0.6)
    [t.set_va('center') for t in ax.get_yticklabels()]
    [t.set_ha('left') for t in ax.get_yticklabels()]
    [t.set_va('center') for t in ax.get_xticklabels()]
    [t.set_ha('right') for t in ax.get_xticklabels()]
    [t.set_va('center') for t in ax.get_zticklabels()]
    [t.set_ha('left') for t in ax.get_zticklabels()]
    ax.xaxis.set_rotate_label(False)
    ax.yaxis.set_rotate_label(False)
    ax.zaxis.set_rotate_label(True)

def get_colorbar(ax):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="3%")
    return cax
