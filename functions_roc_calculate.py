import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle

from sklearn import svm, datasets
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp

from matplotlib import rc

rc={'xtick.labelsize': 20, 'ytick.labelsize': 20,"font.weight":'bold', 'xtick.color':'black','ytick.color':'black','xtick.top':False, 'xtick.labeltop':False, 'xtick.bottom':True, 'xtick.labelbottom':True}

def plot_roc(fpr,tpr,filename='testroc.png',title='',str_score='',show=False):
  plt.rcParams.update(**rc)
  fig,ax = plt.subplots(figsize=(12,12))
  ax.spines['bottom'].set_color('black')
  ax.spines['top'].set_color('black')
  ax.spines['right'].set_color('black')
  ax.spines['left'].set_color('black')
  fig.suptitle(title,fontsize=25)
  plt.plot([0.0,1.0],[0.0,1.0],color='black',linestyle='dashed',linewidth=3)
  plt.plot(fpr,tpr,color='mediumvioletred',linewidth=7)
  plt.text(0.5,0.1,str_score,size=42,color='black')
  plt.ylim(0.0,1.0)
  plt.xlim(0.0,1.0)
  axis_font = { 'size':48}
  plt.xlabel("False Positive Rate (FPR)",**axis_font)
  plt.ylabel("True Positive Rate (TPR)",**axis_font)
  plt.tick_params(axis='both',labelsize=35)
  plt.tight_layout()
  plt.subplots_adjust(top=0.85)
  plt.savefig(filename,transparent=True,dpi=800)
  if show:
    plt.show()

def plot_rocs(fprs,tprs,filename='testroc.png',title='',str_scores=[],colors=[],show=False):
  plt.rcParams.update(**rc)
  fig,ax = plt.subplots(figsize=(12,12))
  ax.spines['bottom'].set_color('black')
  ax.spines['top'].set_color('black')
  ax.spines['right'].set_color('black')
  ax.spines['left'].set_color('black')
  fig.suptitle(title,fontsize=25)
  plt.plot([0.0,1.0],[0.0,1.0],color='black',linestyle='dashed',linewidth=3)
  for i,fpr in enumerate(fprs):
    if len(colors) == len(fprs):
      plt.plot(fprs[i],tprs[i],color=colors[i],linewidth=7)
    else:
      plt.plot(fprs[i],tprs[i],color='black',linewidth=7)
    if len(str_scores) == len(fprs):
      plt.text(0.5,0.05 + 0.09*i,str_scores[i],size=35,color=colors[i])
  plt.ylim(0.0,1.0)
  plt.xlim(0.0,1.0)
  axis_font = { 'size':48}
  plt.xlabel("False Positive Rate (FPR)",**axis_font)
  plt.ylabel("True Positive Rate (TPR)",**axis_font)
  plt.tick_params(axis='both',labelsize=35)
  plt.tight_layout()
  plt.subplots_adjust(top=0.85)
  plt.savefig(filename,transparent=True,dpi=800)
  if show:
    plt.show()
  plt.close()

def twovariable_roc_curve(true_labels,scores):
  xmax = 0.2
  xmin = 2.0
  ymax = 1.0
  ymin = 0.0
  delta_x = 0.01 #rmsd less is better
  Nx = (xmax - xmin)/delta_x + 1
  xarray =  [ xmin + delta_x*t for t in range(0,Nx)]
  delta_y = 0.01 #fraction more is better
  Ny = (ymax - ymin)/delta_y + 1
  yarray =  [ ymin + delta_y*t for t in range(0,Ny)]
  fprs=[]
  tprs=[]
  for x in xarray:
    for y in yarray:
      for i in range(scores.shape[0]):
        ax = scores[i,0]
        by = scores[i,1]
        if true_labels[i]>0:
          if ax > x and bx > y :
            tpr +=1
          else:
            fpr+=1
        else:
          if ax > x and bx > y :
            fpr +=1
          else:
            tpr+=1
      fpr = fpr/float(scores.shape[0])
      tpr = tpr/float(scores.shape[0])
      fprs.append(fpr)
      tprs.append(tpr)
  plot_roc(fprs,tprs)
  

def simplistic_roc_curve(true_labels,scores,pos_labels=1): #pos_labels 0 and 1
  '''
  Parameters:
    parameters are exactly the same as parameters for roc curve
  Returns:
    fpr , tpr, thresholds
  '''
  fpr , tpr, thresholds = roc_curve(true_labels,scores,pos_labels)  
  auc_score = roc_auc_score(true_labels,scores)
  return fpr,tpr,thresholds,auc_score

def feature_selection(true_labels=None,scores_matrix=None):
  print(__doc__)

  from sklearn.svm import SVC
  from sklearn.datasets import load_digits
  from sklearn.feature_selection import RFE
  import matplotlib.pyplot as plt

  # Load the digits dataset
  digits = load_digits()
  X = digits.images.reshape((len(digits.images), -1))
  y = digits.target
  print(y.shape)
  print('X',X.shape)
  #exit()
  Y = true_labels.ravel()
  X = scores_matrix.reshape(Y.shape[0],-1)
  print(Y.shape)
  print('X',X.shape)
  #exit()
  # Create the RFE object and rank each pixel
  svc = SVC(kernel="linear", C=1)
  rfe = RFE(estimator=svc, n_features_to_select=1, step=1)
  rfe.fit(X, Y)
  ranking = rfe.ranking_
  print(ranking)

  # Plot pixel ranking
  #plt.matshow(ranking, cmap=plt.cm.Blues)
  #plt.colorbar()
  #plt.title("Ranking with RFE")
  #plt.show()

def feature_selection_trees(true_labels=None, scores_matrix=None):
  
  # Feature Importance
  from sklearn import metrics
  from sklearn.ensemble import ExtraTreesClassifier
  # load the iris datasets
  # fit an Extra Trees model to the data
  model = DecisionTreeClassifier()
  Y = true_labels.ravel()
  X = scores_matrix.reshape(Y.shape[0],-1)
  print(Y.shape)
  print('X',X.shape)
  model.fit(X,Y)
  # display the relative importance of each attribute
  print('trees: importances',model.feature_importances_)
  #print('trees: estimators ',model.estimators_)
  #print('trees: classes ',model.classes_)
  Z = model.score(X, Y)
  print('scores: ',Z)
  #Z = Z.reshape(xx.shape)
  #cs = plt.contourf(xx, yy, Z, cmap=cmap)

def feature_selection_rfe(true_labels=None, scores_matrix=None):
  
  # Recursive Feature Elimination
  from sklearn import datasets
  from sklearn.feature_selection import RFE
  from sklearn.linear_model import LogisticRegression
  # load the iris datasets
  # create a base classifier used to evaluate a subset of attributes
  model = LogisticRegression()
  # create the RFE model and select 3 attributes
  Y = true_labels.ravel()
  X = scores_matrix.reshape(Y.shape[0],-1)
  print(Y.shape)
  print('X',X.shape)
  rfe = RFE(model, 3)
  rfe = rfe.fit(X,Y)
  # summarize the selection of the attributes
  print('support: ',rfe.support_)
  print('ranking: ',rfe.ranking_)


def decision_tree_regression(true_labels=None,scores_matrix=None,list_names=None,max_depth=3):
  from sklearn import tree
  from sklearn.tree.export import export_text
  import graphviz
  clf = tree.DecisionTreeClassifier(random_state=0, max_depth=max_depth)
  Y = true_labels.ravel()
  X = scores_matrix.reshape(Y.shape[0],-1)
  clf = clf.fit(X,Y)

  r = export_text(clf, feature_names=list_names)
  f = open('decisiontree_'+'_'.join(list_names)+'.out','w')
  f.write('max_depth %d:\n'%max_depth)
  f.write(str(r))
  f.close()
  tree.plot_tree(clf.fit(X,Y),filled=True)
  plt.show()
  f1 = open('dotfile.dot','w')
  dot_data = tree.export_graphviz(clf,out_file=f1,
                                  feature_names=list_names,
                                  class_names = ['low','medium','high'],
                                  filled=True,
                                  rounded=True,
                                  special_characters=True)
  import os  
  os.unlink('dotfile.dot')  

  graph = graphviz.Source(dot_data)
  graph

if __name__=="__main__":
  # Import some data to play with
  iris = datasets.load_iris()
  X = iris.data
  y = iris.target

  # Binarize the output
  print("shape before binarize: ",y.shape)
  y = label_binarize(y, classes=[0, 1, 2])
  n_classes = y.shape[1]
  print("shape of y: ",y.shape)
  print(n_classes)
  print('feature_selection')
  feature_selection()
  
