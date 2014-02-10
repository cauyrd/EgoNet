# draw selected subnetwork
# usage: python drawgraph.py import_gene_list diff_gene_list networkobj genemat nodename 
import sys,pickle
from matplotlib import pyplot as plt
import networkx as nx
from sklearn.ensemble import ExtraTreesClassifier
import numpy as np
from Bio import Entrez
Entrez.email = "yourname@domain.com"     # Always tell NCBI who you are

def data_input(file_matrix):
	# read matrix file
	data_dict = {}
	ifp = open(file_matrix)
	for line in ifp:
		item = line.rstrip().split()
		data_dict[item[0]] = map(float,item[1:])
	label = data_dict['outcome']
	del data_dict['outcome']
	return data_dict,label

def select_mat(data_dict,sub_fea,G):
	mat = []
	for each in sub_fea:
		try:
			mat.append(data_dict[each])
		except KeyError:
			print 'node:'+str(each)+' not found in gene expression data'
			G.remove_node(each)
	return np.array(mat).transpose()

ifp = open(sys.argv[1])
mark_gene = set([line.rstrip() for line in ifp])
ifp.close()
ifp = open(sys.argv[2])
diff_gene = set([line.rstrip() for line in ifp])
ifp.close()
ifp = open(sys.argv[3])
net_dict = pickle.load(ifp)
G = net_dict[sys.argv[5]][1]
features = sorted(G.nodes())
data_dict, label = data_input(sys.argv[4])
X = select_mat(data_dict, features, G)
new_features = sorted(G.nodes())

#####################################################################
# computing feature importance using forests of trees
clf = ExtraTreesClassifier(random_state=0, compute_importances=True)
clf.fit(X,label)
importance = clf.feature_importances_
nodesize = 1e4*importance
myset = set(new_features)
tumorlist = list(myset.intersection(mark_gene))
unknownlist = list(myset.difference(mark_gene))
difflist = list(myset.intersection(diff_gene))
nodifflist = list(myset.difference(diff_gene))
node_dict = dict(zip(new_features,nodesize))

#####################################################################
# draw network
#pos=nx.spring_layout(G)
pos=nx.graphviz_layout(G)

# getting gene symbol from entrez id
idlist = ','.join(new_features)
handle = Entrez.esummary(db='gene',id=idlist)
record = Entrez.read(handle)
names = [n['Name'] for n in record]
labels = dict(zip(new_features, names))
plt.figure(figsize=(8,8))

# import gene highlight subfigure
plt.subplot(221)
nx.draw_networkx_nodes(G,pos,nodelist=tumorlist, node_size=[node_dict[v] for v in tumorlist], node_color='#FF7676')
nx.draw_networkx_nodes(G,pos,nodelist=unknownlist, node_size=[node_dict[v] for v in unknownlist], node_color='#A0CBE2')
nx.draw_networkx_edges(G,pos,alpha=0.5,width=4)
nx.draw_networkx_labels(G,pos,labels,font_size=13)
plt.title("score = %6.3f"%(net_dict[sys.argv[5]][0]))
plt.axis('off')

# different expressed gene highlight subfigure
plt.subplot(222)
nx.draw_networkx_nodes(G,pos,nodelist=difflist, node_size=[node_dict[v] for v in difflist], node_color='#00B100')
nx.draw_networkx_nodes(G,pos,nodelist=nodifflist, node_size=[node_dict[v] for v in nodifflist], node_color='#A0CBE2')
nx.draw_networkx_edges(G,pos,alpha=0.5,width=4)
nx.draw_networkx_labels(G,pos,labels,font_size=13)
plt.title("score = %6.3f"%(net_dict[sys.argv[5]][0]))
plt.axis('off')
plt.tight_layout()
plt.savefig(sys.argv[5]+".png") # save as png
plt.show()
