# usage: link_mat num_sample degree import sys,os import networkx as nx import numpy as np from operator import itemgetter from sklearn.preprocessing import Binarizer from sklearn.utils import shuffle

import sys,os
import networkx as nx
import numpy as np
from operator import itemgetter
from sklearn.preprocessing import Binarizer
from sklearn.utils import shuffle
from sklearn import cross_validation, svm
from sklearn.metrics import roc_curve, auc

def count_net(outfile, datafile): 
	"""count AUC for top identified subnet by classifying the expression data"""
	ifp = open(outfile)
	ifp.readline()
	tar = ifp.readline().rstrip().split()[2:]
	ifp.close()
	ifp = open(datafile)
	data_dict = {}
	for line in ifp:
		item = line.rstrip().split()
		data_dict[item[0]] = map(float,item[1:])
	ifp.close()
	mat = []
	for each in tar:
		mat.append(data_dict[each])
	mat = np.array(mat).transpose()
	X_train, X_test, y_train, y_test = cross_validation.train_test_split(mat, data_dict['outcome'], test_size=0.2, random_state=0)
	clf = svm.SVC(kernel='linear', probability=True)
	probas_ = clf.fit(X_train, y_train).predict_proba(X_test)
	fpr, tpr, thresholds = roc_curve(y_test, probas_[:, 1])
	roc_auc = auc(fpr, tpr)
	return roc_auc

itertime = 10 
num_gene = 500
num_sample = 100
np.random.seed(1234567890)
random_state = np.random.RandomState(0)
output_fp = open('linear.sim.txt','w')

for head in ['in','out']:
	print >> output_fp, 'head '+ head
	for pec in [1, 0.9, 0.8, 0.7, 0.6]:
		ego_count = []
		chuang_count = []
		for i in range(itertime):
			
			# generating the network
			G = nx.scale_free_graph(num_gene, seed=i)
			Gsim = nx.Graph(G)
			Gsim.remove_edges_from(Gsim.selfloop_edges())

			# generating the data matrix
			data_mat = np.random.randn(num_sample,num_gene)

			# generating the Y which has linear reliationship with selected genes
			node_and_degree=Gsim.degree()
			sort_node = sorted(node_and_degree.items(),key=itemgetter(1))
			deg = np.random.randint(5,20) # degree between [5,20)
			for each in sort_node:
				if each[1] >= deg:
					head_node = each[0]
					break
			try:
				hub_ego=nx.ego_graph(Gsim,head_node,radius = 1) # step 1 only
			except NameError:
				head_node = each[0]
				hub_ego=nx.ego_graph(Gsim,head_node,radius = 1)
			index = shuffle(hub_ego.nodes(),random_state=random_state)
			head_set = set([head_node])
			subidx = set(index[:int(pec*len(index))])
			if head == 'in':
				subidx |= head_set
			else:
				subidx -= head_set
			Y = np.zeros(num_sample)
			Y[::5] += 3 * (0.5-np.random.rand(num_sample/5)) # add noise to targets
			for each in subidx:
				Y += data_mat[:,each]
			binarizer = Binarizer()
			label = binarizer.transform(Y)

			# output the gene expression matrix
			ofp = open('linear.genemat','w')
			for each in sorted(Gsim.nodes()):
				print >> ofp, str(each)+'\t'+'\t'.join(map(str,data_mat[:,each]))
			print >> ofp, 'outcome\t'+'\t'.join(map(str,label))
			ofp.close()
	
			# output the network adjlist file
			nx.write_adjlist(Gsim,"linear.adjlist")
			os.system('epd_python egonet.py -n linear.adjlist -g linear.genemat -o linear.egonet.txt -s 0')
			ego_count.append(count_net('linear.egonet.txt','linear.genemat'))
			os.system('epd_python chuang.py -n linear.adjlist -g linear.genemat -o linear.chuang.txt')
			chuang_count.append(count_net('linear.chuang.txt','linear.genemat'))
		print >> output_fp, str(pec*100)+'%'+'\t%.2f (%.2f)\t%.2f (%.2f)' % (np.mean(ego_count), np.std(ego_count), np.mean(chuang_count), np.std(chuang_count))
