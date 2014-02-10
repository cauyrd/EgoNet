# usage: link_mat num_sample degree 
import sys,os
import networkx as nx
import numpy as np
from operator import itemgetter
import random
from sklearn.preprocessing import Binarizer

def count_net(myfile, nodes):
	"""docstring for count_net"""
	ifp = open(myfile)
	ifp.readline()
	tar = map(int,ifp.readline().rstrip().split()[2:])
	ifp.close()
	if sorted(tar) == sorted(nodes):
		return 1
	else:
		return 0


itertime = 100
num_gene = 500
num_sample = 100
np.random.seed(1234567890)
random.seed(1234567890)

svm_count = 0
rf_count = 0
knn_count = 0
for i in range(itertime):
	G = nx.scale_free_graph(num_gene, seed=i)
	Gsim = nx.Graph(G)
	Gsim.remove_edges_from(Gsim.selfloop_edges())

	# generating the data matrix
	data_mat = np.random.randn(num_sample,num_gene)

	# generating the Y which has linear reliationship with selected genes
	node_and_degree=Gsim.degree()
	sort_node = sorted(node_and_degree.items(),key=itemgetter(1))
	deg = np.random.randint(5,20) # degree between [5,10)
	for each in sort_node:
		if each[1] >= deg:
			head_node = each[0]
			break
	try:
		hub_ego=nx.ego_graph(Gsim,head_node,radius = 1) # step 1 only
	except NameError:
		head_node = each[0]
		hub_ego=nx.ego_graph(Gsim,head_node,radius = 1)

	index = hub_ego.nodes()
	#pec = random.uniform(0.5,0.8) # percentage of nodes selected between [0.5,0.8]
	pec = 0.8
	random.shuffle(index)
	subidx = index[:int(pec*len(index))]
	Y = np.zeros(num_sample)
	Y[::5] += 3 * (0.5-np.random.rand(num_sample/5)) # add noise to targets
	for each in subidx:
		Y += data_mat[:,each]
	binarizer = Binarizer()
	label = binarizer.transform(Y)

	# output the gene expression matrix
	ofp = open('linear.'+str(i)+'.genemat','w')
	for each in sorted(Gsim.nodes()):
		print >> ofp, str(each)+'\t'+'\t'.join(map(str,data_mat[:,each]))
	print >> ofp, 'outcome\t'+'\t'.join(map(str,label))
	ofp.close()
	#print 'significant network',index
	nx.write_adjlist(Gsim,"linear."+str(i)+".adjlist")
	os.system('epd_python svmnet.py -n linear.'+str(i)+'.adjlist -g linear.'+str(i)+'.genemat -o linear.svm.'+str(i)+'.txt -s 0')
#	os.system('epd_python ../rfnet.py -n linear.adjlist -g linear.genemat -o linear.rf.txt -s 0 -r 20')
	os.system('epd_python knnnet.py -n linear.'+str(i)+'.adjlist -g linear.'+str(i)+'.genemat -o linear.knn.'+str(i)+'.txt -s 0')
	svm_count += count_net('linear.svm.'+str(i)+'.txt',index)
	#rf_count += count_net('linear.rf.txt',index)
	knn_count += count_net('linear.knn.'+str(i)+'.txt',index)
	print 'network:',i
	print 'svm:',svm_count
#	print 'rf:',rf_count
	print 'knn:',knn_count
	#For testing...
	print 'index',index
	print 'subidx',subidx
	print 'head node',head_node
	
