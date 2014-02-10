# Input: network_file gene_matrix_file
# Python version 2.7 or later

from __future__ import division
import numpy as np
from collections import Counter
from scipy.stats import describe
import networkx as nx
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn import cross_validation
from sklearn import svm
from sklearn.utils import resample
from operator import itemgetter
import sys, getopt
import time,pickle
from Bio import Entrez
Entrez.email = "yourname@domain.com"     # Always tell NCBI who you are

del_nodeset = set()

def usage():
	"""showing help information"""
	print 'Usage:'
	print '	python svmnet.py -n <network_file> -g <gene_matrix_file> -o <output_file> [opts]'
	print 'Opts:'
	print ' -m <int>	:method of classification or regression (default: class)'
	print ' -t <float>	:percentage of top selected gene for searching (default: 1, no sort)'
	print ' -s <float>	:score cutoff for printing selected subnetwork (default: 0.6)'
	print ' -f <pickle>	:saved subnetwork python object used for visualization (default: subnetwork.py)'
	print ' -r <txt>	:saved gene list ranked by two measuring methods (default: gene_rank.txt)'
	print ' -h      	:produce this menu'

def data_input(file_graph, file_matrix):
	# read graph file
	G = nx.read_adjlist(file_graph)
	G.remove_edges_from(G.selfloop_edges())
	# read matrix file
	data_dict = {}
	ifp = open(file_matrix)
	for line in ifp:
		item = line.rstrip().split()
		if item[0] in G.nodes():
			data_dict[item[0]] = map(float,item[1:])
		elif item[0] == 'outcome':
			label = map(float,item[1:])
		else:
			continue
	return G,data_dict,label

def select_mat(data_dict,sub_fea):
	mat = []
	for each in sub_fea:
		try:
			mat.append(data_dict[each])
		except KeyError:
			del_nodeset.add(each)
	return np.array(mat).transpose()

def sorted_map(map):
	ms = sorted(map.iteritems(),key=lambda (k,v):(-v,k))
	return ms

def rename(tuple_list):
	idlist = ','.join([v[0] for v in tuple_list])
	handle = Entrez.esummary(db="gene", id=idlist)
	record = Entrez.read(handle)
	return [(n['Name'],v[1]) for n,v in zip(record,tuple_list)],idlist.split(',')

def get_net_significance(netdict_value, data_mat, label, method):
	random_state = np.random.RandomState(0)
	num_trials = 100
	null_dist = []
	sub_fea = tuple(sorted(netdict_value[1].nodes()))
	sub_mat = select_mat(data_mat,sub_fea)
	for i in range(num_trials):
		rand_label = resample(label,random_state=random_state)
		if method == 'class':
			clf = svm.SVC(cache_size=1000)
			clf.fit(sub_mat,rand_label)
			scores = cross_validation.cross_val_score(clf,sub_mat,rand_label,cv=5)
			null_dist.append(scores.mean())
		else:
			clf = svm.SVR(cache_size=1000)
			clf.fit(sub_mat,rand_label)
			null_dist.append(clf.score(sub_mat,rand_label))
	null_dist.append(netdict_value[0])
	index = np.argsort(null_dist)
	pvalue = 1 - index[-1]/float(num_trials)
	return pvalue

class Subnetwork():
	'''
	Random Forest Learner

	Attributes
	----------
	rf : integer, optional
		Nr of trees to learn (default: 101)
	frac : float, optional
		fample fraction
	'''
	def __init__(self, G=nx.Graph(), method='class'):
		self.G = G
		self.netdict = {}
		self.method = method
		self.nodes = G.nodes()
		self.depth = []
	
	def sort_gene(self,data_mat,label, top=0.1):
		'''
		computing the gene importance using randomforest and select top percentage genes for searching in divide_net
		'''
		feature = sorted(data_mat.keys())
		X = select_mat(data_mat, feature)
		if self.method == 'class':
			forest = RandomForestClassifier(n_estimators=len(feature)*2, compute_importances=True, random_state=0, n_jobs=-1)
		else:
			forest = RandomForestRegressor(n_estimators=len(feature)*2, compute_importances=True, random_state=0, n_jobs=-1)
		forest.fit(X,label)
		importance = forest.feature_importances_
		sort_index = np.argsort(importance)[::-1][:int(top*len(feature))]
		self.nodes = [feature[i] for i in sort_index]
		
	def divide_net(self,data_mat,label):
		if self.method == 'class':
			clf = svm.SVC(cache_size=1000)
		else:
			clf = svm.SVR(cache_size=1000)
		subnet_dict = {}
		#print 'searching '+str(len(self.nodes))+' nodes in all'
		for vex_num, node in enumerate(self.nodes):
			#print str(vex_num)+':node '+str(node)
			if node not in self.G or self.G.degree(node) < 2:
				continue
			score = -np.inf
			depth = 1 
			while True:
				hub_ego = nx.ego_graph(self.G,node,radius=depth)
				sub_fea = tuple(sorted(hub_ego.nodes()))
				if sub_fea in subnet_dict:
					mean_score = subnet_dict[sub_fea]
				else:
					sub_mat = select_mat(data_mat,sub_fea)
					if sub_mat.shape[0] == 0:
						break
					newclf = clf.fit(sub_mat,label)
					if self.method == 'class':
						scores = cross_validation.cross_val_score(newclf,sub_mat,label,cv=5)
						mean_score = scores.mean()
					else:
						mean_score = clf.score(sub_mat,label)
				if mean_score > score:
					score = mean_score
					depth += 1
					self.netdict[node] = [score,hub_ego]
				else:
					nodes = tuple(sorted(self.netdict[node][1].nodes()))
					if nodes in subnet_dict:
						try:
							del self.netdict[node]
						except ValueError:
							#print 'del netdict error for node:',node
							pass
					else:
						subnet_dict[nodes] = self.netdict[node][0]
						self.depth.append(depth-1)
					break
		return

	def print_net(self, output, data_mat, label, cutoff = 0.6, filename = 'subnetwork.pk'):
		"""sorting and print subnetwork based on the score in decreased order"""
		all_net = [(key, value[0]) for key, value in self.netdict.items()]
		subnet_index = [item[0] for item in sorted(all_net,key=itemgetter(1),reverse=True)]
		ofp = open(output,'w')
		print >> ofp, 'score,pvalue\thead_node\tall_nodes'
		for each in subnet_index:
			score = self.netdict[each][0]
			if score > cutoff:
				pvalue = get_net_significance(self.netdict[each], data_mat, label, self.method)
				nodes = self.netdict[each][1].nodes()
				print >> ofp, str(score)+','+str(pvalue)+'\t'+each+'\t'+'\t'.join(nodes)
		
		ofp = open(filename,'w')
		pickle.dump(self.netdict,ofp)
		ofp.close()
		return

	def rank_gene(self, data_mat, cutoff, output):
		"""nodes rank based on two metric"""
		count_dict = {}
		count2_dict = {}
		for each in self.netdict:
			if self.netdict[each][0] < cutoff:
				continue
			X = select_mat(data_mat,sorted(self.netdict[each][1].nodes()))
			if self.method == 'class':
				forest = RandomForestClassifier(compute_importances=True, random_state=0, n_jobs=-1)
			else:
				forest = RandomForestRegressor(compute_importances=True, random_state=0, n_jobs=-1)
			forest.fit(X,label)
			importance = forest.feature_importances_
			for node,imscore in zip(sorted(self.netdict[each][1].nodes()),importance):
				try:
					count_dict[node] += 1
					count2_dict[node] += self.netdict[each][0] * imscore
				except KeyError:
					count_dict[node] = 1
					count2_dict[node] = self.netdict[each][0] * imscore
		
		metric1,id1 = rename(sorted_map(count_dict)[:1000]) # network count
		metric2,id2 = rename(sorted_map(count2_dict)[:1000]) # weighted network count * importance
		ofp = open(output,'w')
		print >> ofp, 'EntrezID\tGeneName\tNetCount\tEntrezID\tGeneName\tM-value'
		for a,b,aid,bid in zip(metric1,metric2,id1,id2):
			print >> ofp, aid+'\t'+str(a[0])+'\t'+str(a[1])+'\t'+bid+'\t'+str(b[0])+'\t'+str(b[1])
		ofp.close()
		return

	def summary(self):
		"""summary basic statistic for identified subnetwork"""
		print str(len(self.netdict))+' subnetwork generated:'
		n, (smin, smax), sm, sv, ss, sk = describe([self.netdict[key][0] for key in self.netdict])
		print 'Subnet socre ['+str(smin)+', '+str(smax)+'] of mean '+str(sm)+' and var '+str(sv)
		n, (smin, smax), sm, sv, ss, sk = describe([len(self.netdict[key][1].nodes()) for key in self.netdict])
		print 'Subnet nodes size ['+str(smin)+', '+str(smax)+'] of mean '+str(int(sm))+' and var '+str(int(sv))
		counter = Counter(self.depth)
		print 'Subnet depth summary:'
		for each in sorted(counter.keys()):
			print 'depth '+str(each)+': '+str(counter[each])

if __name__=='__main__':

	# parameters parsing
	net_file = None
	gene_file = None
	output = None
	objfile = 'subnetwork.pk'
	mymethod = 'class'
	percent = 1 
	score = 0.6
	gene_rank_output = 'gene_rank.txt'
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'n:g:o:m:t:s:f:r:h')
	except getopt.GetoptError as err:
		print str(err)
		usage()
		sys.exit(2)
	for o,a in opts:
		if o == '-n':
			net_file = a
		elif o == '-g':
			gene_file = a
		elif o == '-o':
			output = a
		elif o == '-m':
			mymethod = a
		elif o == '-t':
			percent = float(a)
		elif o == '-s':
			score = float(a)
		elif o == '-f':
			objfile = a
		elif o == '-r':
			gene_rank_output = a
		elif o == '-h':
			usage()
			sys.exit()
		else:
			assert False, "unhandled option"
	if not net_file or not gene_file or not output:
		usage()
		sys.exit(2)

	G, data_dict, label = data_input(net_file, gene_file)
	start = time.clock()
	myclf = Subnetwork(G, mymethod)
	if percent < 1:
		print 'Using randomforest computing gene importance...'
		myclf.sort_gene(data_dict, label,percent)
		end1 = time.clock()
		print 'finish computing feature importance using '+str((end1-start)/60.0)+' mins.'
	print 'searching sub-network...'
	myclf.divide_net(data_dict, label)
	myclf.print_net(output,data_dict,label,score,objfile)
	myclf.rank_gene(data_dict, score, gene_rank_output)
	#print 'In all '+str(len(del_nodeset))+' nodes from network not inclued in gene expression data.'
	#print 'deleted node:',del_nodeset
	myclf.summary()
	end = time.clock()
	print 'program finished using '+str((end-start)/60.0)+' mins.'
