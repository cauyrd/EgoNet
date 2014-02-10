# This program is the implementation of subnetwork selection method described in t    he paper:
# Chuang, Han-Yu, et al. "Network-based classification of breast cancer metastasis    ." Molecular systems biology 3.1 (2007).

import numpy as np
from sklearn import metrics
import networkx as nx
from operator import itemgetter
import sys,getopt
import time
from scipy import stats
from scipy.stats import describe

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
			data_dict[item[0]] = stats.zscore(np.array(map(float,item[1:])))
		elif item[0] == 'outcome':
			label = map(float,item[1:])
		else:
			continue
	return G,data_dict,label

def get_score(node_set, data_dict, label):
	'''
	Calculate Mutual information score
	'''
	mat = []
	for each in node_set:
		mat.append(data_dict[each])
	mat = np.array(mat)
	mean_col = mat.mean(axis=0)
	count = int(np.log2(len(mean_col))+1)
	bins = np.histogram(mean_col,count)[1]
	index = np.digitize(mean_col,bins)
	values = [bins[i-1] for i in index]
	score = metrics.mutual_info_score(label,values)
	return score

class Subnetwork():
	'''
	Identify the subnetwork based on the method from Chuang et al.
	'''
	def __init__(self, G=nx.Graph(), depth=2, rate=0.05):
		self.G = G
		self.subnet_index = []
		self.netdict = {}
		self.nodes = G.nodes()
		self.depth = depth
		self.rate = rate
	
	def divide_net(self,data_mat,label):
		for seed in self.nodes:
			if self.G.degree(seed) < 2:
				continue
			try:
				data_mat[seed]
			except KeyError:
				continue
			score = get_score([seed],data_mat,label)
			innet_set = set([seed])
			outnet_set = set(self.G.neighbors(seed))
			# Extend the seed node
			while True:
				ext_node = None
				score0 = 0
				for node in outnet_set:
					try:
						data_mat[node]
					except KeyError:
						continue
					newset = list(innet_set)+[node]
					newscore = get_score(newset,data_mat,label)
					if newscore > score0:
						score0 = newscore
						ext_node = node
				if score0/score - 1 <= self.rate or nx.shortest_path_length(self.G,ext_node,seed) > self.depth:
					break
				innet_set.add(ext_node)
				outnet_set |= set(self.G.neighbors(ext_node))
				outnet_set -= innet_set
				score = score0
			self.netdict[seed] = [score,innet_set]
		return

	
	def print_net(self, output, significance = 0):
		"""sorting and print subnetwork based on the score in decreased order"""
		all_net = [(key, value[0]) for key, value in self.netdict.items()]
		self.subnet_index = [item[0] for item in sorted(all_net,key=itemgetter(1),reverse=True)]
		ofp = open(output,'w')
		print >> ofp, 'score\tsead_node\tall_nodes'
		for each in self.subnet_index:
			score = self.netdict[each][0]
			nodes = list(self.netdict[each][1])
			if score > significance:
				print >> ofp, str(score)+'\t'+each+'\t'+'\t'.join(nodes)
		ofp.close()
		return

	def summary(self):
		"""summary basic statistic for identified subnetwork"""
		print str(len(self.netdict))+' subnetwork generated:'
		n, (smin, smax), sm, sv, ss, sk = describe([self.netdict[key][0] for key in self.netdict])
		print 'Subnet MI socre range ['+str(smin)+', '+str(smax)+'] of mean '+str(sm)+' and var '+str(sv)
		n, (smin, smax), sm, sv, ss, sk = describe([len(self.netdict[key][1]) for key in self.netdict])
		print 'Subnet nodes size ['+str(smin)+', '+str(smax)+'] of mean '+str(int(sm))+' and var '+str(int(sv))

if __name__=='__main__':

	# parameters parsing
	net_file = None
	gene_file = None
	output = None
	score = 0
	depth = 2
	rate = 0.05
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'n:g:d:r:s:o:h')
	except getopt.GetoptError as err:
		print str(err)
#		usage()
		sys.exit(2)
	for o,a in opts:
		if o == '-n':
			net_file = a
		elif o == '-g':
			gene_file = a
		elif o == '-d':
			depth = int(a)
		elif o == '-r':
			rate = float(a)
		elif o == '-o':
			output = a
		elif o == '-s':
			score = float(a)
		elif o == '-h':
#			usage()
			sys.exit()
		else:
			assert False, "unhandled option"
	if not net_file or not gene_file or not output:
#		usage()
		sys.exit(2)

	start = time.clock()
	print 'program start.'
	G, data_dict, label = data_input(net_file, gene_file)
	print 'data loading done!'
	myclf = Subnetwork(G,depth,rate)
	print 'searching subnetwork...'
	myclf.divide_net(data_dict, label)
	print 'netwrok searching done!'
	myclf.print_net(output,score)
	myclf.summary()
	end = time.clock()
	print 'program finished using '+str((end-start)/60.0)+' mins.'
