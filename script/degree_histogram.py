#!/usr/bin/env python
"""
Random graph from given degree sequence.
Draw degree histogram with matplotlib.
Usage:program.py network.adjlist gene_rank.txt
"""
__author__ = """Rendong Yang (cauyrd@gmail.com)"""

try:
    import matplotlib.pyplot as plt
    import matplotlib
except:
    raise

import networkx as nx
import sys

bn = nx.read_adjlist(sys.argv[1])
bn.remove_edges_from(bn.selfloop_edges())

bn_degree_seq=sorted(nx.degree(bn).values(),reverse=True) # degree sequence
#print "Degree sequence", degree_sequence
ifp = open(sys.argv[2])
ifp.readline()
scale = 100
sn_degree_seq=[]
for i,line in enumerate(ifp):
	if i > scale:
		break
	id = line.rstrip().split()[3]
	sn_degree_seq.append(bn.degree(id))
ifp.close()


plt.loglog(bn_degree_seq,'b-',marker='o', label='all nodes in PPI')
plt.loglog(sn_degree_seq,'r-',marker='o', label='top ranked nodes by M-value')
plt.legend()
plt.title("Degree rank plot")
plt.ylabel("degree")
plt.xlabel("rank")

plt.savefig("degree_histogram.pdf")
plt.show()

