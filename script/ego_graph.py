#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example using the NetworkX ego_graph() function to return the main egonet of 
the largest hub in a Barabási-Albert network.
"""
__author__="""Drew Conway (drew.conway@nyu.edu)"""

from operator import itemgetter
import networkx as nx
import matplotlib.pyplot as plt


if __name__ == '__main__':
    # Create a BA model graph
    n=100
    m=2
    G=nx.generators.barabasi_albert_graph(n,m,seed=0)
    # find node with largest degree
    node_and_degree=G.degree()
    (largest_hub,degree)=sorted(node_and_degree.items(),key=itemgetter(1))[-2]
    # Create ego graph of main hub
    hub_ego=nx.ego_graph(G,largest_hub)
    # Draw graph
    pos=nx.spring_layout(hub_ego)
    nx.draw(hub_ego,pos,node_color='r',node_size=100,with_labels=False)
    # Draw ego as large and red
    nx.draw_networkx_nodes(hub_ego,pos,nodelist=[largest_hub],node_size=600,node_color='w')
    #nx.draw_networkx_nodes(hub_ego,pos,nodelist=[73,75],node_size=100,node_color='w')
    labels = {largest_hub:'EGO'}
    nx.draw_networkx_labels(hub_ego,pos,labels,nodelist=[largest_hub],font_size=16)
    plt.savefig('ego_graph2.png')
    plt.show()
