Introduction
------------
Subnet is implemented by Python and it is designed to detecting disease related subnetwork from a large biological network (PPI, metabolic network) combined with gene expression data.

Pre-installtalation
-------------------
Python version 2.7 or later (http://www.python.org/)

Python packages:

* Biopython (http://biopython.org/wiki/Main_Page)
* NetworkX (http://networkx.github.io/)
* SciPy (http://www.scipy.org/)
* Numpy (http://www.numpy.org/)
* scikit-learn (http://scikit-learn.org/stable/)

Running Subnet
------------------
#### Command-line usage:
    python egonet.py -n <network_file> -g <gene_matrix_file> -o <output_file> [opts]
#### Options:
	-m <int>	:method of classification or regression (default: class)
	-t <float>	:percentage of top selected gene for searching (default: 1, no sort)
	-s <float>	:score cutoff for printing selected subnetwork (default: 0.6)
	-f <pickle>	:saved subnetwork python object used for visualization (default: subnetwork.py)
	-r <txt>	:saved gene list ranked by two measuring methods (default: gene_rank.txt)
	-h      	:produce this menu
#### Example:
    python egonet.py -n sample_data/input/network.adjlist -g sample_data/input/gene_expression.txt -o TNBC.txt -f svm_net.pk
Example data are provide in the directory *sample_data/*

Input and Output File
------
1. network file is the adjacency list format
2. gene matrix file starts with gene name or entrize id as first column and expression values for other columns, the last row starts with "outcome" and labels for each sample
3. the output files contain the ranked subnetwork by predicting accuracy and ranked genes names using M-value.
 M = m*s*i
where m is total number of subnetwork contained gene, s is the score of each subnetwork and i is the importance of gene.

Visualize subnetwork
-------------------
Selected subnetwork can be plotted using as followed:

    python script/drawnet.py mark_gene diff_gene network_obj gene_matrix_file node

#### Example:
    python script/drawnet.py sample_data/visualization/breastcancer.gene sample_data/visualizaiont/diffexpress.gene svm_net.pk sample_data/input/gene_expression.txt 675

Contact us
----------
Questions, suggestions, comments, etc?

Author: Rendong Yang

Send email to cauyrd@gmail.com

Referencing Subnet 
----------------------
updating soon!
