 ## For application of MERLIN for Gene Regulator Network (GRN) Inference and downstream visualization refer to this github page: https://github.com/Roy-lab/merlin-p
 ## For visualization of GRNs follow the steps:
Steps of generating module specific network files

# Step1: Generatin adjacency matrices for module specific subnetwork for each sample
script:network_visualization_scripts/convertNet.py
Example run: convertNet.py expression_matrices/WT_Iep1.txt module3339_subnetwork.txt matlab_inputs/adj_WT_Ire1.txt matlab_inputs/WT_Ire1.txt

# Step2: Calculate regression and correlation coeefficient using outputs in step1 with following matlab scripts:
network_visualization_scripts/runls_beta.m  #for regression coefficients
network_visualization_scripts/runcc_beta.m  #for correlation coefficients

# Step3: Generating input files for Cytoscape:
script: network_visualization_scripts/makeCyto2.py

This script takes 7 imput files:

reg_files.txt: a file containing path to all files contains regression coefficient produced in previous step by runls.m script (this should also include path to the module specific subnetwork).
cc_files.txt: a file containing path to all files containing correlation coefficient produced in previous step by runcc.m script (this should also include path to the module specific subnetwork)
zeromean_expression_files.txt: a file containing path to all genotype specific files which contain mean expression of zeromean normalized matrix per genotype.
deg.txt: two column file showing all the genes in module specific subnetwork as first column and number of outgoing edges they have in subnetwork as second columns
module.txt: two column file showing all the genes in module specific subnetwork as first column and module number in second column
condition_names.txt:a file containing name of the cell types/conditions
in: output files prefix
Example run:

python makeCyto2.py matlab_outputs/reg_files.txt matlab_outputs/cc_files.txt expression_matrices/zeromean_expression_files.txt matlab_outputs/deg.txt matlab_outputs/module.txt condition_names.txt cytoscape_inputs/in; done

Example output per sample:
in_WT_Ire1_cc.txt: edges and correlation coefficient
in_WT_Ire1_reg.txt: edges and correlation coefficient
in_WT_Ire1_att.txt: attribute file which contains columns for gene, its status (eg.TF), weather it is Module gene or not, number of edges it has and mean of zeromean expression

# Step4: 
Note: modifiny att.txt files to change the TF (circle) in column2 to TF2 (diamond) if gene is not an enriched regulators in module, TAR (rectangle) if it is a module gene and will notchange the TF if the gene is an enriched regulators in module using update_attribute_files.R script

# Step5: loading files to Cytoscape.
Cytoscape version 3.10.2
First you have to load the cytoscape_inputs/network_vizualisation.xml file to the cytoscape which contain the certain beautificaitons like edge colors red, thin, node colors rectangle for module genes, circle for enriched TFs and diamond for non-enriched TFs etc.

To load it form cytoscape app go: File -> import -> styles from file Loaded input files to Cyoscape using cytoscape executive files:

To load the executive files from the app select : Tools -> Execute command files. Note: Each time before loading the executive for new module you have to adjust the min-max node size from the app setting it to minimum and maximum number of edges for nodes in the module specific subnetwork.


