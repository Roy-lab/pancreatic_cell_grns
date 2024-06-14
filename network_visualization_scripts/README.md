**Steps of generating module specific network files**

**Step1: Generatin adjacency matrices for module specific subnetwork for each sample using convertNet.py script:**

Example run: *convertNet.py input/WT_rep1_expression_matrix.txt input/module1_subnetwork.txt output/adj_WT_rep1_expression_matrix.txt output/WT_rep1_expression_matrix.txt*

**Step2: Running matlab script runls_beta.m and runcc_beta.m to calculate regression and correlation coeefficient using outputs in step1:**

**Step3: Generating input files for Cytoscape using makeCyto2.py script.** 

This script takes 7 imput files:
1. reg_files.txt: a file containing path to all files contains regression coefficient produced in previous step by runls.m script (this should also include path to the module specific subnetwork).
2. cc_files.txt: a file containing path to all files containing correlation coefficient produced in previous step by runcc.m script (this should also include path to the module specific subnetwork)
3. zeromean_expression_files.txt: a file containing path to all genotype specific files which contain mean expression of zeromean normalized matrix per genotype.
4. deg.txt: two column file showing all the genes in module specific subnetwork as first column and number of outgoing edges they have in subnetwork as second columns
5. module.txt: two column file showing all the genes in module specific subnetwork as first column and module number in second column
6. condition_names.txt:a file containing name of the cell types/conditions
7. in: output files prefix

Example run:

*python makeCyto2.py reg_files.txt  input/cc_files.txt input/zeromean_expression_files.txt input/deg.txt input/module.txt input/condition_names.txt cytoscape_inputs/in; done*
