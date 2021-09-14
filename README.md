## Biological Random Walks: Integrating Gene Expression for Tissue Specific Prediction

### Citation 
[1] M. Gentili, L. Martini, L. Becchetti, M. Sponziello, Biological Random Walks: integrating gene expression in tissue-specific prediction. (2021)

[2] M. Gentili, L. Martini, M. Petti, L. Farina, L. Becchetti, Biological Random Walks: integrating heterogeneous data in disease gene prioritization. IEEE Symposium on Computational Intelligence and Bioinformatics and Computational Biology (CIBCB) (2019)

### Author: 

Michele Gentili, Leonardo Martini, Manuela Petti, Luca Becchetti, Lorenzo Farina, Marialuisa Sponziello

Contact:

Michele Gentili (gentili@diag.uniroma1.it),  Leonardo Martini (martini@diag.uniroma1.it)

Department of Computer, Control, and Management Engineering Antonio Ruberti, Sapienza University of Rome, Rome, Italy

Translational and Precision Medicine Department Sapienza University of Rome, Rome, Italy

### Algorithm Description

The  Biological  Random  Walk  heuristic  provides  a  framework  to  integrateheterogeneous  biological  data  sources  within  diffusion-based  prioritizationmethods that are based on the well known Random Walk with restart algo-rithm (RWR). For the sake of exposition, in the remainder we refer to thebiological information associated to a gene i (e.g., the set of its annotations)as  the  set  of  its  labels,  denoted  by  labels(i). BRW ranks genes according to the following steps:
- We compute the set of statistically significant annotations of known disease gene
- Rather than using the standard method, we compute individual tele-porting  probabilities  for  all  nodes  of  the  PPI.  In  particular,  the  Biological Teleporting Probability (BTP) of a node increases with thesimilarity between its labels and the enriched set
- In a similar fashion,  we  weigh  PPI  network  interactions  using  nodeannotations and the enriched set.  This results in a modified randomwalk, namely the Biological Random Walk (BRW), in which flow prop-agation  is  biased  toward  genes  that  are  functionally  closer  to  thoseforming the seed set
- Finally,  we  rank  genes  according  to  their  Biological  Random  Walk(BRW) score.


![alt text](https://github.com/LeoM93/BiologicalRandomWalks/blob/master/imgs/BRW_flow.png?raw=true)

Biological Random  Walks flow propagation: given the seed nodes (star nodes), the flow propagates to his neighbors. The BRW not only propagatethe  flow  around  them  but  also  teleports  the  flow  to  the  target  of  the  BTP  nodes  (blue  arrows).  So  it  discovers  nodes  that  are  biologically  correlated  to  theseed nodes (just through the BTP, left-lower test node) and those nodes that aren’t reached directly to the BTP but are close to many related nodes.

### Python Libraries
Imported libraries and their version:

- numpy, version 1.19.1
- networkx, version 2.4
- sklearn, version 0.23.1 

### Toy Examples

Please go the the directory BiologicalRqndomWalks/toy_example to see the inputs format of BRW. 

### Input files and formats

The directory BiologicalRqndomWalks/toy_example contains input files as a toy dataset for BRW.

They are:
 - ppi_network.tsv: A Protein-Protein Interaction (PPI) network in tab separated format 
 ```
 < node_id_1 > \t < node_id_2 > 
 ```
 - co_expression_network.tsv: A Weighted Co-Expression network in tab separated format
 ```
 < node_id_1> \t < node_id_2 > \t < score >
 ```
 - seed_set.tsv: A list of node ids (each row contains one id)
 - de_genes.tsv: A list of Differentially Expressed (DE) gene (each row contains one id)
 - annotations.tsv: Gene-annotation associantions dataset in tab separated format
 ```
 < node_id_1 > \t < annotation_id > \t < dataset_name >
 ```
 - disease_ontology.tsv: Statistically significant Disease-annotation associantions dataset in tab separated format
 ```
 < annotation_id > \t < dataset_name >
 ```

### How to Install and run Biological Random Walks


To install the framework download the zip file and decompress it. Then, Iin order to run the algorithm, from terminal go in the BRW directory and type "python3 main.py" followed by option:

 -p < protein_interaction_network_file_path > (required)

 -c < co_expression_network_file_path > (optional)

 -s < seed_set_file_path > (required) 

 -de < differentially_expressed_genes_file_path > (optional, default: None) 

 -do < disease_ontologies_file_path > (required)
 
 -a < gene_ontology_assotiation_dataset_file_path> (required)
 
 -r < restart probability > (optional, default: 0.75)
 
 -o < output_file_path > (required)

To Run the algorithm described in [1], all the options are required,to choose the other one described in [2] do not use options -c and -de

### Example

 ```
 python3 main.py -p toy_example/ppi_network.tsv -c toy_example/co_expression_network.tsv -s toy_example/seed_set.tsv -de toy_example/de_genes.tsv -a toy_example/annotations.tsv -do toy_example/disease_ontology.tsv -o toy_example/output.tsv
 ```


