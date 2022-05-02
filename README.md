# Biological Random Walks: Integrating Gene Expression for Tissue Specific Prediction

## Citation 
[1] M. Gentili, L. Martini, L. Becchetti, M. Sponziello, Biological Random Walks: integrating gene expression in tissue-specific prediction. (2021)

[2] M. Gentili, L. Martini, M. Petti, L. Farina, L. Becchetti, Biological Random Walks: integrating heterogeneous data in disease gene prioritization. IEEE Symposium on Computational Intelligence and Bioinformatics and Computational Biology (CIBCB) (2019)

## Author: 

Michele Gentili, Leonardo Martini, Manuela Petti, Luca Becchetti, Lorenzo Farina, Marialuisa Sponziello

Contact:

Michele Gentili (gentili@diag.uniroma1.it),  Leonardo Martini (martini@diag.uniroma1.it)

Department of Computer, Control, and Management Engineering Antonio Ruberti, Sapienza University of Rome, Rome, Italy

Translational and Precision Medicine Department Sapienza University of Rome, Rome, Italy

## Algorithm Description

The  Biological  Random  Walk  heuristic  provides  a  framework  to  integrateheterogeneous  biological  data  sources  within  diffusion-based  prioritizationmethods that are based on the well known Random Walk with restart algo-rithm (RWR). For the sake of exposition, in the remainder we refer to thebiological information associated to a gene i (e.g., the set of its annotations)as  the  set  of  its  labels,  denoted  by  labels(i). BRW ranks genes according to the following steps:
- We compute the set of statistically significant annotations of known disease gene
- Rather than using the standard method, we compute individual tele-porting  probabilities  for  all  nodes  of  the  PPI.  In  particular,  the  Biological Teleporting Probability (BTP) of a node increases with thesimilarity between its labels and the enriched set
- In a similar fashion,  we  weigh  PPI  network  interactions  using  nodeannotations and the enriched set.  This results in a modified randomwalk, namely the Biological Random Walk (BRW), in which flow prop-agation  is  biased  toward  genes  that  are  functionally  closer  to  thoseforming the seed set
- Finally,  we  rank  genes  according  to  their  Biological  Random  Walk(BRW) score.


![alt text](https://github.com/LeoM93/BiologicalRandomWalks/blob/master/imgs/BRW_flow.png?raw=true)

Biological Random  Walks flow propagation: given the seed nodes (star nodes), the flow propagates to his neighbors. The BRW not only propagatethe  flow  around  them  but  also  teleports  the  flow  to  the  target  of  the  BTP  nodes  (blue  arrows).  So  it  discovers  nodes  that  are  biologically  correlated  to  theseed nodes (just through the BTP, left-lower test node) and those nodes that aren’t reached directly to the BTP but are close to many related nodes.

## Python Libraries
Imported libraries and their version:

- numpy, version 1.19.1
- networkx, version 2.4
- sklearn, version 0.23.1 



## Data Preprocessing

### Computing the Ontology Graph

To compute the ontology graph, run the following commands:

 ```
 cd data_preprocessing
 ```
 
 ```
 python3 compute_ontology_graph.py -go <path to .gaf file> -r <path to Reactome file> -k <path to KEGG file> -o <output file path>  
 ```
The file used in the manuscript can be downloaded at: https://drive.google.com/file/d/12oDaaEs1vso82UXsRe2AWeoGqNccZuLM/view?usp=sharing
An update version of the .gaf file can be downloaded from Gene Onotology Consortium at the following link: http://current.geneontology.org/products/pages/downloads.html

The Reactome file can be downloaded at https://reactome.org/download-data (Uniprot to all pathways)
The KEGG file has been downloaded using KEGG rest api.

The ontology graph has the following structre:
```
 < ensembl_id > \t < annotation_id > \t < dataset_name >
```


### Computing the Tumor-Control Table TCGA

To compute Tumor and Control Table for each Tumor taken in consideration in the article, run the following commands:

 ```
 cd data_preprocessing
 ```
 
 ```
 python3 TCGA_analyzer.py -gdc <path to gdc sheet file> -m <path to manifest file> -rna_dir <path rna dir downloaded using cdc-client> -o <output dir path>  
```
The gdc sheet and manifest files that we have used in the manuscript can be found in the data_set/TCGA directory. The rna_seq file downloaded using cdc-client with the gdc sheet and manifest files is to heavy to be uploaded here. Thus, It can be found at the following link: https://drive.google.com/file/d/1f2V6fji8dPshiH6ohV81K0Cxv9q5D4ew/view?usp=sharing


### Computing the Co-expression network and Differentially Expressed Genes

Once Tumor and control Table have been created following the steps described above, we can generate the co-expression network and the DE genes using the following commands:
 ```
 cd data_preprocessing
 ```
 
  ```
 python3 compute_co_expression_and_de_genes.py -T <path to TCGA Tumor Table> -C <path to TCGA Control Table> -de <de output file path> -co < co-expression network output file path>  
```




## Running Biological Random Walks

### Input files and formats

The directory BiologicalRqndomWalks/toy_example contains input files as a toy dataset for BRW.

They are:
 - ppi_network.tsv: A Protein-Protein Interaction (PPI) network in tab separated format 
 ```
 < ensembl_id_1 > \t < ensembl_id_2 > 
 ```
 - co_expression_network.tsv: A Weighted Co-Expression network in tab separated format
 ```
 < ensembl_id_1> \t < ensembl_id_2 > \t < score >
 ```
 - seed_set.tsv: A list genes  (each row contains one ensembl id)
 - de_genes.tsv: A list of Differentially Expressed (DE) gene (each row contains one ensembl id)
 - annotations.tsv: Gene-annotation associantions dataset in tab separated format generated by compute_ontology_graph.py
 ```
 < ensembl_id > \t < annotation_id > \t < dataset_name >
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
 python3 main.py -p <ppi network file path> -c <path to co-expressiion network> -s <path to seed set> -de <path to de genes> -a <path to ontology graph> -do <path to seed enriched ontologies> -o <output file path>
 ```


