## Biological Random Walks: Integrating Gene Expression for Tissue Specific Prediction

### Citation 
[1] M. Gentili, L. Martini, L. Becchetti, M. Sponziello, Biological Random Walks: integrating gene expression in tissue-specific prediction. 

[2] M. Gentili, L. Martini, M. Petti, L. Farina, L. Becchetti, Biological Random Walks: integrating gene expression in tissue-specific prediction. IEEE Symposium on Computational Intelligence and Bioinformatics and Computational Biology (CIBCB) (2019)

### Author: 

Michele Gentili, Leonardo Martini, Luca Becchetti, Marialuisa Sponziello

Contact:
Michele Gentili (gentili@diag.uniroma1.it),  Leonardo Martini (martini@diag.uniroma1.it)

Department of Computer Science
The University of Rome "La Sapienza", Rome, Italy

### Algorithm Description

Like Random Walks with Restart (RWR) and diffusion methods in general, Biological Random Walks assumes knowledge of a seed set of known genes for a disease of interest. Unlike previous ones, BRW consists of two main steps:

- extracting statistically significant features from biological data, using themto compute a personalization vector and a transition matrix used by the algorithm;  - using the stationary distribution of the correspondingrandom walk to rank genes

![alt text](https://github.com/LeoM93/BiologicalRandomWalks/blob/master/imgs/BRW_flow.png?raw=true)


### Python Libraries

### Toy Examples

Please go the the directory BiologicalRqndomWalks/toy_example to see the inputs format of BRW. 

### Input files

The directory BiologicalRqndomWalks/toy_example contains input files as a toy dataset for BRW.

They are:
 - Protein-Protein Interaction (PPI) network: ppi_network.tsv
 - Weighted Co-Expression network: co_expression_network.tsv
 - Seed set: seed_set.tsv
 - Differentially Expressed (DE) gene set: de_genes.tsv
 - Gene-annotation associantions dataset: annotations.tsv
 - Disease-annotation associantions dataset: disease_ontology.tsv


### Input files format

### How to run BRW


### Example

 ```
 bin/compute_HSPs -p toy_example/protein_sequences.seq -h HSP/hsps_toy_example 
 ```


