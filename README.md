## Biological Random Walks: Integrating Gene Expression for Tissue Specific Prediction

### Citation 
M. Gentili, L. Martini, L. Becchetti, M. Sponziello, Biological Random Walks: integrating gene expression in tissue-specific prediction. 

### Author: 

Michele Gentili, Leonardo Martini, Luca Becchetti, Marialuisa Sponziello

Contact:

Michele Gentili (gentili@diag.uniroma1.it),  Leonardo Martini (martini@diag.uniroma1.it)

Department of Computer Science

The University of Rome "La Sapienza", Rome, Italy

### Description

### Libraries

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


