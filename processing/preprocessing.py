import random
from utils.disease import Disease

def check_input_disease(disease,node_list,):
    # use as seed nodes only the nodes that have a mapping in the ppi network
    checked_disease_genes = []
    unmapped_genes = []

    for gene in disease.disease_genes:

        if gene in node_list:
            checked_disease_genes.append(gene)
        else:
            unmapped_genes.append(gene)


    return Disease(disease.disease_id, disease.disease_name,checked_disease_genes, disease.disease_chromosomes,disease.disease_class)



def split_train_test_set(nodes,train_percentuage = 0.7):
    train_set = set()
    test_set = []
    train_length = int(len(nodes) * train_percentuage)
    node_len = len(nodes)
    while len(train_set) != train_length:
        random_number = random.randint(0, node_len -1)

        train_set.add(nodes[random_number])

    for node in nodes:
        if node not in train_set:
            test_set.append(node)

    return list(train_set), test_set
