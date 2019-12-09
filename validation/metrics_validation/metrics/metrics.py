import math

def DCG(ranked_list,test_set, disease_module_size):
    DCG_by_k = []
    test_set = set(test_set)
    ranked_list = ranked_list[:max(disease_module_size)]
    DCG = 0
    for index, row in enumerate(ranked_list):
        rank = index+1
        if row in test_set:
            DCG += 1 / math.log(rank + 1, 2)
        DCG_by_k.append(DCG)
    DCG_by_k = [DCG_by_k[i-1] for i in disease_module_size]
    return DCG_by_k



def recall_at_k(ranked_list,test_set, disease_module_size):
    recall_at_k = []
    test_set = set(test_set)
    ranked_list = ranked_list[:max(disease_module_size)]
    for k in disease_module_size:
        recall_at_k.append(len(test_set.intersection(ranked_list[:k]))/ len(test_set))


    return recall_at_k

