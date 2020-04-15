import unicodecsv as csv
import json
import os
import operator

from utils.preferences import NETWORK_MEDICINE_PATH


class LoadAlgorithmVariables():

    def __init__(self):
        pass


    def load_algorithm_params(self,experiment_variable_directory_path):
        files = os.listdir(experiment_variable_directory_path)

        experiment_variables_dictionary = {}

        for file in files:
            if ".json" in file:

                with open(experiment_variable_directory_path + file ) as f:
                    experiment_variables = json.load(f)

                for experiment_variable in experiment_variables:
                    algorithm = experiment_variable['algorithm']

                    try:
                        experiment_variables_dictionary[algorithm].append(experiment_variable)

                    except KeyError:
                        experiment_variables_dictionary[algorithm] = [experiment_variable]


        return experiment_variables_dictionary

############################################## Disease Feature #########################################################
def load_train_test_file(file_path):
    with open(file_path, "rb") as fp:
        tsv_reader = csv.reader(fp, delimiter=",")
        train_set = [x[0] for x in tsv_reader]

    return train_set
############################################## End Disease Feature #####################################################

############################################## Validation #############################################################
def load_algorithm_ranked_list(file_path, seed_set,delimiter =","):
    algorithm_output = []
    with open(file_path, 'rb') as file:

        tsv_reader = csv.reader(file, delimiter=delimiter)
        for id, row in enumerate(tsv_reader):
            if id == 0:
                continue

            gene_name = row[0]
            if gene_name not in seed_set:
                algorithm_output.append(gene_name)

    return algorithm_output

############################################## Load Enrichment Analysis#################################################
def load_disease_enrichment_analysis(file_path):
    p_value_by_term_id = {}

    with open(file_path, 'rb') as file:

        tsv_reader = csv.reader(file, delimiter="\t")
        for id, row in enumerate(tsv_reader):

            if id == 0:
                continue

            try:
                term_id = int(row[0])
            except ValueError:
                term_id = row[0]

            p_value = float(row[1])
            odds_ratio = float(row[2])


            p_value_by_term_id[term_id] = [p_value,odds_ratio]

    return p_value_by_term_id



def load_validation_data_frame(file_path):

    column_id_to_name = {}
    algorithm_validation = {}
    with open(file_path,"rb") as fp:
        tsv_reader = csv.reader(fp, delimiter=",")
        for id, row in enumerate(tsv_reader):

            if id == 0:
                for index, item in enumerate(row):
                    column_id_to_name[index] = item.split("-")[0]
            else:
                k = int(row[0])
                for i in range(len(row)):

                    if i == 0:
                        continue

                    try:
                        algorithm_validation[column_id_to_name[i]][k] = float(row[i])
                    except KeyError:
                        algorithm_validation[column_id_to_name[i]] = {k:float(row[i])}




    return algorithm_validation










###############################################  Gene Expression  ###########################################
def load_gene_expression_z_score(file_path,gene_to_ensembl = None):

    patient_list = []
    z_score_dict = {}

    with open(file_path,'rb') as file:
        csv_reader = csv.reader(file,delimiter=',')
        for index,row in enumerate(csv_reader):
            if index == 0:
                for index_col,item in enumerate(row):
                    patient_list.append(item)
            else:
                gene_name = ""
                for index_col, item in enumerate(row):
                    if index_col == 0:
                        if gene_to_ensembl is None:
                            gene_name = item
                        else:
                            try:
                                gene_name = gene_to_ensembl[item]
                            except KeyError:
                                pass
                    else:
                        try:
                            z_score_dict[gene_name].append(float(item))
                        except KeyError:
                            z_score_dict[gene_name]=[float(item)]

    return patient_list, z_score_dict


def filter_co_expression_network(file_path, node_list, pearson_correlation_threshold):

    co_expression_adjacency_matrix = []

    with open(file_path,'rb') as file:
        csv_reader = csv.reader(file,delimiter=',')

        for index,row in enumerate(csv_reader):
            if index == 0:
                continue

            gene_1 = row[0]
            gene_2 = row[1]

            module_pearson_correlation = abs(float(row[2]))

            if gene_1 in node_list and gene_2 in node_list:

                if module_pearson_correlation > pearson_correlation_threshold:
                    co_expression_adjacency_matrix.append([gene_1,gene_2,module_pearson_correlation])



    return co_expression_adjacency_matrix


def load_co_expression_network(file_path):

    gene_to_gene_ge_weight = {}

    with open(file_path,'rb') as file:
        csv_reader = csv.reader(file,delimiter=',')

        for index,row in enumerate(csv_reader):
            if index == 0:
                continue

            gene_1 = row[0]
            gene_2 = row[1]
            pearson_correlation = float(row[2])



            if row[2] == "nan":
                pearson_correlation = 0.0

            try:
                gene_to_gene_ge_weight[gene_1][gene_2] = pearson_correlation
            except KeyError:
                gene_to_gene_ge_weight[gene_1] = { gene_2: pearson_correlation}

            try:
                gene_to_gene_ge_weight[gene_2][gene_1] = pearson_correlation
            except KeyError:
                gene_to_gene_ge_weight[gene_2] = { gene_1: pearson_correlation}

    return gene_to_gene_ge_weight

def load_independent_t_test(file_path):

    independent_t_test_dict = {}
    ranking_by_gene_name = {}

    with open(file_path,'rb') as file:
        csv_reader = csv.reader(file,delimiter=',')
        for index,row in enumerate(csv_reader):

            if index == 0:
                continue

            gene_name = row[0]
            p_value = float(row[1])
            mean = float(row[2])

            independent_t_test_dict[gene_name] = p_value

    sorted_independent_t_test = sorted(independent_t_test_dict.items(), key=operator.itemgetter(1))

    for index, record in enumerate(sorted_independent_t_test):
        gene_name = record[0]

        ranking_by_gene_name[gene_name] = index + 1

    return independent_t_test_dict,ranking_by_gene_name


############################################### End  Gene Expression  ##################################################


def load_gene_ensembl_id_from_disk( delimiter=','):

    file_path = NETWORK_MEDICINE_PATH + "db/ensembl_mapping/ensembl_ppi_barabasi_data_set.txt"
    gene_to_ensembl_id = {}
    ensembl_id_to_gene = {}
    with open(file_path, 'rb') as fp:
        csv_reader = csv.reader(fp, delimiter=delimiter)

        for row in csv_reader:
            gene_name = row[0]
            ensembl_id = row[1]

            gene_to_ensembl_id[gene_name] = ensembl_id
            ensembl_id_to_gene[ensembl_id] = gene_name

    return gene_to_ensembl_id,ensembl_id_to_gene