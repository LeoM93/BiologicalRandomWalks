import unicodecsv as csv
import json
import os



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
    print(file_path)
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



###############################################  Gene Expression  ###########################################
def load_gene_expression_z_score(file_path):

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
                        gene_name = item
                    else:
                        try:
                            z_score_dict[gene_name].append(float(item))
                        except KeyError:
                            z_score_dict[gene_name]=[float(item)]

    return patient_list, z_score_dict


def load_co_expression_network(file_path, is_direct):

    co_expression_network = {}

    with open(file_path,'rb') as file:
        csv_reader = csv.reader(file,delimiter=',')
        for index,row in enumerate(csv_reader):
            if index == 0:
                continue

            gene_1 = row[0]
            gene_2 = row[1]

            pearson_correlation = float(row[2])
            p_value = float(row[3])

            try:
                co_expression_network[gene_1][gene_2] = (pearson_correlation,p_value)
            except KeyError:
                co_expression_network[gene_1] = {gene_2:(pearson_correlation,p_value)}


            if is_direct is False:

                try:
                    co_expression_network[gene_2][gene_1] = (pearson_correlation, p_value)
                except KeyError:
                    co_expression_network[gene_2] = {gene_1: (pearson_correlation, p_value)}

    return co_expression_network


############################################### End  Gene Expression  ##################################################
