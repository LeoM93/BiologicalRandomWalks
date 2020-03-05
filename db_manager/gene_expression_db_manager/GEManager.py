from db_manager.db_manager import DBManager
from utils.preferences import NETWORK_MEDICINE_PATH
from file_manager.file_writer import write_data_on_disk,write_row_on_disk
from file_manager.file_loader import load_gene_expression_z_score,filter_co_expression_network,load_co_expression_network,load_independent_t_test

import math
from scipy.stats import pearsonr
import scipy

import os
class GEManager():

    def __init__(self,table,ppi_name, pearson_correlation_threshold = 0.8):

        self.db_table = table
        self.name = table
        self.gene_expression_db_path = NETWORK_MEDICINE_PATH + "db/gene_expression.db"
        self.ge_db_manager = DBManager(self.gene_expression_db_path)

        self.ppi_name = ppi_name

        self.z_score_file_path = NETWORK_MEDICINE_PATH + "db/gene_expression/" + self.db_table + "_z_score.csv"

        self.differential_expression_file_path = NETWORK_MEDICINE_PATH + "db/gene_expression/" + self.db_table + "_differential_expression_independent_t_test.csv"

        self.co_expression_by_z_score_file_path = NETWORK_MEDICINE_PATH + "db/gene_expression/" + self.db_table + "_co_expression_network_pearson_correlation_by_z_score.csv"

        self.pearson_correlation_threshold = pearson_correlation_threshold

    def print_schema(self):
        self.ge_db_manager.print_db_schema()


    def init_db(self,node_list):

        self.load_z_score()
        self.load_differential_expression_statistics()
        self.write_co_expression_network(node_list)


    def close_connection(self):
        self.ge_db_manager.close_connection()
        del self.ge_db_manager


    def load_z_score(self,):

        if os.path.exists(self.z_score_file_path):
            self.patients, self.z_score = load_gene_expression_z_score(self.z_score_file_path)
        else:
            self._generate_z_score()
            self.patients, self.z_score = load_gene_expression_z_score(self.z_score_file_path)


    def load_differential_expression_statistics(self):

        if os.path.exists(self.differential_expression_file_path):
            self.independent_t_test,self.rank_by_t_test = load_independent_t_test(self.differential_expression_file_path)
        else:
            self._generate_independent_t_test()




    def write_co_expression_network(self,node_list):

        if os.path.exists(self.co_expression_by_z_score_file_path):

            output_file_path = NETWORK_MEDICINE_PATH + "db/gene_expression/co_expression_network_" + self.ppi_name +"_" + self.db_table + "_threshold_" + str(self.pearson_correlation_threshold) +".csv"
            if os.path.exists(output_file_path):
                pass
            else:
                filtered_co_expression_adjacency_matrix = filter_co_expression_network(self.co_expression_by_z_score_file_path,node_list,self.pearson_correlation_threshold)
                write_data_on_disk(output_file_path,rows=filtered_co_expression_adjacency_matrix, headers=["Gene_1","Gene_2","correlation"],write_mode="wb")

        else:
            self._generate_total_co_expression_network()



    def _generate_total_co_expression_network(self):


        keys = set()

        for index, gene_1 in enumerate(self.z_score):
            print(index)
            for gene_2 in self.z_score:
                if gene_1 != gene_2:
                    key = (gene_1, gene_2)

                    sorted_key = sorted(key)
                    tuple_sorted_key = (sorted_key[0],sorted_key[1])

                    if tuple_sorted_key in keys:
                        continue

                    keys.add(tuple_sorted_key)
                    pearson = pearsonr(self.z_score[tuple_sorted_key[0]], self.z_score[tuple_sorted_key[1]])

                    if math.isnan(pearson[0]):
                        pass

                    else:
                        record = [tuple_sorted_key[0],tuple_sorted_key[1], pearson[0],pearson[1]]
                        write_row_on_disk(self.co_expression_by_z_score_file_path,
                                           row=record, write_mode="ab")




    def load_weighted_gene_expression_adjacency_matrix(self, threshold):

        file_path = NETWORK_MEDICINE_PATH + "db/gene_expression/co_expression_network_" + self.ppi_name + "_" + self.db_table + "_threshold_" + str(
            threshold) + ".csv"

        if os.path.exists(file_path):
            return load_co_expression_network(file_path)
        else:
            print("THRESHOLDED CO EXPRESSION NETWORK DOES NOT EXIST! ")
            exit(12)




    def find_gene_z_score(self,gene):
        if gene in self.z_score:

            return self.z_score[gene]
        else:
            return -1

    def find_gene_ranking_by_independent_t_test(self):

        return -1


    def _generate_z_score(self):


        if self.db_table != 'TGCA_2014':
            print()
            print("IT IS POSSIBLE TO USE THIS OPTION ONLY WITH " + self.db_table + " TABLE")
            exit(12)

        gene_expression_query_selection = "SELECT * FROM " + self.db_table
        tcga_table = self.ge_db_manager.comupute_query(gene_expression_query_selection,record=None)

        gene_expression_by_tumor_patients = {}
        gene_expression_by_sane_patients = {}


        for item in tcga_table:

            gene_name = item[0]
            gene_expression = item[1]
            cell_type = item[2]
            patient_id = item[3]


            if cell_type == "TUMOR":
                try:
                    gene_expression_by_tumor_patients[patient_id][gene_name] = gene_expression
                except KeyError:
                    gene_expression_by_tumor_patients[patient_id] = {gene_name: gene_expression}
            else:
                try:
                    gene_expression_by_sane_patients[patient_id][gene_name] = gene_expression
                except KeyError:
                    gene_expression_by_sane_patients[patient_id] = {gene_name: gene_expression}

        z_score_by_patient_id = {}

        # mean and standard deviation of the control group
        mean_by_gene_name = {}
        standard_deviation_by_gene_name = {}

        num_sane_patients = len(gene_expression_by_sane_patients.keys())

        for patient, gene_expression in gene_expression_by_sane_patients.items():

            for gene,value in gene_expression.items():
                try:
                    mean_by_gene_name[gene] += value
                except KeyError:
                    mean_by_gene_name[gene] = value


        for k,v in mean_by_gene_name.items():
            mean_by_gene_name[k] = mean_by_gene_name[k]/num_sane_patients

        for patient, gene_expression in gene_expression_by_sane_patients.items():
            for gene,value in gene_expression.items():

                try:
                    standard_deviation_by_gene_name[gene] += (value - mean_by_gene_name[gene])**2
                except KeyError:
                    standard_deviation_by_gene_name[gene] = (value - mean_by_gene_name[gene])**2

        for k,v in mean_by_gene_name.items():
            standard_deviation_by_gene_name[k] = (standard_deviation_by_gene_name[k]/num_sane_patients)**0.5

        gene_names = set()
        patient_ids = list(gene_expression_by_tumor_patients.keys())

        for patient, gene_expression in gene_expression_by_tumor_patients.items():

            for gene,value in gene_expression.items():


                if standard_deviation_by_gene_name[gene] != 0.0:
                    z_score = (value - mean_by_gene_name[gene]) / standard_deviation_by_gene_name[gene]
                    gene_names.add(gene)

                    try:
                        z_score_by_patient_id[patient][gene] = z_score


                    except KeyError:
                        z_score_by_patient_id[patient] = {gene: z_score}



        gene_names = list(gene_names)
        z_score_table = []

        for gene in gene_names:
            record = [gene]

            for patient in patient_ids:
                try:
                    record.append(z_score_by_patient_id[patient][gene])
                except KeyError:
                    record.append(0.0)

            z_score_table.append(record)

        write_data_on_disk(self.z_score_file_path,headers=patient_ids,rows=z_score_table)


    def _generate_independent_t_test(self):

        if self.db_table != 'TGCA_2014':
            print()
            print("IT IS POSSIBLE TO USE THIS OPTION ONLY WITH " + self.db_table + " TABLE")
            exit(12)

        gene_expression_query_selection = "SELECT * FROM " + self.db_table
        tcga_table = self.ge_db_manager.comupute_query(gene_expression_query_selection, record=None)

        gene_expression_by_tumor_patients = {}
        gene_expression_by_sane_patients = {}

        for item in tcga_table:

            gene_name = item[0]
            gene_expression = item[1]
            cell_type = item[2]
            patient_id = item[3]

            if cell_type == "TUMOR":
                try:
                    gene_expression_by_tumor_patients[patient_id][gene_name] = gene_expression
                except KeyError:
                    gene_expression_by_tumor_patients[patient_id] = {gene_name: gene_expression}
            else:
                try:
                    gene_expression_by_sane_patients[patient_id][gene_name] = gene_expression
                except KeyError:
                    gene_expression_by_sane_patients[patient_id] = {gene_name: gene_expression}

        mean_by_sane_gene_name = {}
        mean_by_tumor_gene_name = {}

        healthy_group_gene_expression_by_gene_name = {}
        tumor_group_gene_expression_by_gene_name = {}

        x_i = {}

        num_sane_patients = len(gene_expression_by_sane_patients.keys())
        num_tumor_patients = len(gene_expression_by_tumor_patients)

        for patient, gene_expression in gene_expression_by_sane_patients.items():

            for gene, value in gene_expression.items():
                try:

                    mean_by_sane_gene_name[gene] += value
                    healthy_group_gene_expression_by_gene_name[gene].append(value)

                except KeyError:

                    healthy_group_gene_expression_by_gene_name[gene] = [value]
                    mean_by_sane_gene_name[gene] = value


        for patient, gene_expression in gene_expression_by_tumor_patients.items():

            for gene, value in gene_expression.items():
                try:
                    mean_by_tumor_gene_name[gene] += value
                    tumor_group_gene_expression_by_gene_name[gene].append(value)
                except KeyError:
                    tumor_group_gene_expression_by_gene_name[gene] = [value]
                    mean_by_tumor_gene_name[gene] = value


        for gene in mean_by_tumor_gene_name.keys():
            mean_by_tumor_gene_name[gene] = mean_by_tumor_gene_name[gene]/num_tumor_patients

        gene_table = []
        for gene, value in tumor_group_gene_expression_by_gene_name.items():
            p_value = scipy.stats.ttest_ind(healthy_group_gene_expression_by_gene_name[gene], tumor_group_gene_expression_by_gene_name[gene], axis=0, equal_var=False)[1]
            gene_table.append([gene, p_value, mean_by_tumor_gene_name[gene]])

        write_data_on_disk(self.differential_expression_file_path,headers=["Gene","P-Value","Mean"],rows=gene_table)
        

