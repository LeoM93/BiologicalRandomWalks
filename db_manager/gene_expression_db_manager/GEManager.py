from db_manager.db_manager import DBManager
from utils.preferences import NETWORK_MEDICINE_PATH
from file_manager.file_writer import write_data_on_disk,write_row_on_disk
from file_manager.file_loader import load_gene_expression_z_score,filter_co_expression_network,load_co_expression_network,load_independent_t_test

import math
from scipy.stats import pearsonr
import numpy as np
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
        self.sub_sample_co_expression_by_z_score_file_path = NETWORK_MEDICINE_PATH + "db/gene_expression/" + self.db_table + "_" + ppi_name + "_sub_sampled_co_expression_network_pearson_correlation_by_z_score.csv"
        self.pearson_correlation_threshold = pearson_correlation_threshold

    def print_schema(self):
        self.ge_db_manager.print_db_schema()


    def init_db(self,node_list):
        self.load_z_score()


    def close_connection(self):
        self.ge_db_manager.close_connection()
        del self.ge_db_manager


    def load_z_score(self,):

        if os.path.exists(self.z_score_file_path):
            self.patients, self.z_score = load_gene_expression_z_score(self.z_score_file_path)
        else:
            self._generate_z_score()
            self.patients, self.z_score = load_gene_expression_z_score(self.z_score_file_path)












    def load_weighted_gene_expression_adjacency_matrix(self):


        file_path = NETWORK_MEDICINE_PATH + "db/gene_expression/co_expression_network_" + self.ppi_name + "_" + self.db_table + "_down_sampled.csv"

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
