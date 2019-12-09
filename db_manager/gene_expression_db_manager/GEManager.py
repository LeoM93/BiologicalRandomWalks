from db_manager.db_manager import DBManager
from utils.preferences import NETWORK_MEDICINE_PATH
from file_manager.file_writer import write_data_on_disk,write_row_on_disk
from file_manager.file_loader import load_gene_expression_z_score,load_co_expression_network

import math
from scipy.stats import pearsonr

import os
class GEManager():

    def __init__(self,table,ppi_name):

        self.db_table = table
        self.name = table
        self.gene_expression_db_path = NETWORK_MEDICINE_PATH + "db/gene_expression.db"
        self.ge_db_manager = DBManager(self.gene_expression_db_path)

        self.ppi_name = ppi_name

        self.z_score_file_path = NETWORK_MEDICINE_PATH + "db/gene_expression/" + self.db_table + "_z_score.csv"
        self.co_expression_by_z_score_file_path = NETWORK_MEDICINE_PATH + "db/gene_expression/" + self.db_table + "_co_expression_network_pearson_correlation_by_z_score.csv"

        self.ppi_filtered_co_expression_by_z_score_file_path = NETWORK_MEDICINE_PATH + "db/gene_expression/" + self.db_table +"_" + self.ppi_name + "_pearson_correlation_by_z_score.csv"


    def print_schema(self):
        self.ge_db_manager.print_db_schema()


    def init_db(self,edge_list):

        self.load_z_score()
        self.load_co_expression_network(edge_list)


    def close_connection(self):
        self.ge_db_manager.close_connection()
        del self.ge_db_manager


    def load_z_score(self):

        if os.path.exists(self.z_score_file_path):
            self.patients, self.z_score = load_gene_expression_z_score(self.z_score_file_path)
        else:
            self._generate_z_score()
            self.patients, self.z_score = load_gene_expression_z_score(self.z_score_file_path)



    def load_co_expression_network(self,edge_list):

        if os.path.exists(self.ppi_filtered_co_expression_by_z_score_file_path):
            self.ppi_filtered_co_expression_network = load_co_expression_network(self.ppi_filtered_co_expression_by_z_score_file_path,is_direct=False)
        else:
            self._generate_partial_co_expression_network(edge_list)

        if os.path.exists(self.co_expression_by_z_score_file_path):
            pass
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




    def _generate_partial_co_expression_network(self,edge_list):

        co_expression_network_table = []

        edge_set = set(edge_list)

        for index,gene_1 in enumerate(self.z_score):
            for gene_2 in self.z_score:
                    if gene_1 != gene_2:
                        key = (gene_1,gene_2)

                        if key in edge_set:
                            pearson = pearsonr( self.z_score[gene_1],self.z_score[gene_2])

                            if math.isnan(pearson[0]):
                                pass

                            else:
                                co_expression_network_table.append([key[0],key[1],pearson[0],pearson[1]])

        write_data_on_disk(self.ppi_filtered_co_expression_by_z_score_file_path,headers=["Gene_1","Gene_2","Pearson","P_Value"],rows=co_expression_network_table)




    def find_gene_z_score(self,gene):
        if gene in self.z_score:

            return self.z_score[gene]
        else:
            return -1

    def find_gene_pearson_correlation_by_ppi_network(self,gene,neighbor):

        if gene in self.ppi_filtered_co_expression_network:
            if neighbor in self.ppi_filtered_co_expression_network[gene]:
                return self.ppi_filtered_co_expression_network[gene][neighbor][0]

        return 0.0



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
        

