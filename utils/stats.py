from utils.preferences import NETWORK_MEDICINE_PATH
from file_manager.file_loader import load_train_test_file,load_algorithm_ranked_list,load_validation_data_frame
from file_manager.file_writer import write_data_on_disk

import os
import pandas as pd
class Stats():

    def __init__(self,diseases,environment_params):

        self.diseases = diseases

        self.ppi_name = environment_params["ppi_name"]
        self.experiment_id = environment_params["experiment_id"]
        self.current_experiment_path = NETWORK_MEDICINE_PATH + "exps/" + self.experiment_id + "/" + self.ppi_name + "/"



    def get_induced_sub_network_from_runs(self):
        for disease in self.diseases:

            self.current_experiment_path = NETWORK_MEDICINE_PATH + "exps/" + self.experiment_id + "/" + self.ppi_name + "/" + str(disease.disease_id) + "/"
            self.algorithm_outputs = self.current_experiment_path + "algo_outs/"
            self.train = self.current_experiment_path + "train/"
            self.stats_path = self.current_experiment_path + "stats/"

            if not os.path.exists(self.stats_path):
                os.mkdir(self.stats_path)

            self.induced_sub_network_path = self.stats_path + "induced_sub_networks/"

            if not os.path.exists(self.induced_sub_network_path):
                os.mkdir(self.induced_sub_network_path)


            algorithm_output_directory = self._compute_directory_algorithm_file_names_dictionary(self.algorithm_outputs)

            counter = 0

            for directory_father, files in algorithm_output_directory.items():

                retrieve_gene_distribution = {}

                table = []

                for file in files:

                    train_set = load_train_test_file(self.train + file)

                    algorithm_ranked_list = load_algorithm_ranked_list(directory_father + "/" +file,train_set)[:200]

                    for node in algorithm_ranked_list:
                        if node not in retrieve_gene_distribution:
                            retrieve_gene_distribution[node] = 1
                        else:
                            retrieve_gene_distribution[node] += 1

                for node,value in retrieve_gene_distribution.items():
                    if value == len(files):
                        table.append([node])

                if len(table) != 0:

                    file_name = directory_father.split("/")[17]
                    counter+= 1

                    write_data_on_disk(self.induced_sub_network_path + file_name + "_" + str(counter),table, headers=None)





    def get_algorithm_comparison_over_axis(self,algorithm_names,file_name,disease_module_size = 60):

        table = []

        for disease in self.diseases:

            self.current_experiment_path = NETWORK_MEDICINE_PATH + "exps/" + self.experiment_id + "/" + self.ppi_name + "/" + str(disease.disease_id) + "/"
            self.validation = self.current_experiment_path + "validation_outputs/"

            algorithm_validation = load_validation_data_frame(self.validation + file_name)

            for i in range(len(algorithm_names)):
                for j in range(i+1, len(algorithm_names)):

                    x = algorithm_validation[algorithm_names[i]][disease_module_size]
                    y = algorithm_validation[algorithm_names[j]][disease_module_size]

                    table.append([x,y])

        self.current_experiment_path = NETWORK_MEDICINE_PATH + "exps/" + self.experiment_id + "/"

        if not os.path.exists(self.current_experiment_path + "aggregated_stats_diseases/"):
            os.mkdir(self.current_experiment_path + "aggregated_stats_diseases/")

        metric = file_name.split("_")[-1].replace(".csv","")

        write_data_on_disk(self.current_experiment_path + "aggregated_stats_diseases/algorithm_comparison_" + "_".join(algorithm_names)+"_" + metric +".tsv",table,headers=algorithm_names,delimiter="\t")



    def get_metrics_comparison_over_axis(self,algorithm_names,file_names, disease_module_size = 200):

        table = []
        for disease in self.diseases:

            self.current_experiment_path = NETWORK_MEDICINE_PATH + "exps/" + self.experiment_id + "/" + self.ppi_name + "/" + str(
                disease.disease_id) + "/"
            self.validation = self.current_experiment_path + "validation_outputs/"

            for i in range(len(algorithm_names)):

                record = [algorithm_names[i]]

                for file_name in file_names:
                    algorithm_validation = load_validation_data_frame(self.validation + file_name)
                    x = algorithm_validation[algorithm_names[i]][disease_module_size]
                    record.append(x)

                table.append(record)

        self.current_experiment_path = NETWORK_MEDICINE_PATH + "exps/" + self.experiment_id + "/"

        if not os.path.exists(self.current_experiment_path + "aggregated_stats_diseases/"):
            os.mkdir(self.current_experiment_path + "aggregated_stats_diseases/")


        metrics = [i.split("_")[-1].replace(".csv","") for i in file_names]

        write_data_on_disk(self.current_experiment_path + "aggregated_stats_diseases/metrics_comparison_" + "_".join(metrics) +".tsv",table,headers=metrics,delimiter="\t")







    def _get_node_intersection_between_algorithm(self):
        for disease in self.diseases:

            self.current_experiment_path = NETWORK_MEDICINE_PATH + "exps/" + self.experiment_id + "/" + self.ppi_name + "/" + str(disease.disease_id) + "/"
            self.algorithm_outputs = self.current_experiment_path + "algo_outs/"
            self.train = self.current_experiment_path + "train/"
            self.stats_path = self.current_experiment_path + "stats/"

            if not os.path.exists(self.stats_path):
                os.mkdir(self.stats_path)

            self.intersection_path = self.stats_path + "intersected_proteins/"

            if not os.path.exists(self.intersection_path):
                os.mkdir(self.intersection_path)

            algorithm_output_directory = self._compute_directory_algorithm_file_names_dictionary(self.algorithm_outputs)

            dictionary = {}

            for directory_father, files in algorithm_output_directory.items():

                file_name = directory_father.split("/")[17]

                for file in files:

                    train_set = load_train_test_file(self.train + file)
                    algorithm_ranked_list = load_algorithm_ranked_list(directory_father + "/" +file,train_set)[:200]
                    dictionary[file_name] = set(algorithm_ranked_list)

            table = []
            algorithm_pairs = set()

            for algorithm_1 in dictionary.keys():
                for algorithm_2 in dictionary.keys():
                    if algorithm_1 != algorithm_2:

                        algorithm_pair_list = sorted([algorithm_1,algorithm_2])
                        algorithm_pair = (algorithm_pair_list[0],algorithm_pair_list[1])

                        if algorithm_pair not in algorithm_pairs:

                            intersection = dictionary[algorithm_pair[0]].intersection(dictionary[algorithm_pair[1]])
                            alg_1_not_in_alg_2 = dictionary[algorithm_pair[0]].difference(intersection)
                            alg_2_not_in_alg_1 = dictionary[algorithm_pair[1]].difference(intersection)

                            table.append([algorithm_pair[0],algorithm_pair[1],"intersection",",".join(intersection)])
                            table.append([algorithm_pair[0],algorithm_pair[1],"alg_1_not_in_alg_2",",".join(alg_1_not_in_alg_2)])
                            table.append([algorithm_pair[0],algorithm_pair[1],"alg_2_not_in_alg_1",",".join(alg_2_not_in_alg_1)])

                            algorithm_pairs.add(algorithm_pair)

            write_data_on_disk(self.intersection_path + "intersection_between_pairs", table, headers=None, delimiter="\t")


    def _compute_directory_algorithm_file_names_dictionary(self,directory_root_path):


        result = [os.path.join(dp, f) for dp, dn, filenames in os.walk(directory_root_path) for f in filenames]

        algorithm_output_directories = {}

        for item in result:
            item = item.replace("\\","//")
            directory_vector = item.split("/")
            file_name = directory_vector[-1]
            directory_father = "/".join(directory_vector[:len(directory_vector) - 1])

            if ".DS_Store" not in file_name and '_t' not in file_name:

                try:
                    algorithm_output_directories[directory_father].append(file_name)

                except KeyError:
                    algorithm_output_directories[directory_father] = [file_name]

        sorted_algorithm_output_directories ={}
        for k,v in algorithm_output_directories.items():
            sorted_algorithm_output_directories[k] = sorted(v)

        return sorted_algorithm_output_directories
