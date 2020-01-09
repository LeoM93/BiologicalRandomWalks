import os
from multiprocessing import Pool

from utils.preferences import NETWORK_MEDICINE_PATH
import itertools

from algorithm_frameworks.random_walk_framework import RandomWalkFramework
from algorithm_frameworks.biological_true_positive_framework import BiologicalTruePositiveFramework
from algorithm_frameworks.diamond_framework import DiamondFramework
from algorithm_frameworks.biological_random_walk_framework import BiologicalRandomWalkFramework

from file_manager.file_loader import load_train_test_file,load_disease_enrichment_analysis
from file_manager.file_loader import LoadAlgorithmVariables

from utils.utils import print_experiment_stats

class NetworkAlgorithmFramework():

    def __init__(self,environment_params,enriched_analysis_params,diseases,ppi_network,db_managers,gene_expression_manager):

        self.experiment_id = environment_params['experiment_id']
        self.ppi_name = environment_params['ppi_name']
        self.disease_list = diseases
        self.ppi_network = ppi_network
        self.number_fold = environment_params['num_of_fold']
        self.enrichment_analysis_params = enriched_analysis_params
        self.db_managers = db_managers
        self.gene_expression_manager = gene_expression_manager

        if environment_params['ppi_randomization']:
            self.experiment_path = NETWORK_MEDICINE_PATH + "exps/" + self.experiment_id + "/" + environment_params['ppi_name'] + "_randomized/"
        else:
            self.experiment_path =  NETWORK_MEDICINE_PATH + "exps/" + self.experiment_id + "/" + environment_params['ppi_name'] + "/"


    def _run_chosen_algorithm(self,output_directory_path,seed_set_path,experiment_num,disease_id):
        if self.algorithm_params['algorithm'] == 'rwr':
            self._run_random_walk_with_restart(output_directory_path,seed_set_path,experiment_num)
        elif self.algorithm_params['algorithm'] == 'btp':
            self._run_biological_true_positive(output_directory_path,seed_set_path,experiment_num,disease_id)
        elif self.algorithm_params['algorithm'] == 'DIAMOnD':
            self._run_DIAMOnD(output_directory_path,seed_set_path,experiment_num)
        elif self.algorithm_params['algorithm'] == 'BRW':
            self._run_BRW(output_directory_path,seed_set_path,experiment_num,disease_id)



    def _run_random_walk_with_restart(self,output_file_path,seed_set_path,experiment_number):

        seed_set = load_train_test_file(seed_set_path)
        rwr_framework = RandomWalkFramework(experiment_number, seed_set,
                                            self.ppi_network, output_file_path,self.algorithm_params)

        rwr_framework.run()

    def _run_biological_true_positive(self,output_file_path,seed_set_path,experiment_num,disease_id):

        self._compute_enrichment_analysis_path(str(disease_id), experiment_num,self.algorithm_params['term_manager'].name)
        p_value = load_disease_enrichment_analysis(self.enrichment_analysis_path)
        self.algorithm_params['p_value'] = p_value

        seed_set = load_train_test_file(seed_set_path)
        btp_framework = BiologicalTruePositiveFramework(experiment_num,seed_set,self.ppi_network,output_file_path,self.algorithm_params)
        btp_framework.run()


    def _run_BRW(self,output_file_path,seed_set_path,experiment_num,disease_id):
        seed_set = load_train_test_file(seed_set_path)

        if self.algorithm_params["teleporting_parameters"]["name"] == "biological_teleporting" and self.algorithm_params["walking_parameters"]["name"] == "default":

            self._compute_enrichment_analysis_path(str(disease_id), experiment_num, self.algorithm_params["teleporting_parameters"]['term_manager'].name)
            p_value = load_disease_enrichment_analysis(self.enrichment_analysis_path)
            self.algorithm_params["teleporting_parameters"]['p_value'] = p_value



        if self.algorithm_params["teleporting_parameters"]["name"] == "biological_teleporting" and self.algorithm_params["walking_parameters"]["name"] == "biological_walking":

            self._compute_enrichment_analysis_path(str(disease_id), experiment_num,
                                                   self.algorithm_params["teleporting_parameters"]['term_manager'].name)
            teleporting_p_value = load_disease_enrichment_analysis(self.enrichment_analysis_path)
            self.algorithm_params["teleporting_parameters"]['p_value'] = teleporting_p_value


            self._compute_enrichment_analysis_path(str(disease_id), experiment_num,
                                                   self.algorithm_params["walking_parameters"]['term_manager'].name)
            walking_p_value = load_disease_enrichment_analysis(self.enrichment_analysis_path)
            self.algorithm_params["walking_parameters"]['p_value'] = walking_p_value

        brw_framework = BiologicalRandomWalkFramework(experiment_num,seed_set,self.ppi_network,output_file_path,self.algorithm_params)
        brw_framework.run()


    def _run_DIAMOnD(self,output_file_path,seed_set_path,experiment_num):

        seed_set = load_train_test_file(seed_set_path)
        diamond_framework = DiamondFramework(experiment_num,seed_set,self.ppi_network,output_file_path,self.algorithm_params)
        diamond_framework.run()


    def _compute_enrichment_analysis_path(self,disease_id,experiment_num,db_name):
        if self.enrichment_analysis_params['filter_flag']:
            self.enrichment_analysis_path = self.experiment_path + disease_id + "/enrichment_analysis/" + \
                                            self.enrichment_analysis_params["db_name"] + "/"+\
                                            str(self.enrichment_analysis_params['min_num_of_genes_for_term_id']) + "_" + \
                                            str(self.enrichment_analysis_params['max_num_of_genes_for_term_id']) + "/" + experiment_num

        else:

            self.enrichment_analysis_path = self.experiment_path + disease_id + "/enrichment_analysis/" + \
                                            db_name + "/no_filter/" + experiment_num


    def run(self,algorithm_params):

        self.algorithm_params = algorithm_params


        for disease in self.disease_list:

            disease_id = disease.disease_id

            seed_path = self.experiment_path + str(disease_id) + "/train/"

            seed_files = os.listdir(seed_path)

            output_directory_path = self.set_up_chosen_algorithm(disease_id)

            for experiment_num in seed_files:

                current_seed_file_path = seed_path + experiment_num

                self._run_chosen_algorithm(output_directory_path,current_seed_file_path,experiment_num,disease_id)


    def _get_one_combination_of_brw_nested_dictionary(self,dictionary):

        dictionary_list = []

        list_parameter_name = []
        list_parameter_value_values = []

        for par, values in dictionary.items():
            list_parameter_name.append(par)
            list_parameter_value_values.append(values)

        for product in list(itertools.product(*list_parameter_value_values)):
            new_dictionary = {}
            for i, parameter_name in enumerate(list_parameter_name):
                new_dictionary[parameter_name] =product[i]
            dictionary_list.append(new_dictionary)
        return dictionary_list



    def get_all_combination_of_brw_parameters(self,fdr_corrections,
                                          teleporting_parameter_combination,
                                          teleporting_score_function_combination, teleporting_score_function_name,
                                          walking_combinations):

        dictionary_list = []
        algorithm = "BRW"


        for fdr_correction in fdr_corrections:
            for teleporting in teleporting_parameter_combination:
                for scoring_function in teleporting_score_function_combination:
                    for walking in walking_combinations:

                        if teleporting['gene_expression_coefficient'] + teleporting['term_manager_coefficient'] != 1.0:
                            continue

                        if teleporting["name"] == "biological_teleporting" and walking["name"] == "default":

                            teleporting['score_function'] = {
                                'name': teleporting_score_function_name,
                                'parameters': scoring_function
                            }

                            if teleporting['term_manager'] == "go":
                                teleporting['term_manager'] = self.db_managers["go"]

                            elif teleporting['term_manager'] == "mirna":
                                teleporting['term_manager'] = self.db_managers["mirna"]

                            if teleporting["gene_expression"] == "TCGA":
                                teleporting["gene_expression"] = self.gene_expression_manager

                            dictionary = {'teleporting_parameters': teleporting, 'fdr_correction': fdr_correction,
                                          'algorithm': algorithm, 'walking_parameters': {'name': walking["name"]}}
                            dictionary['walking_parameters']['edge_relevance'] = -1

                            if fdr_correction == 0:
                                dictionary['fdr_correction'] = "false"
                            else:
                                dictionary['fdr_correction'] = "true"



                            dictionary_list.append(dictionary)


                        else:

                            teleporting['score_function'] = {
                                'name': teleporting_score_function_name,
                                'parameters': scoring_function
                            }

                            if teleporting['term_manager'] == "go":
                                teleporting['term_manager'] = self.db_managers["go"]

                            elif teleporting['term_manager'] == "mirna":
                                teleporting['term_manager'] = self.db_managers["mirna"]

                            if teleporting["gene_expression"] == "TCGA":
                                teleporting["gene_expression"] = self.gene_expression_manager

                            if walking['term_manager'] == "go":
                                walking['term_manager'] = self.db_managers["go"]

                            elif walking['term_manager'] == "mirna":
                                walking['term_manager'] = self.db_managers["mirna"]



                            dictionary = {'teleporting_parameters': teleporting,
                                          'walking_parameters': walking,
                                          'fdr_correction': fdr_correction,
                                          'algorithm': algorithm

                                          }
                            if fdr_correction == 0:
                                dictionary['fdr_correction'] = "false"
                            else:
                                dictionary['fdr_correction'] = "true"

                            dictionary_list.append(dictionary)
        return dictionary_list

    def get_all_combinations_of_brw(self,brw):


        teleporting_parameters = brw['teleporting_parameters']
        teleporting_score_function = teleporting_parameters['score_function']['parameters']
        teleporting_score_function_name = teleporting_parameters['score_function']['name']
        walking_parameters = brw['walking_parameters']
        teleporting_parameters.pop('score_function', None)

        fdr_corrections = brw["fdr_correction"]

        teleporting_parameter_combination = self._get_one_combination_of_brw_nested_dictionary(teleporting_parameters)
        teleporting_score_function_combination = self._get_one_combination_of_brw_nested_dictionary(teleporting_score_function)
        walking_combinations = self._get_one_combination_of_brw_nested_dictionary(walking_parameters)

        dictionary_list = self.get_all_combination_of_brw_parameters(fdr_corrections, teleporting_parameter_combination,
                                                                 teleporting_score_function_combination,
                                                                 teleporting_score_function_name,
                                                                 walking_combinations)
        return dictionary_list


    def run_multiprocess(self, num_pool = 3):

        load_algorithm_params = LoadAlgorithmVariables()

        algorithm_params_dict = load_algorithm_params.load_algorithm_params( NETWORK_MEDICINE_PATH + "exps/" + self.experiment_id + "/experiment_variables/")

        algorithm_params_list = []

        for algorithm_name,algorithm_list in algorithm_params_dict.items():
            if algorithm_name == "BRW":
                for item in algorithm_list:
                    algorithm_params_list += self.get_all_combinations_of_brw(item)

            elif algorithm_name == "btp":

                for btp_dict in algorithm_list:
                    btp_algorithm = self.get_btw_algorithm(btp_dict)
                    algorithm_params_list.append(btp_algorithm)

            else:
                algorithm_params_list += algorithm_list


        for item in algorithm_params_list:
            print()
            print(item)
            print()


        processing = Pool(num_pool)
        processing.map(self.run, algorithm_params_list)
        processing.close()
        processing.join()

    def get_btw_algorithm(self, btp_dict):
        algorithm_params = btp_dict.copy()

        if btp_dict["gene_expression_manager"] == "TCGA":
            algorithm_params["gene_expression_manager"] = self.gene_expression_manager
        else:
            algorithm_params["gene_expression_manager"] = None

        if btp_dict["term_manager"] == "go":
            algorithm_params["term_manager"] = self.db_managers["go"]
        elif btp_dict["term_manager"] == "mirna":
            algorithm_params["term_manager"] = self.db_managers["mirna"]
        elif btp_dict["term_manager"] == "random":
            algorithm_params["term_manager"] = self.db_managers["random"]
        else:
            algorithm_params["term_manager"] = None

        return algorithm_params

    def set_up_chosen_algorithm(self,disease_id):

        if self.algorithm_params['algorithm'] == 'rwr':
            output_file_path = self._set_up_rwr(disease_id)
        elif self.algorithm_params['algorithm'] == 'btp':
            output_file_path = self._set_up_btp(disease_id)
        elif self.algorithm_params['algorithm'] == 'DIAMOnD':
            output_file_path = self._set_up_DIAMOnD(disease_id)
        elif self.algorithm_params['algorithm'] == 'BRW':
            output_file_path = self._set_up_brw(disease_id)

        return output_file_path



    def _set_up_brw(self, disease_id):

        t_par = self.algorithm_params['teleporting_parameters']


        term_manager_teleporting_path = "teleporting_term_manager_"+ self.algorithm_params['teleporting_parameters']['term_manager'].name + "/fdr_" + \
                           self.algorithm_params['fdr_correction'] + "/teleporting_p_val_" + \
                           str(self.algorithm_params['teleporting_parameters']['p_value_threshold'] * 100)[
                           :3] + "/rel_" + t_par['node_relevance'] + '/f_' + str(
            t_par['score_function']['name']) + "/" + \
                           '/'.join([key + '_' + str(value) if not isinstance(value, str) else value for key, value in
                                     t_par['score_function']['parameters'].items()]) + "/t_" + str(
            t_par['clip_t']) +"/teleporting_term_manager_coefficient_" + str(t_par["term_manager_coefficient"])




        gene_expression_teleporting_path = "/teleporting_ge_" + t_par["gene_expression"].name +"/teleporting_ge_normalization_"+t_par["gene_expression_normalization"]+ \
                                    "/teleporting_ge_f_" + t_par["gene_expression_function"] +"/teleporting_ge_coefficient_"+ str(t_par["gene_expression_coefficient"])+\
                                                   "/teleporting_aggregation_function_"+ t_par["teleporting_aggregation_function"]+"/r_" + str(
                    t_par['restart_probability'] * 100)[:3] + "/"



        try:


            walking_path = "walking_term_manager_" +self.algorithm_params['walking_parameters']['term_manager'].name + "/walking_p_val_" + str(
            self.algorithm_params['walking_parameters']['p_value_threshold'] * 100)[:3] \
                       + "/edge_" + str(self.algorithm_params['walking_parameters']['edge_relevance']) \
                       + "/edge_score_" + str(self.algorithm_params['walking_parameters']['edge_score']) + '/walking_aggregation_function_' + self.algorithm_params['walking_parameters']['walking_aggregation_function'] + "/"

            if  self.algorithm_params['walking_parameters']['walking_aggregation_function'] != "None":
                walking_path += "pearson_threshold_" + str(self.algorithm_params['walking_parameters']['threshold_pearson_correlation']) + "/"


        except KeyError:
            walking_path = ""

        if self.algorithm_params['teleporting_parameters']['name'] == 'biological_teleporting' and self.algorithm_params['walking_parameters']['name'] == 'default':

            starting_path = self.experiment_path + str(disease_id) + "/algo_outs/" + "BRW/" + "t_b_w_d/"

            output_file_path = starting_path + term_manager_teleporting_path + gene_expression_teleporting_path

        else:
            starting_path = self.experiment_path + str(disease_id) + "/algo_outs/" + "BRW/" + "t_b_w_b/"

            output_file_path = starting_path + term_manager_teleporting_path +  gene_expression_teleporting_path + walking_path



        if not os.path.exists(output_file_path):
            os.makedirs(output_file_path)
        return output_file_path

    def _set_up_DIAMOnD(self, disease_id):
        output_file_path = self.experiment_path + str(disease_id) + "/algo_outs/" + "DIAMOnD/"
        if not os.path.exists(output_file_path):
            os.makedirs(output_file_path)
        return output_file_path



    def _set_up_rwr(self, disease_id):
        try:
            restart_probability = self.algorithm_params['restart_probability']
        except KeyError:
            print("RESTART PROBABILITY NOT SET")
            exit(1)

        output_file_path = self.experiment_path + str(disease_id) + "/algo_outs/" + "rwr/" + str(
            restart_probability * 100) + "/"

        if not os.path.exists(output_file_path):
            os.makedirs(output_file_path)
        return output_file_path



    def _set_up_btp(self, disease_id):

        try:
            metrics = self.algorithm_params['metrics']
        except KeyError:

            print("ALGORITHM OPTIONS ARE NOT DEFINED")
            exit(1)

        output_file_path = ""

        term_manager_file_path =  self.experiment_path + str(disease_id) + "/algo_outs/" + "btp/" + \
                           self.algorithm_params['term_manager'].name+ "/p_val_" + \
                               str(self.algorithm_params['p_value_threshold']*100)[:3] + "/" + self.algorithm_params[
                               'metrics']

        output_file_path += term_manager_file_path

        if self.algorithm_params["gene_expression_manager"] is not None:
            gene_expression_file_path = "/gene_expression_manager_" + self.algorithm_params["gene_expression_manager"].name + "/" \
                           + "gene_expression_function_" + self.algorithm_params["gene_expression_function"] + "/gene_expression_aggregation_function_" \
                                        + self.algorithm_params["gene_expression_aggregation_function"] + "/"

            output_file_path += gene_expression_file_path
        else:
            output_file_path += "/gene_expression_manager_none/"


        if not os.path.exists(output_file_path):
            os.makedirs(output_file_path)
        return output_file_path