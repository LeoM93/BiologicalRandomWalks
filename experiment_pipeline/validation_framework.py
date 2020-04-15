import os

import numpy
import pandas as pd

from utils.preferences import NETWORK_MEDICINE_PATH
from validation.metrics_validation.metrics_validation import MetricsValidation
from validation.term_validation.term_validation import TermValidation
from validation.GWAS_validation.gwas_validation import GWASValidation
from file_manager.file_loader import load_train_test_file,load_algorithm_ranked_list,load_disease_enrichment_analysis
from file_manager.file_writer import write_data_on_disk

class ValidationFramework():

    def __init__(self,diseases,environment_params,metrics_params,):
        self.diseases = diseases

        self.ppi_name = environment_params["ppi_name"]
        self.experiment_id = environment_params["experiment_id"]

        self.metrics_params = metrics_params

        if environment_params['ppi_randomization']:
            self.current_experiment_path = NETWORK_MEDICINE_PATH + "exps/" + self.experiment_id + "/" + self.ppi_name + "_randomized/"

        else:
            self.current_experiment_path = NETWORK_MEDICINE_PATH + "exps/" + self.experiment_id + "/" + self.ppi_name + "/"



    def _compute_metrics(self, current_algorithm_output_file_path, train_list, current_disease_path,test_id):

        if self.metrics_params['validation_type'] == "train_test_validation":

            current_test_file = current_disease_path + "test/" +test_id
            test_list = load_train_test_file(current_test_file)

            output_algorithm = load_algorithm_ranked_list(current_algorithm_output_file_path,train_list,delimiter=",")

            metrics_validation = MetricsValidation(output_algorithm, test_list, self.metrics_params['disease_module_sizes'])
            recall_at_k = metrics_validation.compute_metrics(self.metrics_params['metrics'])

            return recall_at_k

        elif self.metrics_params['validation_type'] == "term_manager_validation":

            enrichment_analysis_path = current_disease_path +"enrichment_analysis/" + self.metrics_params['term_manager_name'] + "/" + self.metrics_params["filter"] + "/" + test_id
            p_value = load_disease_enrichment_analysis(enrichment_analysis_path)
            output_algorithm = load_algorithm_ranked_list(current_algorithm_output_file_path,train_list,delimiter=",")

            term_manager = self.metrics_params["term_manager"]

            term_validation = TermValidation(output_algorithm,self.metrics_params['disease_module_sizes'],p_value,self.metrics_params["p_value_threshold"],term_manager,self.metrics_params['term_counter'])

            precision =  term_validation.compute_metrics()

            return precision

        elif self.metrics_params['validation_type'] == "gwas_validation":

            gwas_manager = self.metrics_params["gwas_manager"]
            output_algorithm = load_algorithm_ranked_list(current_algorithm_output_file_path,train_list,delimiter=",")

            gwas_validation = GWASValidation(output_algorithm,self.metrics_params['disease_module_sizes'],gwas_manager)
            precision = gwas_validation.compute_metrics()

            return precision









    def validate_all(self,name = ""):
        for disease in self.diseases:
            disease_id = disease.disease_id

            self.current_disease_path = self.current_experiment_path + str(disease_id) +"/"

            self.experiment_path = self.current_disease_path + "algo_outs/"

            self.validation_path = self.current_disease_path + "validation_outputs/"

            train_directory_path = self.current_disease_path + "train/"
            test_directory_path = self.current_disease_path + "test/"

            self.racall_by_trial_id_file_path = self.validation_path

            train_files = os.listdir(train_directory_path)

            algorithm_output_dict =  self._compute_directory_algorithm_file_names_dictionary(self.experiment_path)
            results = {}
            tot_number = len(algorithm_output_dict)
            count = 0

            for k,v in algorithm_output_dict.items():
                count+=1
                if count%100 == 0:
                    print('Done: ',numpy.round(count/tot_number,2)*100,'%')

                result = numpy.zeros(len(self.metrics_params['disease_module_sizes']))

                # define the name of the output validation file
                directory_vector = k.split("/")
                file_name = ""
                name_flag = False
                for id,item in enumerate(directory_vector):
                    if name_flag == True:
                        if id != len(directory_vector) -1:
                            file_name += item +"-"
                        else:
                            file_name += item

                    if item == "algo_outs":
                        name_flag = True

                if len(train_files) != len(v):
                    print("N. OF TRAIN FILES DIFFERENT FROM N. OF OUTPUT OF CURRENT ALGORITHM")
                    continue



                element_to_count = 0


                print(k)
                for item in v:
                    file_path = k +"/" + item

                    algorithm_file_id = item.split("_")[-1]

                    current_train_file = train_directory_path + algorithm_file_id
                    train_set = load_train_test_file(current_train_file)

                    try:
                        recall_at_k = self._compute_metrics(file_path, train_set, self.current_disease_path,algorithm_file_id)
                        result = result + recall_at_k
                        element_to_count += 1
                    except Exception as e:
                        print(e)

                result = result / element_to_count
                print(result)
                results[file_name] = result

            pd.DataFrame(results, index = self.metrics_params['disease_module_sizes']).to_csv(self.validation_path +name+"_"+ self.metrics_params['metrics'] + ".csv",index_label ='k')



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