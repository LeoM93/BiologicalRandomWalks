import os

import numpy
import pandas as pd

from utils.preferences import NETWORK_MEDICINE_PATH
from validation.metrics_validation.metrics_validation import MetricsValidation
from file_manager.file_loader import load_train_test_file,load_algorithm_ranked_list


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



    def _compute_metrics(self, current_algorithm_output_file_path, train_list, test_list):

        if self.metrics_params['validation_type'] == "train_test_validation":

            output_algorithm = load_algorithm_ranked_list(current_algorithm_output_file_path,train_list,delimiter=",")

            metrics_validation = MetricsValidation(output_algorithm, test_list, self.metrics_params['disease_module_sizes'])
            recall_at_k = metrics_validation.compute_metrics(self.metrics_params['metrics'])


        return recall_at_k




    def concatenate_all(self):

        for disease in self.diseases:
            disease_id = disease.disease_id

            self.current_disease_path = self.current_experiment_path + str(disease_id) +"/"

            self.experiment_path = self.current_disease_path + "algo_outs/"

            self.validation_path = self.current_disease_path + "validation_outputs/"

            df = pd.DataFrame()

            for item in os.listdir(self.validation_path):

                if item == ".DS_Store":
                    continue

                current_df = pd.DataFrame.from_csv(self.validation_path + item)
                df = pd.concat([df, current_df],axis=1)


            df.to_csv(self.validation_path+ self.metrics_params['metrics'] + "_all.csv",index_label ='k')






    def validate_all(self,name = ""):
        for disease in self.diseases:
            disease_id = disease.disease_id

            self.current_disease_path = self.current_experiment_path + str(disease_id) +"/"

            self.experiment_path = self.current_disease_path + "algo_outs/"

            self.validation_path = self.current_disease_path + "validation_outputs/"

            train_directory_path = self.current_disease_path + "train/"
            test_directory_path = self.current_disease_path + "test/"

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
                #print(file_name)

                if len(train_files) != len(v):
                    print("N. OF TRAIN FILES DIFFERENT FROM N. OF OUTPUT OF CURRENT ALGORITHM")
                    continue



                element_to_count = 0

                for item in v:
                    file_path = k +"/" + item

                    algorithm_file_id = item.split("_")[-1]

                    current_train_file = train_directory_path + algorithm_file_id

                    current_test_file = test_directory_path + algorithm_file_id

                    train_set = load_train_test_file(current_train_file)
                    test_set = load_train_test_file(current_test_file)
                    try:
                        element_to_count += 1
                        recall_at_k = self._compute_metrics(file_path, train_set, test_set)
                        result = result + recall_at_k
                    except Exception as e:
                        pass

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