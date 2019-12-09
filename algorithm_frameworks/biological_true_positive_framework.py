from file_manager.file_writer import write_data_on_disk
from algorithm.biological_true_positive import BiologicalTruePositive
from algorithm_frameworks.algorithm_framework import AlgorithmFramework

class BiologicalTruePositiveFramework(AlgorithmFramework):
    def __init__(self,id_number,train_set,ppi_network,output_directory_path,algorithm_params):

        AlgorithmFramework.__init__(self,id_number,train_set,ppi_network,output_directory_path)

        self.algorithm_params = algorithm_params


    def _check_parameters(self):
        print()
        print("Checking Input Parameters....")
        if self.algorithm_params["term_manager"] == None:
            print("TERM MANAGER NOT SET")
            return False
        return True


    def run(self):
        file_path = self.output_directory_path +  str(self.id_number)

        if self.check_ran_experiments(file_path) != -1:


            # for each disease compute the btp and save it
            print()
            print("BTP ID Number: ", self.id_number)
            print("Train and Test Stats:")
            print("number of train data: ", len(self.train_set))
            print()

            if self._check_parameters():

                btp_algorithm = BiologicalTruePositive(self.train_set,self.ppi_network, self.algorithm_params)
                ranked_list = btp_algorithm.run()
                write_data_on_disk(file_path, ranked_list, ["Node_Name", "Rank"], delimiter=',', write_mode="wb")

        else:
            print("EXPERIMENT N. " + self.id_number + "ALREADY RAN")
            print("DIRECTORY:", file_path)
