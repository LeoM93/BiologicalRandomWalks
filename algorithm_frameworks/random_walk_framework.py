from file_manager.file_writer import write_data_on_disk
from algorithm.random_walk_with_restart import RandomWalkWithRestartAlgorithm
from algorithm_frameworks.algorithm_framework import AlgorithmFramework

class RandomWalkFramework(AlgorithmFramework):


    def __init__(self,id_number,train_set,ppi_network,output_directory_path,algorithm_params):
        AlgorithmFramework.__init__(self, id_number,train_set,ppi_network,output_directory_path)

        # restart probability taken into account
        self.algorithm_params = algorithm_params



    def run(self,):
        file_path = self.output_directory_path + str(self.id_number)

        if self.check_ran_experiments(file_path) != -1:

            # for each disease compute the random walk and save it
            print()
            print("Random Walk ID Number: ", self.id_number)
            print("Train and Test Stats:")
            print("number of train data: ", len(self.train_set))
            print()

            rwr = RandomWalkWithRestartAlgorithm(source_genes=self.train_set, graph=self.ppi_network,
                                                             algorithm_params= self.algorithm_params)
            ranked_list = rwr.run()
            write_data_on_disk(file_path,ranked_list,[ "Node_Name", "Rank"],delimiter=',', write_mode="wb" )

        else:
            print("EXPERIMENT N. " + self.id_number +"ALREADY RAN")
            print("DIRECTORY:", file_path)












