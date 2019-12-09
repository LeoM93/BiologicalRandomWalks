from algorithm.DIAMOnD import DIAMOnDAlgirithm
from file_manager.file_writer import write_data_on_disk
from algorithm_frameworks.algorithm_framework import AlgorithmFramework


class DiamondFramework(AlgorithmFramework):
    def __init__(self,id_number,train_set,ppi_network,output_diamond_path,algorithm_params):

        AlgorithmFramework.__init__(self,id_number,train_set,ppi_network,output_diamond_path)
        self.algorithm_params = algorithm_params

    def run(self):

        output_file_path = self.output_directory_path + self.id_number

        if self.check_ran_experiments(output_file_path) != -1:

            # for each disease compute the random walk and save it
            print()
            print("DIAMOnD ID Number: ", self.id_number)
            print("Train and Test Stats:")
            print("number of train data: ", len(self.train_set))
            print()

            diamond = DIAMOnDAlgirithm(self.train_set,self.ppi_network,self.algorithm_params)
            ranked_list = diamond.run()
            write_data_on_disk(output_file_path ,ranked_list,["Node_Name", "Rank"],delimiter=',', write_mode="wb" )

        else:
            print("EXPERIMENT N. " + self.id_number + "ALREADY RAN")
            print("DIRECTORY:", output_file_path)




