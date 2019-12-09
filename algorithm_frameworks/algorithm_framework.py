import os

class AlgorithmFramework(object):
    def __init__(self,id_number,train_set,ppi_network,output_directory_path ):
        # set of diseases to analyze
        self.id_number = id_number
        # seed nodes
        self.train_set = train_set

        self.ppi_network = ppi_network

        self.output_directory_path = output_directory_path


    def check_ran_experiments(self, output_file_path):

        # if the path doesn't exist then create it
        if not os.path.exists(output_file_path):
            return 0
        else:
            return -1


    def run(self):
        raise "abstract method"