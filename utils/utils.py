from utils.preferences import NETWORK_MEDICINE_PATH
import os



def get_number_experiments_done(experiment_id, ppi_name, disease_id):
    experiment_path = NETWORK_MEDICINE_PATH+'/exps/'+ experiment_id + '/'+ ppi_name+'/'+ disease_id + '/algo_outs'
    file_count = [len(list(filter(lambda x: '_t' not in x, files))) for r, d, files in os.walk(experiment_path)]
    return sum(file_count)

def print_experiment_stats(experiment_id,ppi_name,disease_id,len_list_params,number_fold):
    already_run = get_number_experiments_done(experiment_id,ppi_name,disease_id)
    total_number = len_list_params * number_fold

    print("Experiment Stats: ")
    print("N. of Experiments: ", total_number)
    print("Experiments already Ran: ", already_run )
    print("Experiments To Run: ", total_number - already_run )


def compute_directory_algorithm_file_names_dictionary(directory_root_path):

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