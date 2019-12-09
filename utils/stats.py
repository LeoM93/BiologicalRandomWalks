from utils.preferences import NETWORK_MEDICINE_PATH
import unicodecsv as csv
class Stats():

    def __init__(self,diseases,environment_params):

        self.diseases = diseases

        self.ppi_name = environment_params["ppi_name"]
        self.experiment_id = environment_params["experiment_id"]
        self.current_experiment_path = NETWORK_MEDICINE_PATH + "exps/" + self.experiment_id + "/" + self.ppi_name + "/"





    def load_validation_output_to_dictionary(self,file_path):
        k = []
        column_id_to_exps = {}
        exps_to_recall = {}
        with open(file_path, "rb") as fp:
            csv_reader = csv.reader(fp, delimiter=",")
            for index_row, row in enumerate(csv_reader):
                if index_row == 0:
                    for index_col, item in enumerate(row):
                        if index_col == 0:
                            continue
                        column_id_to_exps[index_col] = item
                else:
                    for index_col, item in enumerate(row):
                        if index_col == 0:
                            k.append(item)
                            continue

                        try:
                            exp_name = column_id_to_exps[index_col]
                            exps_to_recall[exp_name].append(item)
                        except KeyError:
                            exp_name = column_id_to_exps[index_col]
                            exps_to_recall[exp_name] = [item]
        return k, exps_to_recall


