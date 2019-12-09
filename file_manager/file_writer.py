import unicodecsv as csv
import os
import json
from pathlib import Path

class WriteEnrichmentAnalysis():

    def __init__(self,enirchment_analysis_experiment_path):


        self.enirchment_analysis_experiment_path = enirchment_analysis_experiment_path


    def check_enrichment_analysis(self,train_trial):


        if not os.path.exists(self.enirchment_analysis_experiment_path + "/" + train_trial):
            return False
        else:
            return True


    def write_enrichment_analysis(self,p_value_by_term_id,train_trial):

        file_path = self.enirchment_analysis_experiment_path + "/" + train_trial

        file_path = Path(file_path)
        table_to_write = []

        for k, v in p_value_by_term_id.items():
            record = [k]
            for item in v:
                record.append(item)
            table_to_write.append(record)

        write_data_on_disk(file_path, table_to_write, headers=["Term ID", "P-value", "Odds Ratio"], delimiter="\t")



def write_json_on_disk(file_path,data):
    file_path = Path(file_path)

    with open(file_path, 'w') as fp:
        json.dump(data, fp,indent=3)



def write_row_on_disk(output_file_path,row,delimiter = ',',write_mode = 'wb'):

    with open(output_file_path, write_mode) as file:
        csv_writer = csv.writer(file, delimiter=delimiter)
        csv_writer.writerow(row)

#write csv file on disk
def write_data_on_disk(output_file_path,rows,headers = None,delimiter = ',',write_mode = 'wb'):
    output_file_path = Path(output_file_path)

    if headers is None:
        with open(output_file_path, write_mode) as file:
            csv_writer = csv.writer(file, delimiter=delimiter)
            csv_writer.writerows(rows)
    elif headers is not None:
        with open(output_file_path, write_mode) as file:
            csv_writer = csv.writer(file, delimiter=delimiter)
            csv_writer.writerow(headers)
            csv_writer.writerows(rows)
