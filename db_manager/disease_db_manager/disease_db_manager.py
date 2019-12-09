from db_manager.db_manager import DBManager
from utils.preferences import NETWORK_MEDICINE_PATH
from utils.disease import Disease
class DiseaseManager():
    def __init__(self,table):
        self.db_table = table
        self.disease_db_path =NETWORK_MEDICINE_PATH + "db/OMIM.db"
        self.disease_db_manager = DBManager(self.disease_db_path)

    def select_disease(self, disease_id):

        tuple = (disease_id,)
        query_disease_selection = "SELECT * FROM " + self.db_table + " WHERE " + "disease_id =?"
        disease_sql = self.disease_db_manager.comupute_query(query_disease_selection,tuple)

        disease_id = disease_sql[0][0]
        disease_name = disease_sql[0][1]
        disease_genes = disease_sql[0][2].split(",")
        disease_chromosomes = disease_sql[0][3].split(",")
        disease_class = disease_sql[0][4]


        disease = Disease(disease_id,disease_name,disease_genes,disease_chromosomes,disease_class)
        return disease


    def select_all_disease(self):
        query_disease_selection = "SELECT * FROM " + self.db_table
        diseases_sql = self.disease_db_manager.comupute_query(query_disease_selection,None)

        diseases = []

        for  disease_sql in diseases_sql:

            disease_id = disease_sql[0]
            disease_name = disease_sql[1]
            disease_genes = disease_sql[2].split(",")
            disease_chromosomes = disease_sql[3].split(",")
            disease_class = disease_sql[4]

            disease = Disease(disease_id, disease_name, disease_genes, disease_chromosomes, disease_class)
            diseases.append(disease)
        return diseases



