from db_manager.db_manager import DBManager
from db_manager.term_db_manager.term_manager import TermManager
from utils.preferences import NETWORK_MEDICINE_PATH


class DrugBankManager(TermManager):

    def __init__(self,name,table):
        TermManager.__init__(self,name)

        self.db_table = table
        self.drug_bank_db_path =NETWORK_MEDICINE_PATH + "db/drug_bank_drug_target_interactions.db"
        self.drug_bank_db_manager = DBManager(self.drug_bank_db_path)


    def init_db(self,gene_list):
        self._fetch_drug_target_interactions(gene_list)


    def close_connection(self):
        self.drug_bank_db_manager.close_connection()


    def _fetch_drug_target_interactions(self,gene_list,approval = "'approved'"):
        query_drug_selection = "SELECT * FROM " + self.db_table + " WHERE " + self.db_table + ".status = " + approval
        drug_target_interactions_sql = self.drug_bank_db_manager.comupute_query(query_drug_selection,None)
        self.gene_to_term = {}
        self.term_to_gene = {}

        for  record in drug_target_interactions_sql:

            drug_id = record[0]
            gene_name = record[2]

            if gene_name in gene_list:

                try:
                    self.gene_to_term[gene_name].add(drug_id)

                except KeyError:
                    self.gene_to_term[gene_name] = {drug_id}

                try:
                    self.term_to_gene[drug_id].add(gene_name)

                except KeyError:
                    self.term_to_gene[drug_id] = {gene_name}



    def close_connection(self):

        self.drug_bank_db_manager.close_connection()
        del self.drug_bank_db_manager


    def find_gene_ontology(self,gene_name):
        try:
            return self.gene_to_term[gene_name]
        except KeyError:
            return -1

    def get_gene_symbols_by_ontology_id(self,term_id):
        return self.term_to_gene[term_id]