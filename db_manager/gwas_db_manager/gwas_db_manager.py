from db_manager.db_manager import DBManager
from utils.preferences import NETWORK_MEDICINE_PATH


class GWASManager():
    def __init__(self,table):
        self.db_table = table
        self.gwas_db_path =NETWORK_MEDICINE_PATH + "db/GWAS_Collection.db"
        self.gwas_db_manager = DBManager(self.gwas_db_path)


    def init_db(self):
        self.select_gwas_genes()


    def close_connection(self):
        self.gwas_db_manager.close_connection()


    def select_gwas_genes(self):
        query_gwas_selection = "SELECT * FROM " + self.db_table
        genes_sql = self.gwas_db_manager.comupute_query(query_gwas_selection,None)

        self.gwas_genes = {}

        for  gene in genes_sql:

            gene_name = gene[0]
            gene_p_value = float(gene[1])

            self.gwas_genes[gene_name] = gene_p_value



