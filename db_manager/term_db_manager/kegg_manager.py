from db_manager.db_manager import DBManager
from db_manager.term_db_manager.term_manager import TermManager
from utils.preferences import NETWORK_MEDICINE_PATH

class KeggManager(TermManager):
    def __init__(self,name,table):
        TermManager.__init__(self,name)
        self.db_table = table
        self.kegg_db_path = NETWORK_MEDICINE_PATH + "db/Kegg_pathways.db"

        self.kegg_db_manager = DBManager(self.kegg_db_path)

    def _fetch_mirna_ontology(self,gene_list):
        kegg_query = "SELECT * FROM " + self.db_table
        kegg_table = self.kegg_db_manager.comupute_query(kegg_query,None)

        self.gene_to_term = {}
        self.term_to_gene = {}

        gene_set = set(gene_list)

        for record in kegg_table:
            gene_name = record[1]
            kegg_name = record[0]

            if gene_name in gene_set:

                try:
                    self.gene_to_term[gene_name].add(kegg_name)

                except KeyError:
                    self.gene_to_term[gene_name] = {kegg_name}

                try:
                    self.term_to_gene[kegg_name].add(gene_name)

                except KeyError:
                    self.term_to_gene[kegg_name] = {gene_name}


    def init_db(self,node_list):
        self._fetch_mirna_ontology(node_list)



    def close_connection(self):
        self.kegg_db_manager.close_connection()
        del self.kegg_db_manager


    def find_gene_ontology(self,gene_name):
        try:
            return self.gene_to_term[gene_name]
        except KeyError:
            return -1

    def get_gene_symbols_by_ontology_id(self,term_id):
        return self.term_to_gene[term_id]