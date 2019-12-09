from db_manager.db_manager import DBManager
from db_manager.term_db_manager.term_manager import TermManager
from utils.preferences import NETWORK_MEDICINE_PATH


class GeneOntologyManager(TermManager):
    def __init__(self,name,term_description = "biological_process"):
        TermManager.__init__(self,name)

        self.gene_ontology_db_path = NETWORK_MEDICINE_PATH + "db/GO_consortium_full"
        self.gene_ontology_db_manager = DBManager(self.gene_ontology_db_path)

        self.term_description = term_description



    def init_db(self,gene_list):

        self._fetch_gene_product_table()
        self._fetch_association_table()
        self._fetch_term_table()
        self._fetch_gene_ontology(gene_list)


    def _fetch_gene_product_table(self):
        gene_product_table = "SELECT * FROM gene_products"
        self.gene_product_dictionary = {}

        gene_product = self.gene_ontology_db_manager.comupute_query(gene_product_table,None)

        for record in gene_product:

            gene_name = record[1]
            gene_id = record[0]
            try:
                self.gene_product_dictionary[gene_name].append(gene_id)

            except KeyError:
                self.gene_product_dictionary[gene_name] = [gene_id]


    def _fetch_association_table(self):
        association_table_query = "SELECT * FROM associations"
        association_table = self.gene_ontology_db_manager.comupute_query(association_table_query,None)
        self.association = {}

        for record in association_table:
            association_id = record[0]
            term_id = record[1]
            gene_id = record[2]

            try:
                self.association[gene_id].add(term_id)

            except KeyError:
                self.association[gene_id] = {term_id}


    def _fetch_term_table(self):
        term_table_query = "SELECT * FROM terms"
        term_table = self.gene_ontology_db_manager.comupute_query(term_table_query, None)
        self.term = {}

        for record in term_table:
            term_description = record[1]
            if term_description == self.term_description:
                term_id = record[0]
                self.term[term_id] = 1



    def _fetch_gene_ontology(self,gene_list):

        self.gene_to_term = {}
        self.term_to_gene = {}


        for gene in gene_list:
            try:
                gene_product_ids = self.gene_product_dictionary[gene]
                for gene_product_id in gene_product_ids:

                    try:
                        term_ids = self.association[gene_product_id]

                        for term_id in term_ids:
                            if term_id in self.term:
                                try:
                                    self.gene_to_term[gene].add(term_id)
                                except KeyError:
                                    self.gene_to_term[gene] = {term_id}

                                try:
                                    self.term_to_gene[term_id].add(gene)

                                except KeyError:
                                    self.term_to_gene[term_id] = {gene}
                    except KeyError:
                        pass

            except KeyError:
                self.gene_to_term[gene] = -1


    def close_connection(self):
        self.gene_ontology_db_manager.close_connection()
        del self.gene_ontology_db_manager


    def find_gene_ontology(self,gene_name):
        if gene_name in self.gene_to_term:
            return self.gene_to_term[gene_name]
        else:
            return -1

    def get_gene_symbols_by_ontology_id(self,term_id):
        return self.term_to_gene[term_id]

    def print_stats(self,gene_list):

        print()
        print("Gene Ontology DB stats:")
        print("N. of terms: ",len(self.term))
        print("N. of genes in the ppi network: ", len(gene_list))
        print("N.of genes belonging to both the ppi and the Gene Ontology DB: ", len(self.gene_to_term))
        print()




