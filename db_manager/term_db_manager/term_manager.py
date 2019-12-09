class TermManager():
    def __init__(self,name):
        self.name = name


    def close_connection(self):
        raise "abstract mwthod"


    def find_gene_ontology(self,gene_name):
        raise "abstract method"

    def get_gene_symbols_by_ontology_id(self,term_id):
        raise "abstract method"
