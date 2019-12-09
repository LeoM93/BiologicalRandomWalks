from db_manager.term_db_manager.term_manager import TermManager


class RandomTermManager(TermManager):
    def __init__(self,name,):
        TermManager.__init__(self,name)


    def init_db(self,gene_list, random_gene_to_term):

        self.gene_to_term = {}
        self.term_to_gene = {}


        for gene in gene_list:
            try:
                annotations = random_gene_to_term[gene]

                for annotation in annotations:
                    try:
                        self.gene_to_term[gene].add(annotation)
                    except KeyError:
                        self.gene_to_term[gene] = {annotation}

                    try:
                        self.term_to_gene[annotation].add(gene)
                    except KeyError:
                        self.term_to_gene[annotation] = {gene}


            except KeyError:
                self.gene_to_term[gene] = -1


    def close_connection(self):
        pass


    def find_gene_ontology(self,gene_name):
        if gene_name in self.gene_to_term:
            return self.gene_to_term[gene_name]
        else:
            return -1

    def get_gene_symbols_by_ontology_id(self,term_id):
        return self.term_to_gene[term_id]





