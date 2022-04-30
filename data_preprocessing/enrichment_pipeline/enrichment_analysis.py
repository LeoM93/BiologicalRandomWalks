import scipy.stats as stats


class EnrichmentAnalysis():

    def __init__(
        self,
        gene_to_ontologies,
        ontologies_to_genes,
        disease_genes):

        self.gene_to_ontologies = gene_to_ontologies
        self.ontologies_to_genes = ontologies_to_genes

        self.universe = set(self.gene_to_ontologies.keys())
        
        self.disease_genes = disease_genes.intersection(self.universe)
        
        self.num_of_seed_nodes = len(self.disease_genes)


    def group_disease_genes_by_term_id(self):
        self.disease_genes_by_term_id = {}
        
        print("Setting Up for Enrichment Analysis..............")
        print()

        for gene in self.disease_genes:

            if gene in self.gene_to_ontologies:
                term_descriptions = self.gene_to_ontologies[gene]

                for term_description in term_descriptions:
                    try:
                        self.disease_genes_by_term_id[term_description].append(gene)
                    except KeyError:
                        self.disease_genes_by_term_id[term_description] = [gene]

        if len(self.disease_genes_by_term_id) == 0:
            print("IMPOSSIBLE COMPUTE P VALUE: SET IS EMPTY")
            exit(1)


    def find_biological_processes(self,):

        print()
        print("Retrieving biological processes............")
        biological_processes = {}

        for term_id, disease_genes in self.disease_genes_by_term_id.items():
            print()
            print("Term ID: ", term_id, " Genes in seed nodes: ", disease_genes)
            gene_list_set = set(self.ontologies_to_genes[term_id])

            # instruction added
            gene_list_set = gene_list_set.intersection(self.universe)

            for gene in disease_genes:
                gene_list_set.remove(gene)

            gene_list_str = " "

            for gene in gene_list_set:
                gene_list_str += gene + " "

            print("Gene in the Universe: ",gene_list_str)
            print("N. of genes in the Universe: ",len(gene_list_set))
            print()

            biological_processes[term_id] = [set(disease_genes), gene_list_set]

        print()
        print("Biological processes retrieved..............")

        return biological_processes


    def build_matrix(self,biological_processes):

        print()
        print("Fisher Matrix building Process is starting............")

        fisher_matrix_by_process_id = {}

        for biological_process_id_1, biological_process_1 in biological_processes.items():

            fisher_matrix_by_process_id[biological_process_id_1] = [len(biological_process_1[0]), len(biological_process_1[1]),
                                                                self.num_of_seed_nodes - len(biological_process_1[0]),
                                                                len(self.universe) - self.num_of_seed_nodes - len(biological_process_1[1])]

        for k,v in fisher_matrix_by_process_id.items():
            print()
            print("Term ID: ",k)
            print("N. of seed nodes with ",k," as term ID: ",v[0])
            print("N. of nodes in the Universe with ",k," as term ID: ", v[1])
            print("N. of seed nodes with different Term IDs: ", v[2])
            print("N. of nodes in the Universe with different Term IDs : ", v[3])
            print()

        print("Built Fisher matrix by Term ID..............")
        print()

        return fisher_matrix_by_process_id


    def get_enirchment_analysis(self):

        self.group_disease_genes_by_term_id()
        print("Enrichment Analysis is Starting.............")

        biological_processes = self.find_biological_processes()
        fisher_matrix_by_biological_process_id = self.build_matrix(biological_processes)
        print("P-value computation Process is starting..............")

        p_value_by_biological_process_id = {}

        for biological_processes_id, biological_processes_matrix in fisher_matrix_by_biological_process_id.items():
            oddsratio, pvalue = stats.fisher_exact([[biological_processes_matrix[0], biological_processes_matrix[1]], [biological_processes_matrix[2], biological_processes_matrix[3]]])
            p_value_by_biological_process_id[biological_processes_id] = pvalue

        for k,v in p_value_by_biological_process_id.items():
            print()
            print("Term ID: ",k)
            print("P-value: ",v)

        print()
        print("P-value computation is finished.......")


        return p_value_by_biological_process_id
