import scipy.stats as stats


class EnrichmentAnalysis():
    def __init__(self,db_manager,disease, universe, enrichment_case = 2,min_num_of_genes_for_term_id = 0,max_number_of_gene_for_term_id = 100000):
        self.db_manager = db_manager
        self.disease_genes = disease.disease_genes
        self.universe = universe
        self.num_of_seed_nodes = len(disease.disease_genes)

        # define the max number of genes that a term can be linked with
        self.max_number_of_gene_for_term_id = max_number_of_gene_for_term_id

        # define the min number of genes that a term can be linked with
        self.min_number_of_gene_for_term_id =min_num_of_genes_for_term_id

        # choose between two different ways to compute the fisher test
        self.enrichment_case = enrichment_case


    def group_disease_genes_by_term_id(self,seed_set):
        self.disease_genes_by_term_id = {}
        print("Setting Up for Enrichment Analysis..............")
        print()

        for gene in seed_set:

            term_descriptions = self.db_manager.find_gene_ontology(gene)
            if term_descriptions != -1:
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
            gene_list_set = set(self.db_manager.get_gene_symbols_by_ontology_id(term_id))

            for gene in disease_genes:

                gene_list_set.remove(gene)

            if (len(disease_genes) + len(gene_list_set) > self.max_number_of_gene_for_term_id) or (len(disease_genes) + len(gene_list_set)<self.min_number_of_gene_for_term_id):
                continue

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

        if self.enrichment_case == 1:

            #first case in which we consider the Universe as all the nodes with at least one term ID of the seed nodes
            for biological_process_id_1, biological_process_1 in biological_processes.items():
                #set containing nodes with different biological process
                other_biological_processes = set()

                for biological_process_id_2, biological_process_2 in biological_processes.items():

                    #if the biological process is different from that one taken in consideration
                    if biological_process_id_1 != biological_process_id_2:
                        #do the union of sets
                        other_biological_processes = other_biological_processes.union(biological_process_2[1])

                other_biological_processes_set = set()
                for gene in other_biological_processes:
                    #check if a gene is in the seed nodes or in the biological process taken in consideration
                    if gene in self.disease_genes or gene in biological_process_1[1]:
                        continue
                    other_biological_processes_set.add(gene)
            #end of the first case

                #create the fisher matrix
                fisher_matrix_by_process_id[biological_process_id_1] = [len(biological_process_1[0]),len(biological_process_1[1]),self.num_of_seed_nodes - len(biological_process_1[0]),len(other_biological_processes_set)]

        else:
            # second case in which we consider the Universe as all the nodes of the ppi_network.
            for biological_process_id_1, biological_process_1 in biological_processes.items():

                #second case for computing the p_value
                fisher_matrix_by_process_id[biological_process_id_1] = [len(biological_process_1[0]),
                                                                len(biological_process_1[1]),
                                                                self.num_of_seed_nodes - len(biological_process_1[0]),
                                                                len(self.universe.get_node_list()) - self.num_of_seed_nodes - len(biological_process_1[1])]

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


    def get_enirchment_analysis(self,seed_set):


        self.group_disease_genes_by_term_id(seed_set)

        print("Enrichment Analysis is Starting.............")

        biological_processes = self.find_biological_processes()

        fisher_matrix_by_biological_process_id = self.build_matrix(biological_processes)

        print("P-value computation Process is starting..............")

        p_value_by_biological_process_id = {}

        for biological_processes_id, biological_processes_matrix in fisher_matrix_by_biological_process_id.items():
            oddsratio, pvalue = stats.fisher_exact([[biological_processes_matrix[0], biological_processes_matrix[1]], [biological_processes_matrix[2], biological_processes_matrix[3]]])
            p_value_by_biological_process_id[biological_processes_id] = [pvalue,oddsratio]

        for k,v in p_value_by_biological_process_id.items():
            print()
            print("Term ID: ",k)
            print("Odds ratio: ",v[1])
            print("P-value: ",v[0])
        print()
        print("P-value computation is finished.......")


        return p_value_by_biological_process_id