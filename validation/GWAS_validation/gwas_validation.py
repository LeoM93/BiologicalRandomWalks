class GWASValidation():

    def __init__(self, ranked_list, disease_module_size,gwas_manager):

        self.ranked_list = ranked_list
        self.disease_module_size = disease_module_size
        self.gwas_manager = gwas_manager


    def compute_metrics(self):

        precision_at_k = []
        ranked_list = self.ranked_list[:max(self.disease_module_size)]

        precision = 0
        for index, gene in enumerate(ranked_list):
            rank = index + 1

            if gene in self.gwas_manager.gwas_genes:
                precision +=1

            if rank in self.disease_module_size:
                precision_at_k.append(precision / rank )

        return precision_at_k






