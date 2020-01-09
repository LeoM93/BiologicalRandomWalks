from processing.p_value_correction import fdr_correction

class TermValidation():

    def __init__(self, ranked_list, disease_module_size, p_value, p_value_threshold,term_manager):

        self.ranked_list = ranked_list
        self.disease_module_size = disease_module_size
        self.p_value = p_value
        self.p_value_threshold = p_value_threshold
        self.term_manager = term_manager



    def _compute_p_value_fdr_correction(self,p_value,p_value_threshold):
        p_value_set = set()
        p_value_score = []

        for k, v in p_value.items():
            p_value_score.append(v[0])

        result, scores = fdr_correction(p_value_score,alpha=p_value_threshold)

        for p_value_id, item_1, item_2, item_3 in zip(p_value, p_value_score, scores, result):
            if item_3:
                p_value_set.add(p_value_id)

        return p_value_set


    def compute_metrics(self):

        p_value = self._compute_p_value_fdr_correction(self.p_value,self.p_value_threshold)

        precision_at_k = []
        ranked_list = self.ranked_list[:max(self.disease_module_size)]

        precision = 0
        for index, gene in enumerate(ranked_list):
            rank = index + 1

            term_gene = self.term_manager.find_gene_ontology(gene)

            if term_gene != -1:
                for term in term_gene:
                    if term in p_value:
                        precision += 1
                        break

            if rank in self.disease_module_size:
                precision_at_k.append(precision / rank )

        return precision_at_k




