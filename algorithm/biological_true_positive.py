from algorithm.algorithm import Algorithm
from operator import itemgetter
from processing.p_value_correction import fdr_correction
import numpy as np

class BiologicalTruePositive(Algorithm):
    def __init__(self,source_genes,ppi_network,algorithm_params):

        Algorithm.__init__(self, source_genes, ppi_network,algorithm_params)


        self.term_manager = self.algorithm_params['term_manager']
        self.metrics= self.algorithm_params['metrics']
        self.p_value_threshold = self.algorithm_params['p_value_threshold']
        self.p_value = self.algorithm_params['p_value']
        self.fdr_correction = self.algorithm_params['fdr_correction']


        self.gene_expression_manager = self.algorithm_params["gene_expression_manager"]
        self.gene_expression_function = self.algorithm_params["gene_expression_function"]
        self.gene_expression_normalization = self.algorithm_params["gene_expression_normalization"]
        self.gene_expression_aggregation_function = self.algorithm_params["gene_expression_aggregation_function"]

    def _compute_p_value(self,):
        p_value_set = set()
        for term, value in self.p_value.items():

            if value[0] <= self.p_value_threshold:
                p_value_set.add(term)

        return p_value_set


    def _compute_p_value_fdr_correction(self):
        p_value_set = set()
        p_value_score = []

        for k, v in self.p_value.items():
            p_value_score.append(v[0])

        result, scores = fdr_correction(p_value_score,alpha=self.p_value_threshold)

        for p_value_id, item_1, item_2, item_3 in zip(self.p_value, p_value_score, scores, result):
            if item_3:
                p_value_set.add(p_value_id)

        return p_value_set


    def _intersection_normalized_by_disease(self,gene,p_value_set,gene_term_set):

        intersection_size = len(p_value_set.intersection(gene_term_set))
        normalized_intersection = intersection_size/len(p_value_set)

        return [gene,normalized_intersection]

    def _intersection(self,gene,p_value_set,gene_term_set):

        intersection_size = len(p_value_set.intersection(gene_term_set))

        return [gene,intersection_size]


    def _inclusion(self,gene,p_value_set,gene_term_set):

        intersection_size = len(p_value_set.intersection(gene_term_set))
        inclusion_size = intersection_size/len(gene_term_set)
        return [gene,inclusion_size]



    def _get_biological_score(self,gene,p_value_set):

        gene_terms = self.term_manager.find_gene_ontology(gene)
        biological_score = [gene, 0]

        if gene_terms != -1:
            gene_term_set = set(gene_terms)
            biological_score = self.get_biological_node_relevance(gene, gene_term_set, p_value_set)

        return biological_score

    def _get_gene_expression_score(self,gene):
        gene_expression_score = [gene, 1.0]

        if self.gene_expression_manager is not None:
            gene_expression = self._get_gene_expression_by_gene_name(gene)

            if gene_expression != -1:
                gene_expression_score = self._get_gene_expression_function(gene, gene_expression)

        return gene_expression_score



    def run(self):

        if self.fdr_correction == 1:
            p_value_set = self._compute_p_value_fdr_correction()
        else:
            p_value_set = self._compute_p_value()

        output_algorithm = []

        for gene in self.G.get_node_list():

            biological_score = self._get_biological_score(gene,p_value_set)
            gene_expression_score = self._get_gene_expression_score(gene)


            node_score = self._get_aggregated_node_score(biological_score, gene_expression_score)
            output_algorithm.append(node_score)

        sorted_ranked_list = sorted(output_algorithm, key=itemgetter(1), reverse=True)

        return sorted_ranked_list


    def _get_gene_expression_by_gene_name(self,gene):
        if self.gene_expression_normalization == "z_score":
            return self.gene_expression_manager.find_gene_z_score(gene)


    def _get_gene_expression_function(self,gene,gene_expression):

        if "l_" in self.gene_expression_function:
            p = float(self.gene_expression_function.split("_")[1])
            norm = 0
            for expression in gene_expression:
                norm += abs(expression)**p
            norm = norm**(1/p)
            return [gene,norm]

        if self.gene_expression_function == "max":
            norm = np.amax(gene_expression)
            if norm < 0:
                norm = 0
            return [gene,norm]


    def _get_aggregated_node_score(self,biological_score, gene_expression_score):

        gene_name = biological_score[0]

        if self.gene_expression_aggregation_function == "sum":
            return [gene_name, biological_score[1] + gene_expression_score[1]]
        else:
            return [gene_name, biological_score[1]]




    def get_biological_node_relevance(self, gene, gene_term_set, p_value_set):
        if self.metrics == "intersection":
            record = self._intersection(gene, p_value_set, gene_term_set)
        elif self.metrics == "intersection_normalized_by_disease":
            record = self._intersection_normalized_by_disease(gene, p_value_set, gene_term_set)
        return record




