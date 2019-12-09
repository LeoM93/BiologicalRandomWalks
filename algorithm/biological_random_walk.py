import numpy as np

from algorithm.algorithm import Algorithm
from processing.p_value_correction import fdr_correction
import math

CONV_THRESHOLD = 0.000001

class BiologicalRandomWalk(Algorithm):

    def __init__(self,seed_set,ppi_network,algorithm_params):
        Algorithm.__init__(self,seed_set,ppi_network,algorithm_params)

        self._init_biological_random_walk()



    def _init_biological_random_walk(self):


        if self.algorithm_params['teleporting_parameters']["name"] == "biological_teleporting" and self.algorithm_params['walking_parameters']["name"] == "default":

            # term manager field
            self.teleporting_p_value = self.algorithm_params["teleporting_parameters"]['p_value']
            self.restart_probability = self.algorithm_params['teleporting_parameters']['restart_probability']
            self.t = self.algorithm_params['teleporting_parameters']['clip_t']
            self.relevance_function = self.algorithm_params['teleporting_parameters']['node_relevance']
            self.score_function = self.algorithm_params['teleporting_parameters']['score_function']
            self.teleporting_p_value_threshold = self.algorithm_params['teleporting_parameters']['p_value_threshold']
            self.fdr_correction = self.algorithm_params['fdr_correction']
            self.teleporting_term_manager = self.algorithm_params["teleporting_parameters"]['term_manager']
            self.teleporting_term_manager_coefficient = self.algorithm_params["teleporting_parameters"][
                'term_manager_coefficient']

            # gene expression field
            self.gene_expression_manager = self.algorithm_params["teleporting_parameters"]["gene_expression"]
            self.gene_expression_normalization = self.algorithm_params["teleporting_parameters"][
                "gene_expression_normalization"]
            self.gene_expression_function = self.algorithm_params["teleporting_parameters"]["gene_expression_function"]
            self.gene_expression_coefficient = self.algorithm_params["teleporting_parameters"][
                "gene_expression_coefficient"]
            self.teleporting_aggregation_function = self.algorithm_params["teleporting_parameters"][
                "teleporting_aggregation_function"]


        elif self.algorithm_params['teleporting_parameters']["name"] == "biological_teleporting" and self.algorithm_params['walking_parameters']["name"] == "biological_walking":

            # term manager field
            self.teleporting_p_value = self.algorithm_params["teleporting_parameters"]['p_value']
            self.restart_probability = self.algorithm_params['teleporting_parameters']['restart_probability']
            self.t = self.algorithm_params['teleporting_parameters']['clip_t']
            self.relevance_function = self.algorithm_params['teleporting_parameters']['node_relevance']
            self.score_function = self.algorithm_params['teleporting_parameters']['score_function']
            self.teleporting_p_value_threshold = self.algorithm_params['teleporting_parameters']['p_value_threshold']
            self.fdr_correction = self.algorithm_params['fdr_correction']
            self.teleporting_term_manager = self.algorithm_params["teleporting_parameters"]['term_manager']
            self.teleporting_term_manager_coefficient = self.algorithm_params["teleporting_parameters"]['term_manager_coefficient']


            # gene expression field
            self.gene_expression_manager = self.algorithm_params["teleporting_parameters"]["gene_expression"]
            self.gene_expression_normalization = self.algorithm_params["teleporting_parameters"]["gene_expression_normalization"]
            self.gene_expression_function = self.algorithm_params["teleporting_parameters"]["gene_expression_function"]
            self.gene_expression_coefficient = self.algorithm_params["teleporting_parameters"]["gene_expression_coefficient"]
            self.teleporting_aggregation_function = self.algorithm_params["teleporting_parameters"]["teleporting_aggregation_function"]




            self.walking_p_value = self.algorithm_params["walking_parameters"]['p_value']
            self.walking_term_managet = self.algorithm_params["walking_parameters"]['term_manager']
            self.walking_p_value_threshold = self.algorithm_params['walking_parameters']['p_value_threshold']



            self.gene_expression_biological_walking_aggregation = self.algorithm_params['walking_parameters']['walking_aggregation_function']


    def _compute_p_value(self,p_value,p_value_threshold):
        p_value_set = set()
        for term, value in p_value.items():

            if value[0] <= p_value_threshold:
                p_value_set.add(term)

        return p_value_set

    def _compute_p_value_rank(self,):


        dict_pval = {k:v[0] for k,v in self.teleporting_p_value.items()}
        sorted_by_value = sorted(dict_pval.keys(), key=dict_pval.get)

        p_value_dict = {k:self._inverse_sigmoid(i) for i,k in enumerate(sorted_by_value)}

        self.sum_p_value_rank = sum(p_value_dict.values())
        return p_value_dict

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


    def _initialize_page_rank(self):

        if self.algorithm_params['teleporting_parameters']['name'] == 'default':
            return self._default()

        else:
            return self._get_aggregated_initial_page_rank()


    def _get_aggregated_initial_page_rank(self):

        biological_teleporting = self._biological_teleporting()
        gene_expression_teleporting = self._get_gene_expression_teleporing_probability()

        return self._get_aggregation_function(biological_teleporting,gene_expression_teleporting)



    def _get_aggregation_function(self,pr_1, pr_2):

        aggregated_teleporting = {}

        if self.teleporting_aggregation_function == "sum":

            aggregated_teleporting = {k: self.teleporting_term_manager_coefficient * pr_1[k] +
                                         self.gene_expression_coefficient * pr_2[k]
                                      for k in pr_1.keys()}

        elif self.teleporting_aggregation_function == "product":
            aggregated_teleporting = {k: self.teleporting_term_manager_coefficient * pr_1[k] *
                                         self.gene_expression_coefficient * pr_2[k]
                                      for k in pr_1.keys()}


        total = sum(aggregated_teleporting.values())
        page_rank_dictionary = {k: v / total for k, v in aggregated_teleporting.items()}

        self.initial_pr = page_rank_dictionary
        return page_rank_dictionary



    def _set_graph_weights(self):
        if self.algorithm_params['walking_parameters']['name'] == 'default':
            self._compute_graph_weight()

        else:

            if self.fdr_correction == 'false':
                self.enriched_walking_set = self._compute_p_value(self.walking_p_value, self.walking_p_value_threshold)

            else:
                self.enriched_walking_set = self._compute_p_value_fdr_correction(self.walking_p_value,
                                                                                 self.walking_p_value_threshold)

            if len(self.enriched_walking_set) == 0:
                self.enriched_walking_set = {}

            self._compute_graph_weight()




    def _compute_graph_weight(self):

        self.gene_to_gene_weight = {}

        for gene in self.G.ppi_network:
            self.gene_to_gene_weight[gene] = {}
            neighbors = self.G.get_neighbors(gene)

            total_weight = 0

            for neighbor in neighbors:

                if self.algorithm_params['walking_parameters']['edge_relevance'] == -1:
                    self.gene_to_gene_weight[gene][neighbor] = 1
                    total_weight+=1
                    continue

                gene_terms = self.walking_term_managet.find_gene_ontology(gene)
                neighbor_terms = self.walking_term_managet.find_gene_ontology(neighbor)

                if gene_terms == -1 or neighbor_terms == -1:
                    self.gene_to_gene_weight[gene][neighbor] = 1
                    total_weight += 1
                    continue

                intersection = gene_terms.intersection(neighbor_terms)
                union = gene_terms.union(neighbor_terms)


                if self.algorithm_params['walking_parameters']['edge_relevance'] == 0:
                    intersection = len(intersection)
                    union = len(union)

                elif self.algorithm_params['walking_parameters']['edge_relevance'] == 1:
                    intersection = len(intersection.intersection(self.enriched_walking_set))
                    union = len(union)


                elif self.algorithm_params['walking_parameters']['edge_relevance'] == 2:

                    intersection = len(intersection.intersection(self.enriched_walking_set))
                    union = len(union.intersection(self.enriched_walking_set))

                else:

                    print()
                    print(" Edge Relevance: " + str(self.algorithm_params['walking_parameters']['edge_relevance']) +
                          " is not a admissible option")
                    exit(3)

                if union == 0:
                    score = 0
                else:
                    score = intersection/union


                if self.algorithm_params['walking_parameters']['edge_score'] == 'sum':

                    score = 1 + score
                else:

                    print()
                    print(" Edge Score: " + self.algorithm_params['walking_parameters']['edge_scpore'] +
                          " is not a admissible option")
                    exit(3)


                self.gene_to_gene_weight[gene][neighbor] = self._compute_graph_aggregated_weight(gene,neighbor,score)
                total_weight += self.gene_to_gene_weight[gene][neighbor]

            for neighbor in neighbors:

                self.gene_to_gene_weight[gene][neighbor] = self.gene_to_gene_weight[gene][neighbor]/total_weight


    def _compute_graph_aggregated_weight(self,gene,neighbor,biological_score):

        if self.gene_expression_biological_walking_aggregation == "sum":

            final_score = biological_score + self.gene_expression_manager.find_gene_pearson_correlation_by_ppi_network(gene,neighbor)
            return final_score

        else:
            return biological_score


    def _default(self):
        page_rank_dictionary = {}


        for gene in self.G.get_node_list():
            if gene in self.source_node_list:
                page_rank_dictionary[gene] = 1/len(self.source_node_list)
            else:
                page_rank_dictionary[gene] = 0
        self.initial_pr = page_rank_dictionary
        return page_rank_dictionary


    def _sigmoid(self,x,steep,translation):
        return 1 / (1 + math.exp(-steep*(x-translation)))



    def _get_node_relevance(self,gene,relevance_function):

        if relevance_function == 'inclusion':
            if gene in self.source_node_list:
                return 1
            else:
                gene_terms = self.teleporting_term_manager.find_gene_ontology(gene)
                if gene_terms == -1:
                    return 0
                else:
                    return len(self.enriched_set.intersection(gene_terms))/len(gene_terms)

        elif relevance_function == 'intersection':
            if gene in self.source_node_list:
                return 1

            else:
                gene_terms = self.teleporting_term_manager.find_gene_ontology(gene)
                if gene_terms == -1:
                    return 0
                else:
                    return len(self.enriched_set.intersection(gene_terms))/len(self.enriched_set)

        elif self.relevance_function == 'product':
            return (1+self._get_node_relevance(gene,'inclusion'))*(1+self._get_node_relevance(gene,'intersection'))/4

        else:
            print("NO NODE RELEVANCE FUNCTION FOUND")
            exit(10)



    def _get_node_score(self,relevance):


        if self.score_function['name'] == 'default':
            return  relevance

        elif self.score_function['name'] == 'power':
            return relevance**self.score_function['parameters']['power']

        elif self.score_function['name'] == 'sigmoid':
            return self._sigmoid(relevance,self.score_function['parameters']['steep'],self.score_function['parameters']['translation'])

        elif self.score_function['name'] == 'ranking':
            return relevance


    def _inverse_sigmoid(self,rank):
        steep = self.score_function['parameters']['steep']
        translation = self.score_function['parameters']['translation']
        return 1 - (1/(1+math.exp(-steep*(rank-translation))))


    def _biological_teleporting(self):

        if self.fdr_correction == 'false':
            self.enriched_set = self._compute_p_value(self.teleporting_p_value,self.teleporting_p_value_threshold)

        else:
            self.enriched_set = self._compute_p_value_fdr_correction(self.teleporting_p_value,self.teleporting_p_value_threshold)


        if len(self.enriched_set) == 0:
            return {}

        if self.score_function['name'] == 'default':
            page_rank_dictionary = {gene:  min(self.t,self._get_node_relevance(gene,self.relevance_function)) for gene in self.G.get_node_list() }

        elif self.score_function['name'] != 'ranking':
            page_rank_dictionary = {gene:  min(self.t,self._get_node_score(self._get_node_relevance(gene,self.relevance_function))) for gene in self.G.get_node_list() }


        else:

            sorted_pr = sorted([(gene,  min(self.t,self._get_node_relevance(gene,self.relevance_function))) for gene in
                                self.G.get_node_list()], key=lambda x: x[1], reverse=True)

            if self.score_function['parameters']['function_ranking_score'] == 'inv_sig':
                page_rank_dictionary = {gene: (self._inverse_sigmoid(rank) if gene not in self.source_node_list else self._inverse_sigmoid(1))  for rank,(gene,score) in enumerate(sorted_pr)}
            else:
                print('No FUNCTION FOUND')
                print(self.score_function['parameters']['function_ranking_score'])
                exit(1)

        total = sum(page_rank_dictionary.values())
        page_rank_dictionary = {k: v / total for k, v in page_rank_dictionary.items()}

        return page_rank_dictionary



    def _get_gene_expression_by_gene_name(self,gene):
        if self.gene_expression_normalization == "z_score":
            return self.gene_expression_manager.find_gene_z_score(gene)

        if "log_z_score_base_" in self.gene_expression_normalization:
            return self._get_log_gene_expression_by_gene_name(gene)

    def _get_log_gene_expression_by_gene_name(self, gene):

        base = float(self.gene_expression_normalization.split("_")[-1])
        gene_expression = self.gene_expression_manager.find_gene_z_score(gene)

        if gene_expression == -1:
            return gene_expression
        else:
            log_gene_expression = []

            for expression in gene_expression:

                log_arg = 1 + abs(expression)
                log_gene_expression.append(math.log(log_arg, base))

            return log_gene_expression


    def _get_gene_expression_function(self,gene_expression):

        if "l_" in self.gene_expression_function:

            p = float(self.gene_expression_function.split("_")[1])
            norm = 0


            for expression in gene_expression:
                norm += abs(expression)**p

            norm = norm**(1/p)

            return norm

        if self.gene_expression_function == "max":
            norm = np.amax(gene_expression)

            if norm < 0:
                norm = 0

            return norm

        if self.gene_expression_function == "absolute_max":

            absolute_vector = np.absolute(gene_expression)
            norm = np.amax(absolute_vector)

            return norm


    def _get_gene_expression_teleporing_probability(self):
        page_rank_dictionary = {}

        for gene in self.G.get_node_list():
            gene_expression = self._get_gene_expression_by_gene_name(gene)

            if gene_expression != -1:
                gene_expression_aggregation_value = self._get_gene_expression_function(gene_expression)

                page_rank_dictionary[gene] = gene_expression_aggregation_value

            else:
                page_rank_dictionary[gene] = 0.0

        total = sum(page_rank_dictionary.values())
        page_rank_dictionary = {k: v / total for k, v in page_rank_dictionary.items()}

        return page_rank_dictionary



    def norm_l1(self,p_t_1,p_t):
        sum = 0
        for gene in p_t_1:
            sum += abs(p_t_1[gene] - p_t[gene])

        return sum

    def _generate_ranked_list(self):
        generate_probabilities = []
        for k,v in self.page_rank_vector.items():
            generate_probabilities.append([k,v])
        sorted_list =  sorted(generate_probabilities, key=lambda x: x[1], reverse=True)
        return sorted_list


    def _compute_next_page_rank(self,p_t):


        residual_j = {}
        for current_gene in p_t.keys():
            sum = 0
            for gene in self.G.get_neighbors(node_name=current_gene):
                sum += p_t[gene] * self.gene_to_gene_weight[gene][current_gene]
            residual_j[current_gene] = sum

        p_t_1 = { gene: (1 - self.restart_probability) * residual_j[gene] + self.restart_probability * self.initial_pr[gene] for gene in p_t.keys()}

        return p_t_1


    def run(self):

        p_v = self._initialize_page_rank()
        self._set_graph_weights()

        if len(p_v) == 0:
            return [],[]
        diff_norm = 1

        while diff_norm > CONV_THRESHOLD:

            p_t_1 = self._compute_next_page_rank(p_v)

            diff_norm = self.norm_l1(p_t_1, p_v)

            p_v = p_t_1

        self.page_rank_vector = p_v

        return self._generate_ranked_list()



