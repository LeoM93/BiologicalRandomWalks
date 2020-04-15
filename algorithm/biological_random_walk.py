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


        # term manager field

        self.teleporting_p_value = self.algorithm_params["teleporting_parameters"]['p_value']
        self.restart_probability = self.algorithm_params['teleporting_parameters']['restart_probability']
        self.t = self.algorithm_params['teleporting_parameters']['clip_t']
        self.relevance_function = self.algorithm_params['teleporting_parameters']['node_relevance']
        self.score_function = self.algorithm_params['teleporting_parameters']['score_function']
        self.teleporting_p_value_threshold = self.algorithm_params['teleporting_parameters']['p_value_threshold']
        self.fdr_correction = self.algorithm_params['fdr_correction']
        self.teleporting_term_manager = self.algorithm_params["teleporting_parameters"]['term_manager']


        # gene expression field
        self.gene_expression_manager = self.algorithm_params["teleporting_parameters"]["gene_expression"]
        self.gene_expression_normalization = self.algorithm_params["teleporting_parameters"]["gene_expression_normalization"]
        self.gene_expression_function = self.algorithm_params["teleporting_parameters"]["gene_expression_function"]
        self.teleporting_aggregation_function = self.algorithm_params["teleporting_parameters"]["teleporting_aggregation_function"]




        self.walking_p_value = self.algorithm_params["walking_parameters"]['p_value']
        self.walking_term_managet = self.algorithm_params["walking_parameters"]['term_manager']
        self.walking_p_value_threshold = self.algorithm_params['walking_parameters']['p_value_threshold']


        self.gene_expression_biological_walking_aggregation = self.algorithm_params['walking_parameters']['walking_aggregation_function']
        self.gene_expression_pearson_correlation = self.algorithm_params['walking_parameters']['pearson_correlation']
        self.gene_expression_steep = self.algorithm_params['walking_parameters']['pearson_correlation'][0]
        self.gene_expression_translation = self.algorithm_params['walking_parameters']['pearson_correlation'][1]


    def _compute_p_value(self,p_value,p_value_threshold):
        p_value_set = set()
        for term, value in p_value.items():

            if value[0] <= p_value_threshold:
                p_value_set.add(term)

        return p_value_set



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
        return self._get_aggregated_initial_page_rank()


    def _get_aggregated_initial_page_rank(self):

        biological_teleporting = self._biological_teleporting()
        gene_expression_teleporting = self._get_gene_expression_teleporing_probability()

        return self._get_aggregation_function(biological_teleporting,gene_expression_teleporting)



    def _get_aggregation_function(self,pr_1, pr_2):

        aggregated_teleporting = {}

        if self.teleporting_aggregation_function == "sum":

            aggregated_teleporting = {k:  pr_1[k]+pr_2[k] for k in pr_1.keys()}

        elif self.teleporting_aggregation_function == "product":
            aggregated_teleporting = {k:  pr_1[k] * pr_2[k] for k in pr_1.keys()}


        total = sum(aggregated_teleporting.values())
        page_rank_dictionary = {k: v / total for k, v in aggregated_teleporting.items()}

        self.initial_pr = page_rank_dictionary
        return page_rank_dictionary



    def _set_graph_weights(self):

        if self.fdr_correction == 'false':
            self.enriched_walking_set = self._compute_p_value(self.walking_p_value, self.walking_p_value_threshold)

        else:
            self.enriched_walking_set = self._compute_p_value_fdr_correction(self.walking_p_value,
                                                                                 self.walking_p_value_threshold)

        if len(self.enriched_walking_set) == 0:
            self.enriched_walking_set = {}

        self._compute_graph_weight()





    def _compute_edge_relevance(self,intersection,union):


        if self.algorithm_params['walking_parameters']['edge_relevance'] == 0:
            intersection = len(intersection.intersection(self.enriched_walking_set))
            union = 1

        else:
            print()
            print(" Edge Relevance: " + str(self.algorithm_params['walking_parameters']['edge_relevance']) +
                  " is not a admissible option")
            exit(3)

        if union == 0:

            return 0
        else:
            return intersection / union


    def _compute_edge_score(self,edge_relevance_score):
        if self.algorithm_params['walking_parameters']['edge_score'] == 'sum':

            return 1 + edge_relevance_score

        elif self.algorithm_params['walking_parameters']['edge_score'] == 'None':

            return 1

        else:

            print()
            print(" Edge Score: " + self.algorithm_params['walking_parameters']['edge_score'] +
                  " is not a admissible option")
            exit(3)



    def _compute_term_manager_graph_weight(self):
        gene_to_gene_weight = {}

        for gene in self.G.ppi_network:
            gene_to_gene_weight[gene] = {}
            neighbors = self.G.get_neighbors(gene)

            total_weight = 0

            for neighbor in neighbors:

                gene_terms = self.walking_term_managet.find_gene_ontology(gene)
                neighbor_terms = self.walking_term_managet.find_gene_ontology(neighbor)

                if gene_terms == -1 or neighbor_terms == -1:
                    gene_to_gene_weight[gene][neighbor] = 1
                    total_weight += 1

                else:

                    intersection = gene_terms.intersection(neighbor_terms)
                    union = gene_terms.union(neighbor_terms)

                    edge_relevance_score = self._compute_edge_relevance(intersection,union)

                    edge_function_score = self._compute_edge_score(edge_relevance_score)

                    gene_to_gene_weight[gene][neighbor] = edge_function_score
                    total_weight += gene_to_gene_weight[gene][neighbor]


            for neighbor in neighbors:
                gene_to_gene_weight[gene][neighbor] = gene_to_gene_weight[gene][neighbor]/total_weight

        return gene_to_gene_weight


    def _compute_gene_expression_graph_weight(self):

        gene_to_gene_ge_weight = self.gene_expression_manager.load_weighted_gene_expression_adjacency_matrix()

        gene_to_gene_weight = {}
        for gene in self.G.ppi_network:
    
            gene_to_gene_weight[gene] = {}
            neighbors = self.G.get_neighbors(gene)
    
            total_weight = 0.0
    
            for neighbor in neighbors:
    
                if gene in gene_to_gene_ge_weight:
                    if neighbor in gene_to_gene_ge_weight[gene]:

                        gene_to_gene_weight[gene][neighbor] = self._sigmoid(gene_to_gene_ge_weight[gene][neighbor],steep=self.gene_expression_steep,translation=self.gene_expression_translation)
                        total_weight += total_weight


            '''for neighbor in neighbors:
                if gene in gene_to_gene_ge_weight:

                    if neighbor not in gene_to_gene_ge_weight[gene]:
                        gene_to_gene_weight[gene][neighbor] = total_weight / len(neighbors)
                else:
                    gene_to_gene_weight[gene][neighbor] = total_weight / len(neighbors)

            total_weight = sum(gene_to_gene_weight[gene].values())
            '''


            for neighbor in neighbors:
                if neighbor in gene_to_gene_weight[gene]:
                    if total_weight > 0.0:
                        gene_to_gene_weight[gene][neighbor] = gene_to_gene_weight[gene][neighbor] / total_weight
                    else:
                        gene_to_gene_weight[gene][neighbor] = gene_to_gene_weight[gene][neighbor]

        return gene_to_gene_weight


    def _compute_graph_weight(self):

        if self.gene_expression_biological_walking_aggregation is "None":
            self.gene_to_gene_weight = self._compute_term_manager_graph_weight()
        else:
            self.gene_to_gene_weight = self.compute_aggregated_graph_weight()


    def compute_aggregated_graph_weight(self):

        gene_to_gene_term_weight = self._compute_term_manager_graph_weight()
        gene_to_gene_ge_weight = self._compute_gene_expression_graph_weight()

        gene_to_gene_weight = {}

        for gene in self.G.ppi_network:
            gene_to_gene_weight[gene] = {}

            neighbors_ppi = gene_to_gene_term_weight[gene]
            neighbors_correlation_network = {}

            if gene in gene_to_gene_ge_weight:
                neighbors_correlation_network = gene_to_gene_ge_weight[gene]

            total_neighbors = set(neighbors_ppi.keys()).union(set(neighbors_correlation_network))
            total_weight = 0

            for neighbor in total_neighbors:


                if neighbor in neighbors_ppi and neighbor in neighbors_correlation_network:

                    gene_to_gene_weight[gene][neighbor] = self._compute_aggregated_graph_weight( gene_to_gene_term_weight[gene][neighbor],gene_to_gene_ge_weight[gene][neighbor])
                    total_weight += gene_to_gene_weight[gene][neighbor]

                elif neighbor in neighbors_ppi:

                    gene_to_gene_weight[gene][neighbor] = self._compute_aggregated_graph_weight(gene_to_gene_term_weight[gene][neighbor], 0.0)
                    total_weight += gene_to_gene_weight[gene][neighbor]

                else:
                    gene_to_gene_weight[gene][neighbor] = self._compute_aggregated_graph_weight(0.0, gene_to_gene_ge_weight[gene][neighbor])
                    total_weight += gene_to_gene_weight[gene][neighbor]


            for neighbor in total_neighbors:
                if total_weight > 0.0:
                    gene_to_gene_weight[gene][neighbor] = gene_to_gene_weight[gene][neighbor]/total_weight
                else:
                    gene_to_gene_weight[gene][neighbor] = gene_to_gene_weight[gene][neighbor]


        return gene_to_gene_weight


    def _sigmoid(self,x,steep,translation):

        if steep == "_" and translation == "_":
            return x

        elif steep == "_" and translation != "_":
            if x > translation:
                return x
            else:
                return 0

        else:
            return 1 / (1 + math.exp(-steep*(x-translation)))



    def _compute_aggregated_graph_weight(self,biological_weight, gene_expression_weight):

        if self.gene_expression_biological_walking_aggregation == "sum":

            return biological_weight + gene_expression_weight

        elif self.gene_expression_biological_walking_aggregation == "product":

            return  biological_weight * gene_expression_weight




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



        elif relevance_function == 'intersection_normalized_by_degree':

            node_degree = len(self.G.get_neighbors(gene))

            if gene in self.source_node_list:
                return 1/node_degree

            else:
                gene_terms = self.teleporting_term_manager.find_gene_ontology(gene)
                if gene_terms == -1:
                    return 0
                else:
                    return (len(self.enriched_set.intersection(gene_terms))/len(self.enriched_set))/node_degree


        elif relevance_function == 'weighted_intersection':
            gene_terms = self.teleporting_term_manager.find_gene_ontology(gene)

            if gene_terms == -1:
                return 0
            else:
                term_intersection = self.enriched_set.intersection(gene_terms)
                intersection = 0
                for term in term_intersection:

                    genes_involved_in_terms = set(self.teleporting_term_manager.get_gene_symbols_by_ontology_id(term))
                    seed_set = set(self.source_node_list)
                    gene_intersection = genes_involved_in_terms.intersection(seed_set)
                    intersection += len(gene_intersection)

                return intersection

        elif relevance_function == 'fair_intersection':
            gene_terms = self.teleporting_term_manager.find_gene_ontology(gene)

            if gene_terms == -1:
                return 0
            else:
                return len(self.enriched_set.intersection(gene_terms))

        elif relevance_function == 'None':

            if gene in self.source_node_list:
                return 1
            else:
                return 0

        else:
            print("NO NODE RELEVANCE FUNCTION FOUND")
            exit(10)



    def _get_node_score(self,relevance):

        if self.score_function['name'] == 'default':
            return relevance

        elif self.score_function['name'] == 'power':
            return relevance**self.score_function['parameters']['power']

        elif self.score_function['name'] == 'sigmoid':
            return self._sigmoid(relevance,self.score_function['parameters']['steep'],self.score_function['parameters']['translation'])

        elif self.score_function['name'] == 'ranking':
            return relevance


    def _biological_teleporting(self):

        if self.fdr_correction == 'false':
            self.enriched_set = self._compute_p_value(self.teleporting_p_value,self.teleporting_p_value_threshold)

        else:
            self.enriched_set = self._compute_p_value_fdr_correction(self.teleporting_p_value,self.teleporting_p_value_threshold)


        if len(self.enriched_set) == 0:
            return {}


        page_rank_dictionary = {gene:  min(self.t,self._get_node_score(self._get_node_relevance(gene,self.relevance_function))) for gene in self.G.get_node_list() }



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



