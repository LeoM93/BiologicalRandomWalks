import sys
import numpy as np
from algorithm.algorithm import Algorithm
CONV_THRESHOLD = 0.000001


class RandomWalkWithRestartAlgorithm(Algorithm):
    def __init__(self, source_genes,graph, algorithm_params):

        Algorithm.__init__(self, source_genes,graph,algorithm_params)

        self.normalized_adjacency_matrix = self.G.get_normalized_adjacency_matrix()

        #probability to teleport to one of the source node
        self.restart_probability = self.algorithm_params['restart_probability']



    #the first two methods are related to Kohler et al paper(2008)
    def initialize_page_rank(self):
        #initial page rank vector
        p_v = [0] * self.G.get_number_of_nodes()

        #for each source, the page rank is inizialized to 1 divided the number of sources
        for source_name in self.source_node_list:
            try:
                source_index = self.G.get_node_index_by_name(source_name)
                p_v[source_index] = 1 / float(len(self.source_node_list))

            except ValueError:
                sys.exit("Source node is not in the graph".format(source_name,self.source_node_list))

        return np.array(p_v)

    #compute page rank vector at time t + 1
    def compute_next_page_rank(self,p_t,p_v):
        epsilon = np.squeeze(np.asarray(np.dot(self.normalized_adjacency_matrix, p_t)))

        no_restart = epsilon * (1 - self.restart_probability)
        restart = p_v * self.restart_probability
        return np.add(no_restart, restart)


    def _generate_ranked_list(self):
        generate_probabilities = zip(self.G.get_node_list(), self.page_rank_vector.tolist())
        sorted_list =  sorted(generate_probabilities, key=lambda x: x[1], reverse=True)
        return sorted_list


    def run(self):

        p_v = self.initialize_page_rank()
        p_t = np.copy(p_v)

        diff_norm = 1
        while (diff_norm > CONV_THRESHOLD):
            # first, calculate p^(t + 1) from p^(t)
            p_t_1 = self.compute_next_page_rank(p_t, p_v)

            # calculate L1 norm of difference between p^(t + 1) and p^(t),
            # for checking the convergence condition
            diff_norm = np.linalg.norm(np.subtract(p_t_1, p_t), 1)

            # then, set p^(t) = p^(t + 1), and loop again if necessary
            # no deep copy necessary here, we're just renaming p
            p_t = p_t_1

        self.page_rank_vector = p_t
        return self._generate_ranked_list()
