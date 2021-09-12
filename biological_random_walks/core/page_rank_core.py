import networkx as nx
CONV_THRESHOLD = 0.000001

class RandomWalkWithRestartCore:

	def __init__(self, 

		personalization_vector,
		G,
		restart_prob = 0.75):

		self.restart_prob = restart_prob
		self.personalization_vector = personalization_vector
		self.G =self.__normalize_graph__(G)


	def __compute_next_page_rank__(self,p_t):
		
		residual_j = {}
		
		for i in p_t.keys():
			sum_ = 0
			
			for j in self.G[i]:
				sum_ += p_t[j] * self.G[j][i]["weight"]
			residual_j[i] = sum_

		p_t_1 = { gene: (1 - self.restart_prob) * residual_j[gene] + self.restart_prob * self.personalization_vector[gene] for gene in p_t.keys()}

		return p_t_1

	def __generate_ranked_list__(self,page_rank_vector):
		
		generate_probabilities = []
		for k,v in page_rank_vector.items():
			generate_probabilities.append([k,v])
		sorted_list =  sorted(generate_probabilities, key=lambda x: x[1], reverse=True)
		return sorted_list

	def __norm_l1__(self,p_t_1,p_t):
		
		sum_ = 0
		for gene in p_t_1:
			sum_ += abs(p_t_1[gene] - p_t[gene])
		return sum_

	def __normalize_graph__(self, G):

		G_normalized = nx.DiGraph()

		for node_1 in G:
			total_weight = 0.0
			
			for node_2 in G[node_1]:
				total_weight += G[node_1][node_2]['weight']

			for node_2 in G[node_1]:
				if total_weight != 0.0:
					G_normalized.add_edge(node_1,node_2, weight = G[node_1][node_2]['weight']/total_weight)
				else:
					G_normalized.add_edge(node_1,node_2, weight = 0.0)

		return G_normalized


	def get_adjacency_matrix(self,):
		adjacency_matrix = []

		for source in self.G:
			for target in self.G[source]:
				adjacency_matrix.append([source, target,self.G[source][target]["weight"]])

		adjacency_matrix.sort(key = lambda x: (x[0],x[1]))
		
		return adjacency_matrix

	def get_p_0(self,):

		p_0 = [[k,v] for k,v in self.personalization_vector.items()]
		p_0.sort(key = lambda x: x[1], reverse = True)
		return p_0


	def run(self):

		p_v = self.personalization_vector

		diff_norm = 1

		while diff_norm > CONV_THRESHOLD:

			p_t_1 = self.__compute_next_page_rank__(p_v)

			diff_norm = self.__norm_l1__(p_t_1, p_v)
			p_v = p_t_1

		return self.__generate_ranked_list__(p_v)