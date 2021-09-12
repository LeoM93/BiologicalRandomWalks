import networkx as nx
from biological_random_walks.matrix_creation.matrix_aggregation import MatrixAggregation

class ConvexCombinationMatrixAggregationCreation(MatrixAggregation):
	
	def __init__(self, PPI_network, CO_expression_network):
		
		self.PPI = PPI_network
		self.CO_expression_network = CO_expression_network

		assert self.PPI != None and self.CO_expression_network != None, "PPI or CO-Expression network are None."

		self.choose_policy_exit = 1



	def choose_policy(self,V,U, chosen_policy = "PPI_network"):
		
		if chosen_policy == "Intersection":
			return V.intersection(U)
		elif chosen_policy == "PPI_network":
			return V
		else:
			print("No corrected choosen policy", chosen_policy)
			exit(self.choose_policy_exit)


	def run(self, chosen_policy, alpha = 0.5):

		PPI_nodes = set(self.PPI.nodes())
		CO_expression_nodes = set(self.CO_expression_network.nodes())


		nodes_from_given_policy = self.choose_policy(PPI_nodes,CO_expression_nodes, chosen_policy = chosen_policy)

		PPI_sub_network = self.PPI.subgraph(nodes_from_given_policy)
		CO_expression_sub_network = self.CO_expression_network.subgraph(nodes_from_given_policy)

		PPI_sub_normalized_sub_network = self._normalize_graph(PPI_sub_network)
		CO_expression_normalized_sub_network = self._normalize_graph(CO_expression_sub_network)

		aggregated_graph = self._aggregate_adjacency_matrix(PPI_sub_normalized_sub_network,CO_expression_normalized_sub_network,nodes_from_given_policy,alpha)

		V = set(aggregated_graph.nodes())
		
		return aggregated_graph, V



	def _aggregate_adjacency_matrix(self,PPI_network,CO_expression_network, V, alpha):
		
		final_graph = nx.DiGraph()

		PPI_edges = PPI_network.edges()
		CO_expression_edges = CO_expression_network.edges()

		for edge in PPI_edges:

			source = edge[0]
			target = edge[1]

			score_ppi_network = PPI_network[source][target]['weight']

			if score_ppi_network > 0.0:
				final_graph.add_edge(source, target, weight = alpha*score_ppi_network)

		for edge in CO_expression_edges:

			source = edge[0]
			target = edge[1]

			score_co_expression_network = CO_expression_network[source][target]['weight']

			if final_graph.has_edge(source,target):
				final_graph[source][target]['weight'] = final_graph[source][target]['weight'] + (1 - alpha)*score_co_expression_network

			else:
				if source in V and target in V:
					
					if score_co_expression_network > 0.0:
						final_graph.add_edge(source, target, weight =  (1 - alpha)*score_co_expression_network)


		return final_graph



	def _print_graph(self,G):
		for edge in G.edges():
			print(edge[0],edge[1],G[edge[0]][edge[1]]['weight'])


	def _normalize_graph(self, G):

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
