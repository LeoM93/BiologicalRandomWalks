from biological_random_walks.personalization_vector_creation.pv_creation import PersonalizationVectorCreation


class TopologicalPersonalizationVectorCreation(PersonalizationVectorCreation):
	def __init__(self, source, universe, G, secondary_seed_set):
		
		self.universe = universe
		self.G = G

		self.selected_seed_set = self.universe.intersection(source)
		self.source_not_in_G = source.difference(self.selected_seed_set)

		self.secondary_seed_set = secondary_seed_set

	def run(self, ):
		return self._set_up_topological_personalization_vector()


	def _get_radius_2_neighbors(self,node):
		
		neighbors = set(self.G[node])
		neighbors_radius_2_set = set()

		for neighbor in neighbors:
			neighbors_radius_2 =self.G[neighbor]

			for neighbor_radius_2 in neighbors_radius_2:
				neighbors_radius_2_set.add(neighbor_radius_2)

		neighbors_radius_2_set = neighbors_radius_2_set.difference(node)
		neighbors_radius_2_set = neighbors_radius_2_set.difference(neighbors)

		return neighbors_radius_2_set


	def _compute_topological_node_probability(self, node, phi):
		
		node_degree = len(self.G[node])
		neighbors = set(self.G[node])

		neighbors_radius_2 = self._get_radius_2_neighbors(node)
		
		proportion_of_disease_neighbors_radius_1 = len(neighbors.intersection(self.selected_seed_set))/node_degree
		proportion_of_disease_neighbors_radius_2 = len(neighbors_radius_2.intersection(self.selected_seed_set))/len(neighbors_radius_2)

		return proportion_of_disease_neighbors_radius_1 + proportion_of_disease_neighbors_radius_2


	def _set_up_topological_personalization_vector(self,):
		
		assert type(self.secondary_seed_set) == dict, "Secondary seed set is not a dictionaty with Key (string) and value (float)" 
		assert self.G != None and self.secondary_seed_set != None, "Not enough input parameters to compute topological personalization vector" 

		personalization_vector = {}

		for node in self.universe:
			
			if node in self.secondary_seed_set:

				phi = self.secondary_seed_set[node]
				score = self._compute_topological_node_probability(node, phi)
				personalization_vector[node] = score

			else:
				personalization_vector[node] = 0.0

		l_1_personalization_vector = sum(personalization_vector.values())

		personalization_vector = {k:v/l_1_personalization_vector for k,v in personalization_vector.items()}
		return personalization_vector
