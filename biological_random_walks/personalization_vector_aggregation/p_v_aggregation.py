class PersonalizationVectorAggregation():
	
	def __init__(self, personalization_vectors,universe):

		map__index__p_vs = {}

		self.universe = universe
		self.map__index__p_vs = {i:p_v for i,p_v in enumerate(personalization_vectors)}
		


	def run(self,chosen_policy = "Sum"):
		
		aggregated_personalization_vector = {}
	
		for node in self.universe:
			aggregated_personalization_vector[node] = 0.0

			for i, p_v in self.map__index__p_vs.items():
				
				if chosen_policy == "Sum":
					aggregated_personalization_vector[node] += p_v[node]
				elif chosen_policy == "Product":
					aggregated_personalization_vector[node] *= p_v[node]

		l_1 = sum(aggregated_personalization_vector.values())

		personalization_vectors = {k: v/l_1 for k,v in aggregated_personalization_vector.items()}

		return personalization_vectors










