from biological_random_walks.personalization_vector_creation.pv_creation import PersonalizationVectorCreation

class DefaultPersonalizationVectorCreation(PersonalizationVectorCreation):

	
	def __init__(self, source, universe):
		
		self.universe = universe
		self.selected_seed_set = self.universe.intersection(source)
		self.source_not_in_G = source.difference(self.selected_seed_set)

	
	def run(self, ):
		return self._set_up_default_personalization_vactor()

	
	def _set_up_default_personalization_vactor(self,):
		
		personalization_vector = {}
		size = len(self.selected_seed_set)

		# seed set different from empty set
		assert size != 0, ",".join(list(self.source_not_in_G)) + " are not in G"
		
		for node in self.universe:
			
			if node in self.selected_seed_set:
				personalization_vector[node] = 1.0/size
			else:
				personalization_vector[node] = 0.0

		return personalization_vector