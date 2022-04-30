from biological_random_walks.personalization_vector_creation.pv_creation import PersonalizationVectorCreation


class BiologicalPersonalizationVectorCreation(PersonalizationVectorCreation):
	
	def __init__(self, source, universe, disease_ontology = None, map__gene_name__ontologies = None, universe_ontologies = None):
		
		self.universe = universe

		self.selected_seed_set = self.universe.intersection(source)
		self.source_not_in_G = source.difference(self.selected_seed_set)


		self.disease_ontology = disease_ontology
		self.map__gene_name__ontologies = map__gene_name__ontologies

	
	def run(self, ):
		return self._set_up_biological_personalization_vector()


	def _set_up_biological_personalization_vector(self,discriminant = True):
		
		assert self.map__gene_name__ontologies != None and self.disease_ontology != None, "Not enough input parameters for biological teleporting probability"
		
		assert len(self.selected_seed_set) != 0, "No source gene in PPI network"
		personalization_vector = {}
		for node in self.universe:

			if node in self.selected_seed_set and discriminant:
				personalization_vector[node] = len(self.disease_ontology)
				continue

			if node in self.map__gene_name__ontologies:
				node_relevance = 0.0
				node_ontologies = self.map__gene_name__ontologies[node]
					
				#node relevance contribution over a set of different biological information
				for k,v in self.disease_ontology.items():
					#node relevance
					if k in node_ontologies:
						node_relevance += len(v.intersection(node_ontologies[k]))/len(v)

				personalization_vector[node] = node_relevance
			else:
				personalization_vector[node] = 0.0
	
		l_1_norm = sum(personalization_vector.values())
		assert l_1_norm > 0.0, "personalization vector is the null vector"

		normalized_personalization_vector = {k: v / l_1_norm for k, v in personalization_vector.items()}

		return normalized_personalization_vector
