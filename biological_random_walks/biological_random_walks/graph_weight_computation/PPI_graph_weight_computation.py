import networkx as nx
class ComputePPIGraphWeight():
	
	def __init__(self, G, disease_ontology = None, map__gene__ontologies = None, c = 1):

		self.constant = c
		self.PPI = G

		self.disease_ontology = disease_ontology
		self.map__gene__ontologies = map__gene__ontologies


	def _get_edge_relevance(self, source, target, database_source, current_disease_ontology):

		if source in self.map__gene__ontologies and target in self.map__gene__ontologies:
			
			if database_source in self.map__gene__ontologies[source] and database_source in self.map__gene__ontologies[target]:

				return len(current_disease_ontology
					.intersection(self.map__gene__ontologies[source][database_source])
					.intersection(self.map__gene__ontologies[target][database_source]))

			else:
				return 0

		else:
			return 0

	def compute_weight_on_graph(self,):

		assert self.disease_ontology != None and self.map__gene__ontologies != None, "Not enough input parameter for computing PPI biological weight"

		weighted_PPI = nx.DiGraph()

		for node in self.PPI:
			neighbors = self.PPI[node]

			for neighbor in neighbors:
				weight = self.constant
				
				for database_source, current_disease_ontology in self.disease_ontology.items():
					weight += self._get_edge_relevance(node, neighbor, database_source, current_disease_ontology)

				
				if weight > 0.0:
					weighted_PPI.add_edge(node, neighbor, weight = weight)

		assert len(weighted_PPI.edges()) == len(self.PPI.edges()) and len(weighted_PPI.nodes()) == len(self.PPI.nodes()), "nodes or edges not overlapping between G and weighted G"

		return weighted_PPI











