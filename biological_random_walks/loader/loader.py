import csv
import networkx as nx

class Loader():
	def __init__(self, 
		
		ppi_file_path = None, 
		co_expression_file_path = None, 
		
		seed_file_path = None,
		secondary_seed_file_path = None,
		
		disease_ontology_file_path = None,
		map_gene_ontologies_file_path = None,


		):
		if ppi_file_path != None:
			self.ppi_file_path = ppi_file_path
		else:
			self.ppi_file_path = None

		if co_expression_file_path != None:
			self.co_expression_file_path = co_expression_file_path
		else:
			self.co_expression_file_path = None
		
		if seed_file_path != None:
			self.seed_file_path = seed_file_path
		else:
			self.seed_file_path = None

		if secondary_seed_file_path != None:
			self.secondary_seed_file_path = secondary_seed_file_path
		else:
			self.secondary_seed_file_path = None

		
		if disease_ontology_file_path != None:
			self.disease_ontology_file_path = disease_ontology_file_path
		else:
			self.disease_ontology_file_path = None

		if map_gene_ontologies_file_path != None:
			self.map_gene_ontologies_file_path = map_gene_ontologies_file_path
		else:
			self.map_gene_ontologies_file_path = None



	def run(self,):
		
		assert self.ppi_file_path != None or self.co_expression_file_path != None, "No network as input of Random Walks"
		assert self.seed_file_path != None, "No Seed as input of Random Walks"

		if self.ppi_file_path != None:
			PPI = self.load_graph(self.ppi_file_path)
		else:
			PPI = None

		if self.co_expression_file_path:
			CO_expression = self.load_graph(self.co_expression_file_path)
		else:
			CO_expression = None
		
		if self.seed_file_path != None:
			seed_set = self.load_seed_set(self.seed_file_path)

		if self.map_gene_ontologies_file_path != None:
			map__gene__ontologies = self.load_map__gene__ontologies()
		else:
			map__gene__ontologies = None

		if self.disease_ontology_file_path != None:
			disease_ontology = self.load_disease_ontology()
		else:
			disease_ontology = None

		if self.secondary_seed_file_path != None:
			secondary_seed_set = self.load_seed_set(self.secondary_seed_file_path)
		else:
			secondary_seed_set = None


		return PPI, CO_expression, seed_set, secondary_seed_set, map__gene__ontologies, disease_ontology


	def load_map__gene__ontologies(self):
		
		map_gene_ontologies = {}
		
		with open(self.map_gene_ontologies_file_path,"r") as fp:
			csv_reader = csv.reader(fp, delimiter = "\t")

			for index,row in enumerate(csv_reader):
				if index == 0:
					continue
				
				gene_name = row[0]
				term_id = row[1]
				db = row[2]

				if gene_name in map_gene_ontologies:
					
					if db in map_gene_ontologies[gene_name]:
						map_gene_ontologies[gene_name][db].add(term_id)
					else:
						map_gene_ontologies[gene_name][db] = {term_id}
				
				else:
					map_gene_ontologies[gene_name] = {}
					map_gene_ontologies[gene_name][db] = {term_id}

		
		return map_gene_ontologies

	def load_disease_ontology(self,):
		disease_ontology = {}

		with open(self.disease_ontology_file_path, 'r') as fp:
			csv_reader = csv.reader(fp,delimiter = "\t")

			for index, row in enumerate(csv_reader):

				if index == 0:
					continue

				ontology = row[0]
				db = row[1]

				if db in disease_ontology:
					disease_ontology[db].add(ontology)
				else:
					disease_ontology[db] = {ontology}

		return disease_ontology


	def load_graph(self,file_path, has_header = False, absolute_policy = True):
		G = nx.DiGraph()

		with open(file_path, 'r') as fp:
			csv_reader = csv.reader(fp, delimiter = "\t")

			for index, row in enumerate(csv_reader):

				if has_header:
					if index == 0:
						continue

				if len(row) == 3:

					node_1 = row[0]
					node_2 = row[1]

					try:
						score = float(row[2])
					except ValueError:
						score = 0.0

					if absolute_policy:
						score = abs(score)
					else:
						score = score

					if not G.has_edge(node_1,node_2) and not G.has_edge(node_2,node_1):
						G.add_edge(node_1,node_2,weight = score)
						G.add_edge(node_2,node_1,weight = score)

				elif len(row) == 2:
					
					node_1 = row[0]
					node_2 = row[1]

					if not G.has_edge(node_1,node_2) and not G.has_edge(node_2,node_1):
						
						G.add_edge(node_1,node_2,weight = 1.0)
						G.add_edge(node_2,node_1,weight = 1.0)
		return G


	def load_seed_set(self,file_path):
		
		seed_set = set()
		seed_dict = {}

		column_size = 0

		with open(file_path, 'r') as fp:
			csv_reader = csv.reader(fp,delimiter = "\t")

			for index,row in enumerate(csv_reader):

				if index == 0:
					column_size = len(row)

				if column_size == 1:

					seed = row[0]
					seed_set.add(seed)

				elif column_size == 2:
					seed = row[0]
					score = float(row[1])

					seed_dict[seed] = score


		if column_size == 1:
			return seed_set

		if column_size == 2:
			return seed_dict

