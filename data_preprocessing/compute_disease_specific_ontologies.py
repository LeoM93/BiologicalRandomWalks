import os
import csv
import argparse

from enrichment_pipeline.enrichment_analysis import EnrichmentAnalysis
from enrichment_pipeline.p_value_correction import compute_p_value_fdr_correction

class DiseaseOntologies():

	def __init__(self, ontology_graph_file_path, disease_seed_file_path, output_file_path ):

		self.output_file_path = output_file_path
		self.disease_seed_file_path = disease_seed_file_path

		self.ontology_graph_file_path = ontology_graph_file_path

		self.p_value_threshold = 1e-5

	def __load_ontology_graph__(self,):

		self.map__db__gene_id__term_ids = {}
		self.map__db__term_id__gene_ids = {}

		with open(self.ontology_graph_file_path,"r") as fp:
			csv_reader = csv.reader(fp,delimiter= "\t")
			
			for index, row in enumerate(csv_reader):
				if index == 0:
					continue

				gene_name = row[0]
				term = row[1]
				db_name = row[2]

				if db_name not in self.map__db__gene_id__term_ids:
					self.map__db__gene_id__term_ids[db_name] = {}

				if db_name not in self.map__db__term_id__gene_ids:
					self.map__db__term_id__gene_ids[db_name] = {}

				if gene_name not in self.map__db__gene_id__term_ids[db_name]:
					self.map__db__gene_id__term_ids[db_name][gene_name] = {term}
				else:
					self.map__db__gene_id__term_ids[db_name][gene_name].add(term)

				if term not in self.map__db__term_id__gene_ids[db_name]:
					self.map__db__term_id__gene_ids[db_name][term] = {gene_name}
				else:
					self.map__db__term_id__gene_ids[db_name][term].add(gene_name)


	def __load_seed__(self,file_path):
		csv_reader = csv.reader(open(file_path,"r"),delimiter = "\t")

		disease_genes = set()

		for row in csv_reader:
			disease_genes.add(row[0])

		return disease_genes


	def run(self,):
		
		
		self.__load_ontology_graph__()

		disease_ontologies = [["Term_ID","DB"]]


		if os.path.exists(self.output_file_path):
			return

		disease_genes = self.__load_seed__(self.disease_seed_file_path)

		for db, map__gene__term_id in self.map__db__gene_id__term_ids.items():

			er = EnrichmentAnalysis(self.map__db__gene_id__term_ids[db],self.map__db__term_id__gene_ids[db],disease_genes)
			terms_p_values = er.get_enirchment_analysis()

			fdr = compute_p_value_fdr_correction(terms_p_values,p_value_threshold = self.p_value_threshold)

			for term in fdr:
				disease_ontologies.append([term, db])


		csv_writer = csv.writer(open(self.output_file_path,"w"),delimiter = "\t")
		csv_writer.writerows(disease_ontologies)
			


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	
	parser.add_argument('-a',default = None)
	parser.add_argument('-s',default = None)
	parser.add_argument('-o',default = None)


	args = parser.parse_args()

d = DiseaseOntologies(
	
	ontology_graph_file_path = args.a, 
	disease_seed_file_path = args.s, 
	output_file_path = args.o

	)

d.run()
