import csv
import argparse


class OntologyGraph():

	def __init__(self, GO_file_path, KEGG_file_path, Reactome_file_path, output_file_path):

		self.GO_file_path = GO_file_path
		self.KEGG_file_path = KEGG_file_path
		self.Reactome_file_path = Reactome_file_path

		self.UniprotKB__Ensembl__mapping_file_path = "../data_set/ontology/uniprot_to_ensembl_mapping.tsv"
		self.KEGG__UniprotKB__mapping_file_path = "../data_set/ontology/KEGG_identifier_to_Uniprot.txt"

		self.output_ontology_network_path = output_file_path
	
	def __load_uniprot_mapping__(self,):
		
		csv_reader = csv.reader(open(self.UniprotKB__Ensembl__mapping_file_path,"r"),delimiter = "\t")
		self.map__uniprotkb__ensembl_id = {}
		
		for index, row in enumerate(csv_reader):
			
			if index == 0:
				continue

			uniprot_id = row[0]
			ensembl_id = row[1]

			if uniprot_id not in self.map__uniprotkb__ensembl_id:
				self.map__uniprotkb__ensembl_id[uniprot_id] = set()


			self.map__uniprotkb__ensembl_id[uniprot_id].add(ensembl_id)


	def __load_KEGG_to_uniprot_mapping__(self,):
		csv_reader = csv.reader(open(self.KEGG__UniprotKB__mapping_file_path,"r"),delimiter = "\t")
		self.map__KEGG__uniprot_id = {}
		
		for index, row in enumerate(csv_reader):
			
			if index == 0:
				continue

			KEGG_id = row[0]
			uniprot_id = row[1]

			if KEGG_id not in self.map__KEGG__uniprot_id:
				self.map__KEGG__uniprot_id[KEGG_id] = set()


			self.map__KEGG__uniprot_id[KEGG_id].add(uniprot_id)




	def __load_go__(self,):
		csv_reader = csv.reader(open(self.GO_file_path,"r"),delimiter = "\t")
		self.map__ensembl_id__gene_ontologies = {}
		
		for row in csv_reader:
			protein_id = row[1]
			ontology_id = row[4]
			validation = row[6]
			bp = row[8]

			if bp == "P" and validation != "IEA":
				
				if protein_id not in  self.map__uniprotkb__ensembl_id:
					continue

				ensembl_ids = self.map__uniprotkb__ensembl_id[protein_id]
				
				for ensembl_id in ensembl_ids:
					
					if ensembl_id not in self.map__ensembl_id__gene_ontologies:
						self.map__ensembl_id__gene_ontologies[ensembl_id] = set()
					
					self.map__ensembl_id__gene_ontologies[ensembl_id].add((ontology_id, "GO"))

				

	def __load_reactome__(self,):
		csv_reader = csv.reader(open(self.Reactome_file_path,"r"),delimiter = "\t")
		
		for row in csv_reader:

			ensembl_id = row[0]
			ontology_id = row[1]
			validation = row[4]
			
			if "R-HSA" in ontology_id and "ENSG" in ensembl_id and validation != "IEA":
				
				if ensembl_id not in self.map__ensembl_id__gene_ontologies:
					self.map__ensembl_id__gene_ontologies[ensembl_id] = set()
					
				self.map__ensembl_id__gene_ontologies[ensembl_id].add((ontology_id, "Reactome"))

			

	
	def __load_KEGG__(self,):
		csv_reader = csv.reader(open(self.KEGG_file_path,"r"),delimiter = "\t")
		
		for row in csv_reader:

			kegg_id = row[1]
			ontology_id = row[0]

			if kegg_id in self.map__KEGG__uniprot_id:
				uniprot_ids = self.map__KEGG__uniprot_id[kegg_id]

				for uniprot_id in uniprot_ids:
					if uniprot_id in self.map__uniprotkb__ensembl_id:
						ensembl_ids = self.map__uniprotkb__ensembl_id[uniprot_id]

						for ensembl_id in ensembl_ids:
							
							if ensembl_id not in self.map__ensembl_id__gene_ontologies:
								self.map__ensembl_id__gene_ontologies[ensembl_id] = set()
					
							self.map__ensembl_id__gene_ontologies[ensembl_id].add((ontology_id, "KEGG"))



	def __save_final_ontology__(self,):
		
		annotation_graph = []
		csv_writer = csv.writer(open(self.output_ontology_network_path, "w"),delimiter = "\t")
		csv_writer.writerow(["gene_id","term_id","DB"])

		for k,v in self.map__ensembl_id__gene_ontologies.items():
			for record in v:
				
				ontology_id, DB = record
				annotation_graph.append([k, ontology_id, DB])

		csv_writer.writerows(annotation_graph)


	def run(self,):
		
		self.__load_uniprot_mapping__()
		self.__load_KEGG_to_uniprot_mapping__()
		
		self.__load_go__()
		self.__load_KEGG__()
		self.__load_reactome__()

		self.__save_final_ontology__()

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	
	parser.add_argument('-go',default = None)
	parser.add_argument('-r',default = None)
	parser.add_argument('-k',default = None)
	parser.add_argument('-o',default = None)

	args = parser.parse_args()	
	
	o = OntologyGraph(

		GO_file_path =args.go, 
		KEGG_file_path=args.k, 
		Reactome_file_path = args.r,
		output_file_path = args.o

		)

	o.run()