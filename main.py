from biological_random_walks.BiologicalRandomWalks import BiologicalRandomWalks
import os
import argparse

if __name__ == '__main__':

	code = 1
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-s',default = None)
	parser.add_argument('-de',default = None)


	parser.add_argument('-p',default = None)	
	parser.add_argument('-c',default = None)

	parser.add_argument('-do',default = None)

	parser.add_argument('-a',default = None)

	parser.add_argument('-o',default = None)
	parser.add_argument('-r',default = 0.75)


	args = parser.parse_args()
	personalization_vector_creation_policies = []

	if args.s is None:
		
		print("No Seed dected")
		exit(code)

	if args.p is None and args.c is None:
		print("No Networks dected")
		exit(code)

	if args.p is not None:
		ppi_file_path = args.p
	else:
		ppi_file_path = None
	
	if args.c is not None:
		co_expression_file_path = args.c
	else:
		co_expression_file_path = None

	if args.p is not None and args.c is not None:
		matrix_aggregation_policy = "convex_combination"
	
	elif args.p is not None:
		matrix_aggregation_policy = "only_ppi_network"
	
	else:
		matrix_aggregation_policy = "only_co_expression_network"


	if args.s is not None:
		seed_file_path = args.s

	if args.de is not None:
		secondary_seed_file_path = args.de
		personalization_vector_creation_policies.append("topological")
	else:
		secondary_seed_file_path = None


	if args.do is not None and args.a is not None:
		personalization_vector_creation_policies.append("biological")
		network_weight_flag = True
		
		ontologies_path = args.a
		disease_ontology_path = args.do

	else:

		network_weight_flag = False
		ontologies_path = None
		disease_ontology_path = None
	
	if args.o is not None:
		output_file_path = args.o
	else:
		output_file_path = None

	if len(personalization_vector_creation_policies) == 0:
		personalization_vector_creation_policies.append("default")

	brw = BiologicalRandomWalks(
			
		seed_file_path = seed_file_path,
		secondary_seed_file_path = secondary_seed_file_path,

		ppi_file_path = ppi_file_path,
		co_expression_file_path = co_expression_file_path,

		matrix_aggregation_policy = matrix_aggregation_policy,


		disease_ontology_file_path= disease_ontology_path,
		map__gene__ontologies_file_path = ontologies_path,

		personalization_vector_creation_policies = personalization_vector_creation_policies,
		personalization_vector_aggregation_policy = "Sum",

		restart_prob = args.r,

		network_weight_flag = network_weight_flag,

		output_file_path = output_file_path
	)
