import csv
import numpy as np
import argparse
import math

def __load_df__(file_path, Identifiers):
		
	map__ensembl_id__gene_expression = {}
	csv_reader = csv.reader(open(file_path,"r"),delimiter = "\t")

	for index, row in enumerate(csv_reader):
			
		if index == 0:
			continue

		ensembl_id = row[0].split(".")[0]

		if ensembl_id in map__ensembl_id__gene_expression:
			continue

		if ensembl_id not in Identifiers:
			continue

		v = []

		for i in range(1, len(row)):
				
			score = float(row[i])

			v.append(score)
			
		map__ensembl_id__gene_expression[ensembl_id] = np.array(v)


		
	return map__ensembl_id__gene_expression

def __load_identifier__( file_path):
	csv_reader = csv.reader(open(file_path,"r"),delimiter = "\t")

	Identifiers = set()

	for index, row in enumerate(csv_reader):
			

		Identifiers.add(row[0])

	return Identifiers

def create_de_genes( 
	T_file_path, 
	C_file_path, 
	output_file_path,
	threshold = 2.5,
	
	identifier_file_path = "/Users/leonardomartini/Documents/network_medicine/BRW/data_set/network/HIPPIE_candidate_list.txt",

	):
		
	identifiers = __load_identifier__(identifier_file_path)
		
	tumor_data_frame = __load_df__(T_file_path,identifiers)
	sane_data_frame = __load_df__(C_file_path,identifiers)


	intersection_cols = list(set(tumor_data_frame.keys()).intersection(set(sane_data_frame.keys())))
	table = []
				
	sum_vector = []

	for index,col in enumerate(intersection_cols):
			
		t_v = tumor_data_frame[col]
		s_v = sane_data_frame[col]

		mean = np.mean(s_v)
		std = np.std(s_v)
			
		if std != 0.0:
				
			n_ = np.log(np.absolute((t_v-mean)/std))
			print(n_)
			f = np.where(n_ > threshold, 1, 0)
			sum_ = np.sum(f)
				
			if sum_ > 0:
					
				sum_vector.append(sum_)
				table.append([col, sum_])

	table.sort(key = lambda x: x[1], reverse = True)
	mean = np.mean(sum_vector)

	filtered_table = [record for record in table if record[1] > mean]

	csv_writer = csv.writer(open(output_file_path, "w"),delimiter = "\t")
	csv_writer.writerows(filtered_table)


def __np_pearson_cor__(x, y):
		
	xv = x - x.mean(axis=0)
	yv = y - y.mean(axis=0)
		
	xvss = (xv * xv).sum(axis=0)
	yvss = (yv * yv).sum(axis=0)
		
	result = np.matmul(xv.transpose(), yv) / np.sqrt(np.outer(xvss, yvss))
		
	return np.maximum(np.minimum(result, 1.0), -1.0)[0][0]


def get_top_correlations(
	expression_file_path,
	output_file_path,
	id_file_path ,
	threshold = 0.7):
		
	Identifiers = __load_identifier__(id_file_path)
	df = __load_df__(expression_file_path,Identifiers)
	indeces_df = list(df.keys())
		
	file_name = expression_file_path.split("/")[-1].split("__")[0]

	csv_writer = csv.writer(open(output_file_path,"w"),delimiter = "\t")
	co_expression_network = []

	print("computing pearson's correlation coefficients...")

	for i in range(len(indeces_df)):
		for j in range(i+1, len(indeces_df)):

			v_1 = df[indeces_df[i]]
			v_2 = df[indeces_df[j]]


			pcc_f = __np_pearson_cor__(v_1,v_2)

			if not math.isnan(pcc_f):

				if pcc_f > threshold:
					co_expression_network.append([indeces_df[i],indeces_df[j],abs(pcc_f)])


	csv_writer.writerow(["u","v","score"])
	csv_writer.writerows(co_expression_network)


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	
	parser.add_argument('-T',default = None)
	parser.add_argument('-C',default = None)
	parser.add_argument('-f',default = "../data_set/network/HIPPIE_candidate_list.txt")
	parser.add_argument('-de',default = None)
	parser.add_argument('-co',default = None)

	args = parser.parse_args()

	get_top_correlations(args.T, args.co,args.f)
	create_de_genes(args.T,args.C,args.de,2.5,args.f)



