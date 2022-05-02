import json
import unicodecsv as ucsv
import argparse

import csv
import math
import gzip
import shutil

import pandas as pd
import numpy as np

import os


class TCGAAnalyzer():
	
	def __init__(self, sample_sheet_file_path,manifest_file_path, TCGA_directory_path,output_dir_path):
		
		self.sample_sheet_file_path = sample_sheet_file_path
		self.manifest_file_path = manifest_file_path
		self.TCGA_directory_path = TCGA_directory_path
		self.output_dir_path = output_dir_path


	def __load_manifest_files__(self,):
		csv_reader = csv.reader(open(self.manifest_file_path,"r"),delimiter="\t")
		self.manifest_sample_ids = set()
		
		for row in csv_reader:
			
			id_file = row[0]
			self.manifest_sample_ids.add(id_file)

	def __create_mapping_tumor_sane_samples__(self,):
		
		csv_reader = csv.reader(open(self.sample_sheet_file_path,"r"),delimiter="\t")
		self.TCGA_map__project_id___dictionary = {}
		
		for index,row in enumerate(csv_reader):
			if index == 0:
				continue

			file_id = row[0]
			file_name = row[1]
			project_id = row[4]
			case_id = row[5]
			sample_id = row[6]
			sample_type = row[7]


			type_of_tumor = sample_id.split("-")[-1][:2]
			if project_id not in self.TCGA_map__project_id___dictionary:

				self.TCGA_map__project_id___dictionary[project_id] = {
					"T": {},
					"C": {}
				}

			if type_of_tumor == "01":
				self.TCGA_map__project_id___dictionary[project_id]["T"][file_id] = (file_name,file_id)

			if type_of_tumor == "11":
				self.TCGA_map__project_id___dictionary[project_id]["C"][file_id] = (file_name,file_id)

	

	def __get_map__patient__ensembl_id__expression(self, case_control_dict):
		
		map__patient__ensembl_id__expression = {}
		patient_set = set()
		gene_set = set()

		for file_id,record in case_control_dict.items():

			print(file_id)
				
			file_name = record[0]
			case_id = record[1]

			map__patient__ensembl_id__expression[case_id] = {}
			patient_set.add(case_id)

			patient_file_path = self.TCGA_directory_path + file_id + "/" + file_name

			with gzip.open(patient_file_path, 'rb') as f_in:
				csv_reader = ucsv.reader(f_in,delimiter = "\t")
					
				for row in csv_reader:
					
					ensembl_id = row[0]
					RNA_seq = row[1]
					
					gene_set.add(ensembl_id)
					map__patient__ensembl_id__expression[case_id][ensembl_id] = RNA_seq

		return patient_set,gene_set,map__patient__ensembl_id__expression


	def __write_table__(self,patient_set,gene_set,map__patient__ensembl_id__expression,project_id,label):
		M = []

		csv_writer = csv.writer(open(self.output_dir_path + project_id + "__" + label + ".tsv","w"),delimiter = "\t")

		gene_list = list(gene_set)
		gene_list.sort()

		patient_list = list(patient_set)
		patient_list.sort()

		for gene in gene_list:
			record = [gene]

			for patient in patient_list:
				RNA_seq = map__patient__ensembl_id__expression[patient][gene]
				record.append(RNA_seq)

			M.append(record)

		csv_writer.writerow(patient_list)
		csv_writer.writerows(M)

	

	def __remove_null_columns__(self,df):
		dff = pd.DataFrame()

		for cl in df.columns:
			if df[cl].isnull().values.any() == True:
				pass
			else:
				dff[cl] = df[cl]
		return dff
		

	def create_tumor_control_table(self,):

		self.__load_manifest_files__()

		self.__create_mapping_tumor_sane_samples__()

		for project_id, case_control_dict in self.TCGA_map__project_id___dictionary.items():

			print(project_id)

			tumor_patient_set,gene_tumor_set,map__patient_tumor__ensembl_id__expression = self.__get_map__patient__ensembl_id__expression(case_control_dict["T"])
			sane_patient_set,gene_sane_set,map__patient_sane__ensembl_id__expression = self.__get_map__patient__ensembl_id__expression(case_control_dict["C"])

			self.__write_table__(tumor_patient_set,gene_tumor_set,map__patient_tumor__ensembl_id__expression,project_id, "tumor")
			self.__write_table__(sane_patient_set,gene_sane_set,map__patient_sane__ensembl_id__expression, project_id, "control")

	

	

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	
	parser.add_argument('-gdc',default = None)
	parser.add_argument('-m',default = None)
	parser.add_argument('-rna_dir',default = None)

	parser.add_argument('-o',default = None)

	args = parser.parse_args()


	c = TCGAAnalyzer(
		
		sample_sheet_file_path =args.gdc,
		manifest_file_path = args.m,
		TCGA_directory_path = args.rna_dir,
		output_dir_path = args.o
		
		)

	c.create_tumor_control_table()




