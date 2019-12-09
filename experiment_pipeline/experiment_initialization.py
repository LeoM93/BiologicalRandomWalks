from utils.preferences import NETWORK_MEDICINE_PATH,rwr,brw,btp,diamond
from utils.PPI_network import PPINetwork

from db_manager.ppi_db_manager.ppi_db_manager import PPIManager
from db_manager.disease_db_manager.disease_db_manager import DiseaseManager
from db_manager.gene_expression_db_manager.GEManager import GEManager
from db_manager.term_db_manager.gene_ontology_manager import GeneOntologyManager
from db_manager.term_db_manager.mirna_manager import MirnaManager
from db_manager.term_db_manager.kegg_manager import KeggManager
from db_manager.term_db_manager.random_term_manager import RandomTermManager


from enrichment_analysis.enrichment_analysis import EnrichmentAnalysis

from processing.preprocessing import check_input_disease,split_train_test_set
from processing.markov_chain_randomization import SwitchingAlgorithm
from file_manager.file_writer import write_data_on_disk,WriteEnrichmentAnalysis,write_json_on_disk
from file_manager.file_loader import load_train_test_file
import os


class ExperimentInitialization():

    def __init__(self,environment_params,enriched_analysis_params):
        self.ppi_name = environment_params["ppi_name"]

        self.ppi_randomization = environment_params["ppi_randomization"]

        # disease ids vector
        self.disease_ids = environment_params["disease_ids"]

        # num of trial
        self.num_of_trial = environment_params["num_of_fold"]

        # experiment id
        self.experiment_id = environment_params["experiment_id"]

        # omim name
        self.omim_name = environment_params["disease_dataset_name"]

        #microRNA database to use
        self.microRNA_db = environment_params['microRNA_db']

        #kegg pathways database to use
        self.kegg_db = environment_params['kegg_db']

        # gene expression
        self.gene_expression_db = environment_params["gene_expression_db"]


        self.term_manager_randomization = environment_params["term_manager_randomization"]
        # flag for filtering annotations
        self.filter_annotation = enriched_analysis_params["filter_flag"]

        self.min_number_of_gene_for_term_id = enriched_analysis_params["min_num_of_genes_for_term_id"]

        self.max_number_of_gene_for_term_id = enriched_analysis_params["max_num_of_genes_for_term_id"]

        # experiment path
        self.experiment_directory_path = NETWORK_MEDICINE_PATH + "exps/" + self.experiment_id

        # train percentuage
        self.train_percentuage = environment_params["train_percentuage"]


    def _compute_enrichment_analysis_path(self,disease_id_str,db_manager):
        enrichment_analysis_path = self.ppi_directory_path + "/" + disease_id_str + "/enrichment_analysis/" + db_manager.name + "/"

        if self.filter_annotation:
            enrichment_analysis_final_path = enrichment_analysis_path + "filter/" + str(
                self.min_number_of_gene_for_term_id) + "_" + str(self.max_number_of_gene_for_term_id) + "/"

        else:
            enrichment_analysis_final_path = enrichment_analysis_path + "no_filter/"


        if not os.path.exists(enrichment_analysis_path):
            os.makedirs(enrichment_analysis_final_path)

        return enrichment_analysis_final_path



    def initiazlize_db_manager(self,node_list):

        gene_ontology_db_manager = self._initialize_gene_ontology_db(node_list)
        mirna_db_manager = self._initialize_microRNA_db(node_list)
        kegg_manager = self._initialize_kegg_db(node_list)

        gene_ontology_db_manager.close_connection()
        mirna_db_manager.close_connection()
        kegg_manager.close_connection()

        if self.term_manager_randomization != None:
            randomize_db = SwitchingAlgorithm()
            random_switching_phase = None
            random_term_manager = None

            if self.term_manager_randomization == "go":
                random_term_manager = RandomTermManager("go_random")
                random_switching_phase = randomize_db.randomize_bipartite_graph(gene_ontology_db_manager)

            if self.term_manager_randomization == "mirna":
                random_term_manager = RandomTermManager("mirna_random")
                random_switching_phase = randomize_db.randomize_bipartite_graph(mirna_db_manager)

            if self.term_manager_randomization == "kegg":
                random_term_manager = RandomTermManager("kegg_random")
                random_switching_phase = randomize_db.randomize_bipartite_graph(kegg_manager)

            random_term_manager.init_db(node_list,random_switching_phase)

        else:
            random_term_manager = None

        db_managers = {
            "go": gene_ontology_db_manager,
            "mirna": mirna_db_manager,
            "kegg":kegg_manager,
            "random":random_term_manager
        }

        return db_managers

    def _initialize_gene_ontology_db(self,node_list):
        gene_ontology_manager = GeneOntologyManager(name="go")
        gene_ontology_manager.init_db(node_list)
        gene_ontology_manager.print_stats(node_list)

        return gene_ontology_manager



    def _initialize_microRNA_db(self,node_list):
        mirna_manager = MirnaManager(table = self.microRNA_db,name = self.microRNA_db)
        mirna_manager.init_db(node_list)

        return mirna_manager


    def _initialize_kegg_db(self,node_list):
        kegg_manager = KeggManager(table=self.kegg_db, name= self.kegg_db)
        kegg_manager.init_db(node_list)

        return kegg_manager

    def _initialize_gene_expression(self,ppi_network_name, _list):

        ge = GEManager(table=self.gene_expression_db,ppi_name = ppi_network_name)
        ge.init_db(_list)
        ge.close_connection()

        return ge

    def _initialize_diseases(self,node_list):
        disease_manager = DiseaseManager(self.omim_name)


        mapped_diseases = []

        if len(self.disease_ids) == 0:
            diseases = disease_manager.select_all_disease()

            for disease in diseases:
                mapped_disease = check_input_disease(disease, node_list)
                mapped_diseases.append(mapped_disease)

            return mapped_diseases


        else:

            for disease_id in self.disease_ids:
                disease = disease_manager.select_disease(disease_id)

                mapped_disease = check_input_disease(disease,node_list)

                mapped_diseases.append(mapped_disease)

            return mapped_diseases



    def _initialize_ppi_network(self):
        ppi_network_manager = PPIManager(self.ppi_name)
        adjacency_list = ppi_network_manager.select_ppi_network()

        ppi_network = PPINetwork(ppi_name=self.ppi_name, adjacency_matrix=adjacency_list)


        if self.ppi_randomization:
            randomization = SwitchingAlgorithm()
            randomize_ppi_network = randomization.randomize_by_edge_swaps(ppi_network)

            return randomize_ppi_network
        else:
            return ppi_network




    def _initialize_disease_experiment_tree(self,diseases):


        for disease in diseases:

            disease_id = disease.disease_id

            disease_id_str = str(disease_id)

            train_path = self.ppi_directory_path + "/" + disease_id_str +"/train"

            if not os.path.exists(train_path):
                os.makedirs(train_path)

            test_path = self.ppi_directory_path + "/" + disease_id_str + "/test"

            if not os.path.exists(test_path):
                os.makedirs(test_path)

            validation_output_path = self.ppi_directory_path + "/" + disease_id_str +"/validation_outputs"

            if not os.path.exists(validation_output_path):
                os.makedirs(validation_output_path)


            # output directory of each algorithms
            output_algorithm_path = self.ppi_directory_path + "/" + disease_id_str + "/algo_outs"

            if not os.path.exists(output_algorithm_path):
                os.makedirs(output_algorithm_path)


    def _initialize_experiment_tree(self,diseases):
        if not os.path.exists(self.experiment_directory_path):
            os.makedirs(self.experiment_directory_path)


        if self.ppi_randomization:
            self.ppi_directory_path = self.experiment_directory_path + "/" + self.ppi_name + "_randomized"

        else:
            self.ppi_directory_path = self.experiment_directory_path + "/" + self.ppi_name

        if not os.path.exists(self.ppi_directory_path):
            os.makedirs(self.ppi_directory_path)

        self.default_experiment_variables_path = self.experiment_directory_path + "/experiment_variables/"
        if not os.path.exists(self.default_experiment_variables_path):
            os.makedirs(self.default_experiment_variables_path)


        self._initialize_disease_experiment_tree(diseases)



    def _initialize_train_test_set(self,mapped_diseases,ppi_network, db_manager):
        for disease in mapped_diseases:

            disease_id_str = str(disease.disease_id)

            # train and test path for given disease id
            train_path = self.ppi_directory_path + "/" + disease_id_str + "/train"
            test_path = self.ppi_directory_path + "/" + disease_id_str + "/test"

            # enrichment analysis path for given disease id
            gene_ontology_enrichment_analysis_path = self._compute_enrichment_analysis_path(disease_id_str,db_manager["go"])
            mirna_enrichment_analysis_path = self._compute_enrichment_analysis_path(disease_id_str,db_manager["mirna"])

            if db_manager["random"] != None:
                random_enrichment_analysis_path = self._compute_enrichment_analysis_path(disease_id_str,db_manager["random"])


            # if there are no element in train directory and test directory create a number fo train test file
            # depending on the number of trial
            starting_point = 0

            if len(os.listdir(train_path)) != len(os.listdir(test_path)) or len(os.listdir(train_path)) == 0:
                starting_point = 0

            elif len(os.listdir(train_path)) < self.num_of_trial:
                starting_point  = int(max(os.listdir(train_path)))+1
                print("Starting from trial n: ", starting_point)

            else:
                print("No FILE WILL BE CREATED BECAUSE N. OF FOLDS LOWER/GREATER FROM N. OF TRAIN")
                continue


            for i in range(starting_point, self.num_of_trial):

                # split train and test
                train_set, test_set = split_train_test_set(disease.disease_genes,self.train_percentuage)
                train_file_path = train_path + "/" + str(i)
                test_file_path = test_path + "/" + str(i)

                # write on data train and test
                train_list = [[x] for x in train_set]
                write_data_on_disk(train_file_path, train_list, delimiter=",", write_mode="wb")
                test_list = [[x] for x in test_set]
                write_data_on_disk(test_file_path, test_list, delimiter=",", write_mode="wb")

                # compute enrichment analysis
                seed_set = load_train_test_file(train_file_path)

                self._write_gene_ontology_p_value(disease,gene_ontology_enrichment_analysis_path,db_manager["go"],ppi_network,seed_set,str(i))
                self._write_gene_ontology_p_value(disease,mirna_enrichment_analysis_path,db_manager["mirna"],ppi_network,seed_set,str(i))

                if db_manager["random"] != None:
                    self._write_gene_ontology_p_value(disease,random_enrichment_analysis_path,db_manager["random"],ppi_network,seed_set,str(i))



    def _initialize_algorithm_variables(self):

        rwr_file_path = self.default_experiment_variables_path + "rwr.json"

        if not os.path.exists(rwr_file_path):
            write_json_on_disk(rwr_file_path,rwr )


        diamond_file_path = self.default_experiment_variables_path + "diamond.json"

        if not os.path.exists(diamond_file_path):
            write_json_on_disk(diamond_file_path, diamond)

        btp_file_path = self.default_experiment_variables_path + "btp.json"

        if not os.path.exists(btp_file_path):
            write_json_on_disk(btp_file_path, btp)

        brw_file_path = self.default_experiment_variables_path + "brw.json"

        if not os.path.exists(brw_file_path):
            write_json_on_disk(brw_file_path, brw)



    def initialize_experiment_variable(self):


        # create the ppi network
        ppi_network = self._initialize_ppi_network()

        ppi_node_list = ppi_network.get_node_list()
        ppi_edge_list = ppi_network.get_edge_list()
        ppi_name = ppi_network.ppi_name

        # retrieve the disease taken into account in the current experiment
        mapped_diseases = self._initialize_diseases(ppi_node_list)

        # gene ontology db manager
        db_manager =  self.initiazlize_db_manager(ppi_node_list)

        # create the sub directory of the experiment
        self._initialize_experiment_tree(mapped_diseases)

        # create train test splits and enriched analysis
        self._initialize_train_test_set(mapped_diseases,ppi_network,db_manager)

        # initialize z_score
        gene_expression_manager = self._initialize_gene_expression(ppi_name,ppi_edge_list)

        self._initialize_algorithm_variables()

        return mapped_diseases,ppi_network,db_manager,gene_expression_manager






    def _write_gene_ontology_p_value(self,disease,enrichment_file_path,db_manager, ppi_network, seed_set, train_trial):

        # enrichment analysis
        write_enrichment_analysis = WriteEnrichmentAnalysis(enrichment_file_path)
        # check if there is the file containing the enrichment analysis given a particular disease
        if write_enrichment_analysis.check_enrichment_analysis(train_trial=train_trial) is False:


            
            if self.filter_annotation:

                enrichment_analysis = EnrichmentAnalysis(db_manager=db_manager, universe=ppi_network,
                                                         disease=disease,
                                                         min_num_of_genes_for_term_id=self.min_number_of_gene_for_term_id,
                                                         max_number_of_gene_for_term_id=self.max_number_of_gene_for_term_id)

            else:

                enrichment_analysis = EnrichmentAnalysis(db_manager=db_manager, universe=ppi_network,
                                                         disease=disease, min_num_of_genes_for_term_id=0,
                                                         max_number_of_gene_for_term_id=100000)

            

            p_value_by_term_id = enrichment_analysis.get_enirchment_analysis(seed_set)
            write_enrichment_analysis.write_enrichment_analysis(p_value_by_term_id, train_trial)











