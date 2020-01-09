
from multiprocessing import freeze_support
from experiment_pipeline.experiment_initialization import ExperimentInitialization
from experiment_pipeline.network_algorithm_framework import NetworkAlgorithmFramework
from experiment_pipeline.validation_framework import ValidationFramework
from utils.stats import Stats

def run_all():

    run = False
    validate = False
    stats = False


    environment_params = {

        'experiment_id': "brw_journal_experiments_breast_cancer_gene_expression_final_omim_barabasi",
        'disease_dataset_name':"barabasi_omim",


        'disease_ids':[18],
        'ppi_name': "barabasi_ppi",


        'microRNA_db': "mirDB",
        'kegg_db': "kegg_pathways",

        'gene_expression_db': 'TGCA_2014',


        'num_of_fold':100,
        'train_percentuage': 0.70,  # [0,1]

        'ppi_randomization': False,  # True,False
        'term_manager_randomization':None #go, mirna, kegg

    }


    enriched_analysis_params = {
        "db_name":"go",
        "filter_flag":False,
        "min_num_of_genes_for_term_id":0,
        "max_num_of_genes_for_term_id":10000
    }

    experiment_initialization = ExperimentInitialization(environment_params,enriched_analysis_params)
    diseases, ppi_network, db_managers,gene_expression_manager = experiment_initialization.initialize_experiment_variable()


    if run:

        naf = NetworkAlgorithmFramework(environment_params,enriched_analysis_params,diseases,ppi_network,db_managers,gene_expression_manager)
        naf.run_multiprocess()


    if validate:

        metrics_params = {
                'validation_type' : 'train_test_validation',
                'metrics': 'DCG@k', #recall@k, DCG@k
                'disease_module_sizes': [10,20,40,50,60,70,80,90,100,150,200]
            }

        val_f = ValidationFramework(diseases,environment_params,metrics_params)
        val_f.validate_all(name="exploratory_analysis")

        metrics_params = {
            'validation_type': 'train_test_validation',
            'metrics': 'recall@k',  # recall@k, DCG@k
            'disease_module_sizes': [10, 20, 40, 50, 60, 70, 80, 90, 100, 150, 200]
        }

        val_f = ValidationFramework(diseases, environment_params, metrics_params)
        val_f.validate_all(name="exploratory_analysis")

        metrics_params = {
            'validation_type': 'term_manager_validation',
            'term_manager_name': 'go',
            'term_manager': db_managers["go"],
            'filter':'no_filter',
            'metrics':'precision',
            'p_value_threshold': 0.05,
            'disease_module_sizes': [10, 20, 40, 50, 60, 70, 80, 90, 100, 150, 200]

        }
        val_f = ValidationFramework(diseases, environment_params, metrics_params)
        val_f.validate_all(name="exploratory_analysis")

    if stats:
        stats = Stats(diseases,environment_params)

if __name__ == '__main__':
    freeze_support()
    run_all()