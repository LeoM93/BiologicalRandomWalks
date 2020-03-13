
from multiprocessing import freeze_support
from experiment_pipeline.experiment_initialization import ExperimentInitialization
from experiment_pipeline.network_algorithm_framework import NetworkAlgorithmFramework
from experiment_pipeline.validation_framework import ValidationFramework
from utils.stats import Stats

def run_all(experiment_id= None,_run=False, _validate = True, _stats = False):

    run = _run
    validate = _validate
    stats = _stats

    if experiment_id is None:
        experiment_name = "brw_journal_breast_cancer_gene_expression_gwas_lcc_as_seed_breast"
    else:
        experiment_name = experiment_id

    disease_ids = [1]


    environment_params = {

        'experiment_id': experiment_name,
        'disease_dataset_name':"gwas",


        'disease_ids':disease_ids,
        'ppi_name': "barabasi_ppi",


        'microRNA_db': "mirDB",
        'kegg_db': "kegg_pathways",
        'gwas_collection_table':"breast_cancer",
        'gene_expression_db': 'TGCA_2014',


        'num_of_fold':1,
        'train_percentuage': 1,  # [0,1]



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
    diseases, ppi_network, db_managers, gene_expression_manager, gwas_collection_manager = experiment_initialization.initialize_experiment_variable()


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
            'term_manager_name': 'kegg_pathways',
            'term_manager': db_managers["kegg"],
            'filter':'no_filter',
            'metrics':'precision',
            'p_value_threshold': 0.05,
            'term_counter': 1,
            'disease_module_sizes': [10, 20, 40, 50, 60, 70, 80, 90, 100, 150, 200]

        }
        #val_f = ValidationFramework(diseases, environment_params, metrics_params)
        #val_f.validate_all(name="term_counter_1")


    if stats:
        stats = Stats(diseases,environment_params)
        #stats.get_metrics_comparison_over_axis(["BRW", "rwr","DIAMOnD"],["exploratory_analysis_recall@k.csv","exploratory_analysis_DCG@k.csv"])

        stats.get_algorithm_comparison_over_axis(["BRW", "DIAMOnD"],"exploratory_analysis_recall@k.csv",disease_module_size=200)
        stats.get_algorithm_comparison_over_axis(["BRW", "DIAMOnD"],"exploratory_analysis_DCG@k.csv",disease_module_size=200)

if __name__ == '__main__':
    freeze_support()
    run_all(experiment_id=None)