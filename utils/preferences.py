leo = True

if leo:
    NETWORK_MEDICINE_PATH = "/Users/leo/Documents/Phd/data_science/research_field/Network_Medicine/project/netmed/scripts/gene_prioritization/"
else:
    NETWORK_MEDICINE_PATH = "C:/Users/miche/netmed/scripts/gene_prioritization/"



brw = [

{
        "algorithm": "BRW",

        "fdr_correction": [
           1
        ],

        "teleporting_parameters": {
           "name": [
              "biological_teleporting"
           ],

          "gene_expression":[
            "TCGA"
          ],

           "gene_expression_normalization":[
             "log_z_score_base_2"
           ],

           "gene_expression_function":[
             "absolute_max"
           ],
           "teleporting_aggregation_function":[
             "sum"
           ],
           "gene_expression_coefficient":[
             0.5
           ],

           "term_manager_coefficient":[
            0.5
           ],

           "term_manager": [
              "go"
           ],
           "p_value_threshold": [
              0.05
           ],

           "restart_probability": [
              0.9
           ],
           "node_relevance": [
              "weighted_intersection"
           ],
           "clip_t": [
                1
           ],
           "score_function": {
              "name": "default",

              "parameters": {
              }
           }
        },
        "walking_parameters": {
           "name": [
              "biological_walking"
           ],
           "term_manager": [
              "go"
           ],
           "p_value_threshold": [
              0.05
           ],
           "edge_relevance": [
              3
           ],
           "edge_score": [
              "sum"
           ],
           "threshold_pearson_correlation": [
             0.8
           ],

           "walking_aggregation_function":[
             "sum"
           ]
        }
     }
]



rwr = [
    {
        "algorithm": "rwr",
        "restart_probability": 0.75
    }
]

btp =[

    {
      "algorithm": "btp",
      "metrics": "intersection_normalized_by_disease",
      "p_value_threshold": 0.05,
      "term_manager": "go",
      "gene_expression_manager":"TCGA",
      "gene_expression_function":"l_1",
      "gene_expression_normalization":"seed_specific_co_expression_network",
      "gene_expression_aggregation_function":"None",
      "fdr_correction": 1
   }
]

diamond = [ 
    {
        "algorithm": "DIAMOnD",
        "max_num_added_nodes":200

    }
]