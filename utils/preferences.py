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
             "seed_specific_co_expression_network","z_score"
           ],

           "gene_expression_function":[
             "max","l_2"
           ],
           "teleporting_aggregation_function":[
             "union","product"
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
              0.7
           ],
           "node_relevance": [
              "intersection"
           ],
           "clip_t": [
                0.2
           ],
           "score_function": {
              "name": "ranking",

              "parameters": {

                 "function_ranking_score": [
                    "inv_sig"
                 ],
                 "steep": [
                    0.01
                 ],
                 "translation": [
                    100
                 ]
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
              1
           ],
           "edge_score": [
              "sum"
           ],
           "walking_aggregation_function":[
             "sum","None"
           ]
        }
     }
]





rwr = [
    {
        "algorithm": "rwr",
        "restart_probability": 0.25
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
      "gene_expression_aggregation_function":"sum",
      "fdr_correction": 1
   }
]

diamond = [ 
    {
        "algorithm": "DIAMOnD",
        "max_num_added_nodes":200

    }
]