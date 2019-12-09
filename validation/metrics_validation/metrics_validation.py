from validation.metrics_validation.metrics.metrics import DCG,recall_at_k

class MetricsValidation():
    def __init__(self,ranked_list,test_set,disease_module_size):
        self.ranked_list = ranked_list
        self.test_set = test_set
        self.disease_module_size = disease_module_size



    def switch_metrics(self,metrics,ranked_list):


        if metrics == "DCG@k":
            return DCG(ranked_list, self.test_set, self.disease_module_size)


        if metrics == "recall@k":
            return recall_at_k(ranked_list, self.test_set,self.disease_module_size)


    def compute_metrics(self,metrics):

        return self.switch_metrics(metrics,self.ranked_list)



