import os

class Algorithm(object):
    def __init__(self,source_genes, graph,algorithm_params):
        self.G = graph
        self.source_node_list = source_genes
        self.algorithm_params = algorithm_params

    def run(self):
        raise "abstract method"