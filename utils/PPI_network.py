import networkx as nx

from sklearn.preprocessing import normalize
class PPINetwork():

    def __init__(self,ppi_name,adjacency_matrix,nx_graph = None):
        self.ppi_name = ppi_name
        self.adjacency_matrix = adjacency_matrix


        if nx_graph is not None:
            self.ppi_network = nx_graph
        else:
            self._build_graph_from_adjacency_list()

    def _build_graph_from_adjacency_list(self):
        self.ppi_network = nx.Graph()

        for node, neighbors in self.adjacency_matrix.items():
            for neighbor in neighbors:
                if node != neighbor:
                    self.ppi_network.add_edge(node, neighbor, weight=1.0)



    def get_normalized_adjacency_matrix(self):
        adjacency_matrix_not_normalized = nx.to_numpy_matrix(self.ppi_network)
        return normalize(adjacency_matrix_not_normalized, norm='l1', axis=0)

    def get_adjacency_matrix(self):
        return nx.to_numpy_matrix(self.ppi_network)

    def get_number_of_nodes(self):
        return len(self.ppi_network)

    def get_edge_list(self):
        return self.ppi_network.edges()

    def get_node_index_by_name(self,node_name):
        return self.get_node_list().index(node_name)


    def get_neighbors(self,node_name):
        return self.ppi_network[node_name]


    def get_node_list(self):
        return list(self.ppi_network.nodes())


    #print the node and neighbors of a ppi network
    def print_graph(self):
        for node in self.ppi_network:
            print(str(node) + " -- " + str(self.ppi_network[node]))



