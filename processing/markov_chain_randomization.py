import itertools
import random
from utils.PPI_network import PPINetwork
import networkx as nx
class SwitchingAlgorithm():

    def __init__(self,number_of_iterations = 100):

        self.num_iterations = number_of_iterations

    def randomize_by_edge_swaps(self,ppi_network):
        self.PPI = ppi_network
        newgraph = self.PPI.ppi_network.copy()

        edge_list = list(newgraph.edges())
        num_edges = len(edge_list)

        total_iterations = num_edges * self.num_iterations

        for i in range(total_iterations):

            rand_index1 = int(round(random.random() * (num_edges - 1)))
            rand_index2 = int(round(random.random() * (num_edges - 1)))


            original_edge1 = edge_list[rand_index1]
            original_edge2 = edge_list[rand_index2]

            head1, tail1 = original_edge1
            head2, tail2 = original_edge2

            # Flip a coin to see if we should swap head1 and tail1 for
            # the connections
            if random.random() >= 0.5:
                head1, tail1 = tail1, head1

            # The plan now is to pair head1 with tail2, and head2 with
            # tail1
            #
            # To avoid self-loops in the graph, we have to check that,
            # by pairing head1 with tail2 (respectively, head2 with
            # tail1) that head1 and tail2 are not actually the same
            # node. For example, suppose we have the edges (a, b) and
            # (b, c) to swap.
            #
            #   b
            #  / \
            # a   c
            #
            # We would have new edges (a, c) and (b, b) if we didn't do
            # this check.

            if head1 == tail2 or head2 == tail1:
                continue

            # Trying to avoid multiple edges between same pair of nodes;
            # for example, suppose we had the following
            #
            # a   c
            # |*  |           | original edge pair being looked at
            # | * |
            # |  *|           * existing edge, not in the pair
            # b   d
            #
            # Then we might accidentally create yet another (a, d) edge.
            # Note that this also solves the case of the following,
            # missed by the first check, for edges (a, b) and (a, c)
            #
            #   a
            #  / \
            # b   c
            #
            # These edges already exist.

            if newgraph.has_edge(head1, tail2) or newgraph.has_edge(
                    head2, tail1):
                continue

            # Suceeded checks, perform the swap
            original_edge1_data = newgraph[head1][tail1]
            original_edge2_data = newgraph[head2][tail2]

            newgraph.remove_edges_from((original_edge1, original_edge2))

            new_edge1 = (head1, tail2, original_edge1_data)
            new_edge2 = (head2, tail1, original_edge2_data)
            newgraph.add_edges_from((new_edge1, new_edge2))

            # Now update the entries at the indices randomly selected
            edge_list[rand_index1] = (head1, tail2)
            edge_list[rand_index2] = (head2, tail1)

        assert len(newgraph.edges()) == num_edges

        return PPINetwork(ppi_name= self.PPI.ppi_name + "_randomized", adjacency_matrix=None, nx_graph=newgraph)

    def randomize_bipartite_graph(self, db_manager):
        self.db_manager = db_manager
        edges = []

        for gene,annotations in self.db_manager.gene_to_term.items():
            if annotations != -1:
                for annotation in annotations:
                    edges.append((gene,annotation))

        G = nx.DiGraph()
        G.add_edges_from(edges)


        num_edges = len(edges)
        total_iterations = num_edges * self.num_iterations

        for i in range(total_iterations):

            rand_index1 = int(round(random.random() * (num_edges - 1)))
            rand_index2 = int(round(random.random() * (num_edges - 1)))


            original_edge1 = edges[rand_index1]
            original_edge2 = edges[rand_index2]

            head1, tail1 = original_edge1
            head2, tail2 = original_edge2

            if G.has_edge(head1, tail2) or G.has_edge(
                    head2, tail1):
                continue

            original_edge1_data = G[head1][tail1]
            original_edge2_data = G[head2][tail2]

            G.remove_edges_from((original_edge1, original_edge2))

            new_edge1 = (head1, tail2, original_edge1_data)
            new_edge2 = (head2, tail1, original_edge2_data)
            G.add_edges_from((new_edge1, new_edge2))

            # Now update the entries at the indices randomly selected
            edges[rand_index1] = (head1, tail2)
            edges[rand_index2] = (head2, tail1)

        assert len(G.edges()) == num_edges

        self.random_gene_term = {}

        for gene in G:
            if len(G[gene]) != 0:
                for annotation in G[gene]:
                    try:
                        self.random_gene_term[gene].add(annotation)
                    except KeyError:
                        self.random_gene_term[gene] = {annotation}

        return self.random_gene_term