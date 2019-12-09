from db_manager.db_manager import DBManager
from utils.preferences import NETWORK_MEDICINE_PATH

class PPIManager():
    def __init__(self,table):
        self.db_table = table
        self.ppi_db_path = NETWORK_MEDICINE_PATH + "db//PPI_Network.db"
        self.ppi_db_manager = DBManager(self.ppi_db_path)




    def print_schema(self):
        self.ppi_db_manager.print_db_schema()


    def select_ppi_network(self):

        adjacency_list = {}
        ppi_query_selection = "SELECT * FROM " + self.db_table
        ppi_adjacency_table = self.ppi_db_manager.comupute_query(ppi_query_selection,record=None)

        for item in ppi_adjacency_table:

            protein_1 = item[0]
            protein_2 = item[1]
            try:
                adjacency_list[protein_1].append(protein_2)

            except KeyError:
                adjacency_list[protein_1] = [protein_2]

        return adjacency_list
