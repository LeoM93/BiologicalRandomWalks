import sqlite3

class DBManager:
    def __init__(self,db_path):
        self.db_path = db_path
        self.connection = None
        self.cursor = None

        self._connect_db()

    ######################### Funzioni equivalenti ad ogni singolo DataBase ############################################
    def _connect_db(self):
        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()

    def create_table(self,table_name,field):
        print('''CREATE TABLE ''' + table_name + ''' ''' + field)
        self.cursor.execute('''CREATE TABLE ''' + table_name + ''' ''' + field)
        self.connection.commit()

    def drop_table(self,table_name):
        self.cursor.execute('''DROP TABLE ''' + table_name)
        self.connection.commit()

    def insert_in_table(self,table_name,values):
        self.cursor.execute("INSERT INTO " + table_name + " VALUES " + values)
        self.connection.commit()

    def print_db_schema(self):
        self.cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")

        for id, table in enumerate(self.cursor.fetchall()):
            print("Table N." + str(id) +" "+ table[0])

            self.cursor.execute('select * from ' +table[0])
            descriptions = list(self.cursor.description)
            print()
            for id, description in enumerate(descriptions):

                print("\t" + "Column N." + str(id) + " " + description[0])

            print()

    def comupute_query(self,query,record = None):
        if record is None:
            self.cursor.execute(query)
        else:
            self.cursor.execute(query,record)
        return self.cursor.fetchall()

    def close_connection(self):
        self.connection.close()



    ############################################# Funzioni Diseasome ###################################################
    def select_disease(self,table_name, disease_id):

        tuple = (disease_id,)
        self.cursor.execute("SELECT * FROM " + table_name + " WHERE " + table_name + ".disease_id =?",tuple)
        return list(self.cursor.fetchall()[0])

