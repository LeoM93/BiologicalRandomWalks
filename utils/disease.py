class Disease():

    def __init__(self,disease_id,name,genes,chromosomes,disease_class = None):
        self.disease_id = disease_id
        self.disease_name =name
        self.disease_genes =genes
        self.disease_chromosomes = chromosomes
        self.disease_class = disease_class

    def __str__(self):
        name = self.disease_name


        if type(self.disease_genes[0]) is int:

            disease_genes = ' '.join(str(item) for item in self.disease_genes)
        else:
            disease_genes = ' '.join(self.disease_genes)

        return str(self.disease_id) + "\t" + name + "\t" + disease_genes


