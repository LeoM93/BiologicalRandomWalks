import abc

class PersonalizationVectorCreation(metaclass=abc.ABCMeta):
	
	@abc.abstractmethod
	def run(self):
		print("abstract methods")