import abc

class MatrixAggregation(metaclass=abc.ABCMeta):
	
	@abc.abstractmethod
	def run(self, chosen_policy, alpha = 0.5):
		print("abstract methods")