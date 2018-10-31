
Introduction
============

A key factor in molecular dynamics simulations is the consistency and reliability of the molecular parameters, especially the partial atomic charges. charge_assign efficiently, reliably and automatically assigns partial atomic charges to atoms based on known distributions.

Basic Usage
-----------

charge_assign requires a repository of molecules with previously computed partial charges. The repository used by us (the authors) was built from the `ATB database <https://atb.uq.edu.au/>`_ . Please contact us if you are interested in using this repository.

If you want to use your own repository, it can be esily built from a directory of .lgf files (See :doc:`this </conversion>` on how to convert to .lgf). A word of caution however, this operation might consume a lot of RAM::

	repo = Repository.create_from('/path/to/dir')
	repo.write('repo.zip')

Once you have obtained a repository, the basic usage of charge_assign is simple (See :doc:`this </conversion>` on how to convert a molecule from a variety of different formats)::
	
	repo = Repository.read('repo.zip')
	charger = CDPCharger(repository=repo)
	charger.charge(molecule)
	for u, data in mol.nodes(data=True):
		print(u, data['partial_charge_redist'])

charge_assign offers several different possible charger methods, see the :mod:`charger <charge.chargers>` module.


Batch Computations
------------------
	
If you are using multiple chargers, you might want to consider using one :class:`Nauty <charge.nauty.Nauty>` instance to reduce the number of processes::

	nauty = Nauty()
	repo = Repository.read('repo.zip', nauty=nauty)
	charger1 = CDPCharger(repository=repo, nauty=nauty)
	charger2 = ILPCharger(repository=repo, nauty=nauty)
	charger3 = MeanCharger(repository=repo, nauty=nauty)
	

Processing several molecules in parallel can be done using a :class:`MultiProcessor <charge.multiprocessor.MultiProcessor>`::

	class Worker:
		def __init__(self, repo_location: str):
			repo = Repository.read(repo_location)
			self.__charger = CDPCharger(repository=repo)

		def process(self, molecule: nx.Graph) -> nx.Graph:
			charger.charge(molecule)
			return molecule

	with MultiProcessor(Worker, 'repo.zip') as mp:
		for c in mp.processed(molecules, 'processing molecules'):
			pass # do something with the charged molecules
