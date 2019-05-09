
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

To set a default repository location, set the REPO_LOCATION environment variable or add it to charge/settings.py.

Chargers
--------

The charger classes are the main "workhorses" in charge_assign. For each atom in the query molecule, the charger hashes it's neighborhood graph, queries the repository for this hash and finally processes and assigns the returned partial charge values to the atoms. charge_assign offers several different possible charger classes, see the :mod:`charger <charge.chargers>` module for full details.

The processing of the partial charges is the main difference between the chargers. There are three simple chargers that assign a single statistical value to the atoms:

    * The :class:`MeanCharger <charge.chargers.MeanCharger>` assigns the mean partial charge.
    * The :class:`MedianCharger <charge.chargers.MedianCharger>` assigns the median partial charge.
    * The :class:`ModeCharger <charge.chargers.ModeCharger>` assigns the mode of the partial charges. In case of a multimodal distribution, it chooses the mode closest to the median.

However, using the simple chargers might result in a total charge that is way off. For this reason, charge_assign was designed to solve the charge assignment problem: Given a histogram of partial charges for each atom, assign values such that the frequencies of the selected charges are maximal and the distance to a target total charge is less than a user-defined limit. charge_assign offers three redundant implementations that solve this problem, which differ in speed, but should always yield the same results. Sorted from fast to slow:

    * The :class:`CDPCharger <charge.chargers.CDPCharger>` solves the problem with a dynamic program implemented in C.
    * The :class:`ILPCharger <charge.chargers.ILPCharger>` solves the problem with an integer linear program using `PuLP <https://pythonhosted.org/PuLP/>`_.
    * The :class:`DPCharger <charge.chargers.DPCharger>` solves the problem with a dynamic program implemented in Python.

Using those Chargers will yield partial charges for which the total sum is close to the total charge of the molecule and for which the frequencies have been maximised. When looking at the molecules the assigned charges of two atoms with identical k-neighborhoods might differ as the chargers are not bound to apply the same charges. To get symmetric charges you can use these chargers which will (if possible) assign identical charges to atoms with identical atom neighborhoods:

    * The :class:`SymmetricDPCharger <charge.chargers.SymmetricDPCharger>` solves the problem with a dynamic program implemented in Python and assign symmetric charges.
    * The :class:`SymmetricILPCharger <charge.chargers.SymmetricILPCharger>` solves the problem with an integer linear program using `PuLP <https://pythonhosted.org/PuLP/>`_ and assign symmetric charges.

Batch Computations
------------------

charge_assign uses `nauty <http://users.cecs.anu.edu.au/~bdm/nauty/>`_ to hash an atom's neighborhood graph. If you are using multiple chargers, you might want to consider using one :class:`Nauty <charge.nauty.Nauty>` instance to reduce the number of processes::

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


Caching
-------

To speed up a large number of batch computations, charge_assign offers a caching mechanism. It can be easily activated by setting the caching option::

    repo = Repository.read('repo.zip')
    charger = CDPCharger(repository=repo, caching=True)

To reflect any changes in the repository in the cache, the versioning option needs to be set::

    repo = Repository.read('repo.zip', versioning=True)
    charger = CDPCharger(repository=repo, caching=True)

