# charge_assign

## Requirements

* [networkx](http://networkx.github.io/) ` >= 2.0`
* [msgpack-python](https://pypi.python.org/pypi/msgpack-python) ` >= 0.4.8`
* [numpy](http://www.numpy.org) ` >= 1.14.0`
* [nauty](http://users.cecs.anu.edu.au/~bdm/nauty/) ` >= 26r7` (This modules relies on the `dreadnaut` executable which is part of the `nauty` package.)

* optional: [rdkit](https://pypi.python.org/pypi/rdkit) ` >= v2017.03.3`

## Installation

First, make a virtual environment and activate it, then move into the
charge_assign directory and use

```bash
pip install .
```

This should install most dependencies and charge_assign. The rdkit and nauty
packages are not availably from PyPI, so you will have to install those
manually. We are working on making charge_assign installable in Anaconda, where
its dependencies can be installed along with it automatically.
