# charge_assign
[![DOI](https://zenodo.org/badge/103553886.svg)](https://zenodo.org/badge/latestdoi/103553886)

## Requirements

* [networkx](http://networkx.github.io/) ` >= 2.0`
* [msgpack-python](https://pypi.python.org/pypi/msgpack-python) ` >= 0.4.8`
* [numpy](http://www.numpy.org) ` >= 1.14.0`
* [nauty](http://users.cecs.anu.edu.au/~bdm/nauty/) ` >= 26r7` (This modules relies on the `dreadnaut` executable which is part of the `nauty` package.)

* optional: [rdkit](https://pypi.python.org/pypi/rdkit) ` >= v2017.03.3`

## Installation

Charge assign can be used directly as a Python library, and it has a server
mode in which it assigns charges based on a repository. For the server, see
below under Docker.

Installing the dependencies is easiest by using Anaconda, as it can install all
dependencies for you automatically. It is possible to use a combination of
virtualenv, the system package manager, and manual installation as well, but
this takes more work. Both approaches are detailed below.

### Conda

First, we create a new conda environment:
```bash
conda create -n charge_assign
```

Next, we can install charge_assign and its dependencies into it:
```bash
conda install -c rdkit,conda-forge -n charge_assign rdkit nauty
source activate charge_assign
cd charge_assign
pip install .
```

### Virtualenv

The rdkit and nauty packages are not available from PyPI, so they cannot be
installed in the usual way inside a virtual environment. Instead, they must be
installed manually first.

If you are running Ubuntu Linux, the following command will install rdkit from
the repositories:

```bash
sudo apt-get install python-rdkit librdkit1 rdkit-data
```

On Fedora, CentOS and RHEL, you can use

```bash
sudo yum install rdkit
```

Installation instructions for nauty are available on the [Nauty
homepage](http://pallini.di.uniroma1.it/). You will have to either add its
directory to your PATH variable, set the NAUTY_EXC environment variable to
the directory, or add the location of the dreadnaut program
to charge/settings.py.

Next, we can make a virtual environment and install charge_assign and its
dependencies.

```bash
virtualenv -p python3 </path/to/env>
. </path/to/env>/bin/activate
cd charge_assign
pip install .
```

## Docker server

The charge assign server is a simple WSGI server on top of the charge assign
library. It loads a repository ZIP file (see Repository.write()), and accepts
molecule definitions in LGF format via a REST API, which then returns an LGF
with charges included.

The easiest way to use the service is through Docker:

```bash
docker pull enitram/charge_assign
docker run -p 8080:8080 --mount-type=bind,source=/path/to/repo.zip,destination=/home/charge_assign/repo.zip --name charge_assign_server enitram/charge_assign
```

The server needs a repository in a ZIP file, which is bind-mounted in in the
above example. You can create such a repository using the
`scripts/build_repo.py` utility, see `scripts/build_repo.py -h` for
instructions.

Once the server is running, you can connect to it through HTTP:

```bash
curl --data-binary @molecule.lgf -H 'Content-Type: text/plain' 'http://localhost:8080/charge_assign?total_charge=0'
```

Note that the data must be submitted as content type text/plain, while curl by
default sends the input as a form, so the content type needs to be specified
explicitly here. The --data-binary option is needed (rather than --data) because
LGF is a whitespace-sensitive format, and --data does not preserve whitespace.
The REST API is trivial: just send the LGF to `/charge_assign?total_charge=<c>`
and you'll get a corresponding LGF with charges out, a 400 error if no charges
could be assigned for some reason, and a 400 error if the input was invalid.

## Tests

To run the tests, use

```bash
pip install -e .[dev]
```

to install the required dependencies, and then run

```bash
pytest --cov
```

to run the test suite.

## Documentation

To build the documentation, you need to install sphinx first. Installing sphinx is easiest by using Anaconda.

```
conda install -n charge_assign sphinx sphinx_rtd_theme
```

Then run sphinx.

```
source activate charge_assign
cd charge_assign/doc/
make html
```

You can then find the documentation in `charge_assign/doc/build/html/`.
