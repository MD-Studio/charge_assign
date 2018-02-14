from setuptools import setup
setup(
        name = 'charge_assign',
        packages = ['charge'],
        version = '0.0.1',
        description = 'Assigns charges to atoms in molecules by template matching',
        author = 'Martin Engler',
        author_email = 'martin.engler@cwi.nl',
        url = 'https://github.com/enitram/charge_assign',
        download_url = 'https://github.com/enitram/charge_assign/archive/master.tar.gz',
        license = 'Apache License 2.0',
        python_requires='>=3.5, <4',
        install_requires=[
                'msgpack-python>=0.4.8',
                'networkx>=2.0,<3',
                'numpy>=1.14.0,<2'
            ],
        keywords = ['Molecules', 'Charges', 'Force field', 'Graph theory'],
        classifiers = [
            'Development Status :: 3 - Alpha',
            'License :: OSI Approved :: Apache Software License',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6'],
        )
