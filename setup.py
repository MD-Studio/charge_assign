from setuptools import setup, Extension

dp_module = Extension('_dp',
        sources=['charge/c/dp_wrap.c', 'charge/c/dp.c'])

setup(
        name = 'charge_assign',
        packages = ['charge', 'charge_server'],
        version = '0.0.1',
        description = 'Assigns charges to atoms in molecules by template matching',
        author = 'Martin Engler',
        author_email = 'martin.engler@cwi.nl',
        url = 'https://github.com/enitram/charge_assign',
        download_url = 'https://github.com/enitram/charge_assign/archive/master.tar.gz',
        license = 'Apache License 2.0',
        python_requires='>=3.5, <4',
        install_requires=[
                'connexion',
                'flask',
                'msgpack-python>=0.4.8',
                'networkx==2.0',
                'numpy>=1.14.0,<2',
                'pulp>=1.6.8'
            ],
        ext_modules = [dp_module],
        scripts = ['scripts/build_repo.py'],
        extras_require={
            'dev': [
                'flask_testing',
                'pytest',
                'pytest-pep8',
                'pytest-cov',
                'pytest-xdist'
            ]
        },
        keywords = ['Molecules', 'Charges', 'Force field', 'Graph theory'],
        classifiers = [
            'Development Status :: 3 - Alpha',
            'License :: OSI Approved :: Apache Software License',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6'],
        )
