from setuptools import setup, find_packages

setup(
        name='molli_dsi_ts',
        version=-'1.0',
        packages=find_packages(),
        install_requires=['numpy', 'pandas', 'rdkit'],
        python_requires = '>=3.9',
        )
