from setuptools import setup, find_packages

setup(
    name='substructure',
    version='0.1.0',
    description='additional tools to modify substructures',
    author='Karnjit Parmar',
    author_email='ksparmar@illinois.edu',
    packages=find_packages(),  
    install_requires=['numpy', 'pandas', 'rdkit'],       # Add any dependencies here, e.g., ['numpy']
    python_requires='>=3.7',
)
