from setuptools import setup, find_packages

setup(
    name='nbotools',
    version='0.1.0',
    description='Misc. tools for NBO calculations (parsing, input generation, etc)',
    author='Karnjit Parmar',
    author_email='ksparmar@illinois.edu',
    packages=find_packages(),  # Automatically finds 'nbotools' package
    install_requires=['numpy', 'pandas'],       # Add any dependencies here, e.g., ['numpy']
    python_requires='>=3.6',
)
