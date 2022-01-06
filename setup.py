from setuptools import setup

setup(
    name = "plotLattice",
    packages = ['plotLattice'],
    version = '0.0.1',
    author = 'Chenhuan Wang',
    author_email = 's6cnwang@uni-bonn.de',
    description = 'plotting script for data from CosmoLattice program, possibly also suitable for LatticeEasy',
    install_requires=["numpy", "matplotlib"]
)
