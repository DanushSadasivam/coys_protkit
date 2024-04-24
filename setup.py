from setuptools import setup

setup(
    name="coys",
    version="0.0.1",
    description="Protein toolkit package",
    maintainer="Danush Sadasivam",
    maintainer_email="dsadasiv@andrew.cmu.edu",
    license="MIT",
    packages=["coys"],
    long_description="""A module containing functions 
      for calculating properties of peptide sequences, summarizing sequence properties, 
      and generating sequences using genetic algorithms. 
      The `calculate_properties()` function computes various physicochemical 
      properties of peptide sequences. The `sequence_summary` provides a summary of these 
      properties and visualizes the amino acid distribution. 
      The `generate_seq` function utilizes a genetic algorithm to generate peptide sequences 
      optimized for user-defined properties. This module relies on external libraries such 
      ModlAMP for obtaining the individual protein properties, and pyGAD for implementing the genetic algorithm.""",
)
