aggregator
==========

Combining peptides to predict the relative abundance of proteins

=== Developing Package ===

Input Data:
  * ICabund [n,m] # absolute or relative IC
  * SN [n,m] # signal:noise on the absolute or relative abundance
  * Mapping matrix [n,k] # mapping of n features to k species
    * diagonal for metabolomics


Indicator: relative or absolute abund

Classes:
  * 