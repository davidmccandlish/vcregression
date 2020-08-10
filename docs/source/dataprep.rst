.. _dataprep:

Data preparation
==============
The input data should be a ``.csv`` file containing three columns with no header. The three columns are `sequences`, `phenotype`, and `variance`.

`sequence`
  The sequence column should contain strings of sequences of uniform length.

`phenotype`
  This column should contain the phenotype(fitness) of individual sequences.

`variance`
  This column should contain estimation of the variance of measurement noise for each sequence.


example
^^^^^^^^
See the example data files ``smn1data.csv`` for further formatting guidelines.
Here are the top ten lines of ``smn1data.csv``:

.. csv-table::
   :widths: 10, 6, 6

   "UACGUUGGG",0.407,0.076
   "GGCGCUCAC",1.62,0.407
   "AUGGUAUCG",9.752,1.276
   "ACAGCUGGC",0.151,0.028
   "ACAGUUGAU",0.193,0.036
   "GACGUCCUC",20.191,3.182
   "AGCGUCUUC",0.124,0.023
   "AAGGUUGUG",13.497,1.912
   "CUAGCGUUA",0.208,0.046
   "AAGGCAGGA",0.357,0.105
