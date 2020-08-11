.. vcregression documentation master file, created by
   sphinx-quickstart on Tue Jun 16 17:23:21 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for ``vcregression``
=================================================================

``vcregression`` is command line interface for statistical inference of sequence-function relationships using *Gaussian process(GP) regression*.

Basic functions of ``vc_regression`` include:

* estimation of the **variance components** from partially observed fitness landscapes
* calculating the **maximum a posterior estimate (MAP)** using *GP regression*
* calculating the **posterior variance**
* **posterior sampling** using *Hamiltonian Monte Carlo*

.. _installation:

Installation
-------------
``vc_regression`` is a command line interface written in Python 3. To install the software, simply clone the `repository <https://github.com/davidmccandlish/vcregression>`_ to your local directory. Before running the software, install all dependencies using the ``pip`` package manager with the following command line: ::

  pip install requirements.txt

.. _quickstart:

Quick Start
-----------

Estimation of the variance components
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a quick demonstration of the method using the sample data file ``smn1data.csv`` (Wong et al. 2018) [#wong2018]_, first execute the following command line to estimate the variance components ::

  python3 vc_prep.py 4 8 -data data/Smn1/smn1data.csv


MAP estimate
^^^^^^^^^^^^^

To calculate the maximum a posterior estimate (MAP) using the *lambdas* inferred using the command line above, execute the following command line::

  python3 vc_map_estimate.py 4 8 -data data/Smn1/smn1data.csv -lambdas out/lambdas.txt

.. The arguments are:
..
.. * `-a` size of the alphabet
..
.. * `-l` length of the sequence
..
.. * `-seqs` path to csv file of sequences in the training set
..
.. * `-y` path to vector of phenotypic measurements
..
.. * `-vars` path to vector of noise variances


Posterior variance
^^^^^^^^^^^^^^^^^^


Execute the following command to get the posterior variances of for a specified subset of sequences::

  python3 vc_variance.py 4 8 -seqs data/seqsSmn1.txt -lambdas out/lambdas_star.txt -vars data/varSmn1.txt -seqsvar data/seqspossample.txt





Contents
----------
.. toctree::
   :maxdepth: 1

   howtouse
   dataprep
   smn1
   implementation




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. [#wong2018] Wong et al. 2018. *Quantitative Activity Profile and Context Dependence of All Human 50 Splice Sites.*
