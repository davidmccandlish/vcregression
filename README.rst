
==============================================================================
``vcregression``: variance component regression for sequence-function relationships
==============================================================================
.. image:: https://readthedocs.org/projects/vcregression/badge/?version=latest
    :target: https://vcregression.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

``vcregression`` is a command line interface for statistical inference of sequence-function relationships using *Gaussian process(GP) regression*. See the `vcregression documentation <https://vcregression.readthedocs.io/en/latest/>`_ for information about using ``vcregression``.

Installation
-------------
``vc_regression`` is a command line interface written in Python 3. To install the software, simply clone the `repository <https://github.com/davidmccandlish/vcregression>`_ to your local directory. Before running the software, install all dependencies using the ``pip`` package manager with the following command line: ::

  pip install -r requirements.txt




Quick start
------------

Estimating **variance components**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

  python3 vc_prep.py 4 8 -name smn1  -data data/Smn1/smn1data.csv -cv True

Calculate the **maximum a posterior (MAP)** estimate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

  python3 vc_map_estimate.py 4 8 -name smn1 -data data/Smn1/smn1data.csv -lambdas smn1/lambdas.txt

Calculate **posterior variance** for a list of sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

  python3 vc_pos_var.py 4 8 -name smn1 -data data/Smn1/smn1data.csv -lambdas smn1/lambdas.txt -seqsvar data/Smn1/smn1seqpos.csv

Draw 200 samples from the **posterior distribution** using *Hamiltonian Monte Carlo (HMC)*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

  python3 vc_hmc.py 4 8 -name smn1 -data data/Smn1/smn1data.csv -lambdas smn1/lambdas.txt -MAP smn1/map.txt -step_size 1e-05 -n_steps 10 -n_samples 200 -n_tunes 20 -starting_position random -intermediate_output True -sample_name hmc1 -intermediate_output False
