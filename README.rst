
==============================================================================
``vcregression``: variance component regresion for sequence-function relationships
==============================================================================

``vcregression`` is a command line interface for statistical inference of sequence-function relationships using *Gaussian process(GP) regression*. See the `vcregression documentation`_ for information about using ``vcregression`` and details about the recommended installation process.



Install dependencies
---------------------
```
pip install  -r requirements.txt
```


Quick start
------------

Estimating variance components::

  python3 vc_prep.py 4 8 -name smn1  -data data/Smn1/smn1data.csv -cv True

MAP estimate::

  python3 vc_map_estimate.py 4 8 -name smn1 -data data/Smn1/smn1data.csv -lambdas smn1/lambdas.txt

Analytical posterior variance::

  python3 vc_pos_var.py 4 8 -name smn1 -data data/Smn1/smn1data.csv -lambdas smn1/lambdas.txt -seqsvar data/Smn1/smn1seqpos.csv

HMC sampling::

  python3 vc_hmc.py 4 8 -name smn1 -data data/Smn1/smn1data.csv -lambdas smn1/lambdas.txt -MAP smn1/map.txt -step_size 1e-05 -n_steps 10 -n_samples 1000 -n_tunes 20 -starting_position 'random' -intermediate_output True -sample_name hmc1 -intermediate_output False
