.. _implementations:


Implementation
==============
``vc_prep.py``
---------------------------------------

This *CLI* estimates the prior **lambdas** which is needed for all subsequent analyses.

  .. argparse::
     :filename: vc_prep.py
     :func: parser
     :prog: vc_prep

Output files
>>>>>>>>>>>>
Running this command will output a .csv file ``lambdas.csv`` in the output directory. This file contains the optimal *lambdas*. Here is an example ``lambdas.csv`` file:

.. csv-table::
  :widths: 6, 8

  order_0,3356027.47
  order_1,153814.18
  order_2,44970.01
  order_3,3453.46
  order_4,363.17
  order_5,318.42
  order_6,68.85
  order_7,43.17
  order_8,29.17



This CLI also outputs an ``variance_component.txt`` file. Here is an example:

.. csv-table::
  :widths: 6, 8

  order_1,0.127
  order_2,0.390
  order_3,0.180
  order_4,0.071
  order_5,0.149
  order_6,0.048
  order_7,0.026
  order_8,0.006



``vc_map_estimate.py``
---------------------------------------
This *CLI* calculates the **maximum a posterior (MAP) estimate**

  .. argparse::
     :filename: vc_map_estimate.py
     :func: parser
     :prog: vc_map_estimate

Output file
>>>>>>>>>>>
Executing this command output the MAP estimate ``map.csv``. Here is the top 10 rows of an example:

.. csv-table::
  :widths: 10, 10

  sequence,phenotype
  AAAAAAAA,2.038
  AAAAAAAC,1.887
  AAAAAAAG,3.017
  AAAAAAAU,3.021
  AAAAAACA,4.438
  AAAAAACC,2.237
  AAAAAACG,4.232
  AAAAAACU,2.520
  AAAAAAGA,34.346
  AAAAAAGC,34.681


``vc_pos_var.py``
-----------------------------------------------------
This *CLI* calculates the **analytical posterior variance** for a list of sequences.


  .. argparse::
    :filename: vc_pos_var.py
    :func: parser
    :prog: vc_pos_var

input sequences
>>>>>>>>>>>>>>>
The user needs to specify the sequences for which to estimate the posterior variance. The sequences needs to be provided in a ``.csv`` file with one sequence per row. For example

.. csv-table::
  :widths: 10

  "AAAUGAUA"
  "CAUUCGUC"
  "UAGCGUCU"
  "GCGCGAUC"
  "AACCACGU"
  "CGACUCGA"
  "UCCCGUUU"
  "CAACGCAA"
  "CGCUAGGA"
  "AAAUCGAG"

Output file
>>>>>>>>>>>>
Executing this command will ouput a ``.csv`` file named ``varpos.txt``.

.. csv-table::
  :widths: 10, 6

  sequence,variance
  AAAUGAUA,0.082
  CAUUCGUC,0.576
  UAGCGUCU,0.023
  GCGCGAUC,108.691
  AACCACGU,1.879
  CGACUCGA,9.012
  UCCCGUUU,0.370
  CAACGCAA,0.290
  CGCUAGGA,0.087
  AAAUCGAG,88.807




``vc_hmc.py``
-----------------------------------------------------
This *CLI* is used to perform **posterior sampling** using the *Hamitonian Monte Carlo method*.

  .. argparse::
    :filename: vc_hmc.py
    :func: parser
    :prog: vc_hmc

Output files
>>>>>>>>>>>>

This command output 3 files:

Summary file
++++++++++++
This is a file with the suffix ``_hmc_summary.txt`` containing parameters for the HMC run.
Here is an example of a file's content:

.. csv-table::
  :widths: 10, 6

  initial_step_size,1e-05
  final_step_size,0.04301
  ntunes,100
  n_steps,100
  n_samples,200
  start_time,2020-07-28-16-27
  finish_time,2020-07-28-16-59
  total_time,0:32:00


Sample file
+++++++++++
This is a file with the suffix ``_hmc_sample.txt`` containing all HMC samples of the run. The number of columns is equal to the number of HMC samples specified by the user and the number of rows is equal to the total number of possible sequences.

Variance file
++++++++++++++
This is a file with the suffix ``_hmc_variances.txt`` containing the variance of all sequences calculated using the HMC samples.


``vc_hmc_diagnosis.py``
---------------------------------------

This CLI calculates the potential reduction factors for all sequences given multiple **hmc** samples output by ``vc_hmc.py``.

  .. argparse::
     :filename: vc_hmc_diagnosis.py
     :func: parser
     :prog: vc_hmc_diagnosis.py

This command outputs ``R_hat.txt`` which contains the potential scale reduction factors for all sequences.
