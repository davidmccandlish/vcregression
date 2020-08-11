How to use vc_regression
========================


``vcregression`` provides a suite of tools for Bayesian statistical inference of sequence-function relationships using *Gaussian process(GP) regression*. The most important idea behind ``vcregression`` is to infer a prior corresponding to the proportions of variance due to each order of genetic interaction (the **variance components**) based on the *empirical autocorrelation function* observed in the data. The **variance components** provide a covariance structure over the sequence space, which allows us make various Bayesian statistical inferences.

Unlike most regression methods used in modeling sequence-function relationships that only consider genetic interactions (*epistasis*) up to the *2nd* and occasionally the *3rd* order, ``vcregression`` is capable to model *epistasis* among arbitrary number of sites, which makes it a powerful method for predicting the phenotypes of unsampled sequences.

The 4 major functionalities of ``vc_regression`` include:

* estimation of the **variance components** from partially observed fitness landscapes
* calculating the **maximum a posterior estimate (MAP)** using *GP regression*
* calculating the **posterior variance**
* **posterior sampling** using *Hamiltonian Monte Carlo*

Each functionality listed above correspond to a separate command line interface. For details, see :ref:`implementations`.

Before using ``vcregressoin``, make sure the input datafile is in the correct format. See :ref:`dataprep` for guidelines.

For detailed implementations of ``vcregression``, see :ref:`smn1` for an application to a high throughput dataset of human *5' splice sites*.

Notes on the upper limit of the size sequence space
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``vcregression`` leverages the isotropic property of the covariance matrix defined by the **variance components**. All fast calculations of ``vcregression`` relies on representing the covariance matrix using  the **Laplacian** :math:`\mathbf{L}` of the associated *Hamming graph*. As a result, most of the computational and memory cost comes from constructing and storing :math:`\mathbf{L}`, which is defined for the *whole sequence space*. So currently the largest sequence space ``vcregression`` can accomodate is on the order of several million. Therefore, before using ``vcregression`` on your own dataset, first calculate the size of the sequence space:

.. math:: N = a^l,

and make sure :math:`N` is no larger than several million.

If the sequence space associated with your dataset is much larger, a method using the dense covariance matrix is possible (but with computational complexity :math:`\mathcal{O}(n^3)`, if you have :math:`n` data points). However, this functionality has not been implemented in ``vcregression``. If you are interested in this alternative, please contact Juannan Zhou (`Email: jzhou@cshl.edu <jzhou@cshl.edu>`_).
