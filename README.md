# variance_component_regression
vcregression is a command line interface for statistical inference of sequence-function relationships using *Gaussian process(GP) regression*.

##Install all dependencies:
```
pip install requirements.txt

```


##Quick start
###Estimating variance components
```
python3 vc_prep.py 4 8 -data data/Smn1/smn1data.csv

```
###MAP estimate
```
python3 vc_map_estimate.py 4 8 -data data/Smn1/smn1data.csv -lambdas out/lambdas.txt

```
