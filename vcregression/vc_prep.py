import argparse

# add arguments
parser = argparse.ArgumentParser()
parser.add_argument("a", help="alphabet size", type=int)
parser.add_argument("l", help="sequence length", type=int)
parser.add_argument("-name", help="name of output folder")
parser.add_argument("-data", help="path to input data",
                    type=str, required=True)
parser.add_argument(
    "-cv", help="estimate lambdas using regularization \
    with regularization parameter chosen with 10-fold crossvalidation", default=True)


import numpy as np
import scipy as sp
import itertools
import sys
import time
import scipy as sp
import itertools
import os
import math
import csv
import pandas as pd
import random as rd
import statistics
from scipy.sparse import csr_matrix, dia_matrix
from scipy.optimize import minimize
from scipy.special import comb
from scipy.spatial.distance import hamming
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import cg
import vc_regression as vc


def str2bool(v):
    return v in ("True", "true", "t", "1")


args = parser.parse_args()
args.cv = str2bool(args.cv)


a = args.a
l = args.l


if args.name == None:
    args.name = "my_project"


outdir = args.name


if not os.path.exists(outdir):
    os.makedirs(outdir)


# Check if sequence space is too big
if a**l > 5000000:
    print("sequence space is to big!")
    exit()

vc.preliminary_preparation(a, l)

data = pd.read_csv(args.data, header=None)

babel = ''
for i in range(len(data)):
    babel += data[0][i]

alphabet = set(babel)

AA2N = dict([(sorted(alphabet)[i], i) for i in range(len(alphabet))])
N2AA = {v: k for k, v in AA2N.items()}
AA = list(AA2N.keys())
seqsAll = [''.join(seq) for seq in itertools.product(AA, repeat=l)]
pd.DataFrame(seqsAll).to_csv(
    outdir + "/sequences.txt", header=None, index=None)


def seqAA2num(seq):
    return [AA2N[seq[i]] for i in range(len(seq))]
####

seqs = [seqAA2num(data[0][i]) for i in range(len(data))]
tr = np.array([vc.seq2pos(seqs[i]) for i in range(len(seqs))])


if np.shape(seqs)[1] != l:
    print("seqs file dimension incompatible!")
    exit()


ys = np.array(data[1])
sig2s = np.array(data[2])


vc.initialize_computation(seqs, ys, sig2s)

all_distance_class_Q = all(map(lambda x: x > 0, vc.N_d))

rhod = vc.rho_d.copy()
rhod[0] -= np.mean(sig2s)
lambdas_naive = sp.linalg.inv(vc.W_kd.T).dot(rhod)
lambdas_naive_positive_Q = all(map(lambda x: x > 0, lambdas_naive))


if args.cv is True:
    print("estimating lambdas with regularization (regularization \
        parameter chosen using 10-fold crossvalidation)...")
    cv = True

elif not all_distance_class_Q:
    print("certain distance classes missing from data, estimating \
        lambdas using regularization (regularization parameter \
        chosen with 10-fold crossvalidation)...")
    cv = True
elif lambdas_naive_positive_Q:
    print("estimating lambdas using least squares")
    cv = False
else:
    print("naive lambdas contain nonpositive values, estimating \
        lambdas using regularization (regularization parameter \
        chosen with 10-fold crossvalidation)...")
    cv = True


betas = 10 ** np.arange(-2, 6, .5)

rownames = ["order_" + str(k) for k in range(l + 1)]


# Estimate lambdas using 10 fold crossvalidation
if cv is True:

    out = vc.lambdaCV(ys, tr, sig2s, betas, 10)
    beta_cv = out[1]
    lda = vc.solve_lambda_single_beta(ys, tr, sig2s, beta_cv)

    print("lambdas = ", str(lda))

else:

    lda = lambdas_naive
    print("lambdas = ", str(lda))

pd.DataFrame(lda, index=rownames).to_csv(outdir + "/lambdas.txt", header=None)


mks = [comb(l, k) * (a - 1)**k for k in range(l + 1)]

variance_components = np.array([lda[k] * mks[k] for k in range(1, l + 1)])
variance_components /= np.sum(variance_components)
print("variance components = ", str(variance_components))
pd.DataFrame(variance_components, index=rownames[1:]).to_csv(
    outdir + "/variance_components.txt", header=None)
