import argparse
parser = argparse.ArgumentParser()
parser.add_argument("a", help="alphabet size", type=int)
parser.add_argument("l", help="sequence length", type=int)
parser.add_argument("-name", help="name of output folder")
parser.add_argument("-data", help="path to input data",  type=str, required=True)
parser.add_argument("-lambdas", help="path to lambdas",
                    type=str, required=True)


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


if not os.path.exists('out'):
    os.makedirs('out')


args = parser.parse_args()


a = args.a
l = args.l

if args.name == None:
    args.name = "my_project"

name = args.name
outdir = name


# QC
if a**l > 5000000:
    print("sequence space is to big!")
    exit()

vc.preliminary_preparation(a, l)

data = pd.read_csv(args.data, header=None)

#########
babel = ''
for i in range(len(data)):
    babel += data[0][i]

alphabet = set(babel)

AA2N = dict([(sorted(alphabet)[i], i) for i in range(len(alphabet))])
N2AA = {v: k for k, v in AA2N.items()}


def seqAA2num(seq):
    return [AA2N[seq[i]] for i in range(len(seq))]


def seqnum2AA(seq):
    seqAA = [N2AA[seq[i]] for i in range(len(seq))]
    return ''.join(seqAA)

seqs = [seqAA2num(data[0][i]) for i in range(len(data))]
tr = np.array([vc.seq2pos(seqs[i]) for i in range(len(seqs))])
###########
AA = list(AA2N.keys())
seqsAll = [''.join(seq) for seq in itertools.product(AA, repeat=l)]


if np.shape(seqs)[1] != l:
    print("seqs file dimension incompatible!")
    exit()


ys = np.array(data[1])
sig2s = np.array(data[2])

lda_star = pd.DataFrame(np.array(pd.read_csv(
    args.lambdas, header=None, index_col=0)))
lda_star = np.array(lda_star).flatten()
vc.lda_star = lda_star

print("using lambdas = ", str(lda_star))

vc.set_data_as_global_parameters(seqs, ys, sig2s)
vc.construct_A_sparse()
vc.construct_E_sparse()


time_i = time.perf_counter()

print('Estimating posterior mean')
fstar = vc.compute_posterior_mean()
print('%.2f sec' % (time.perf_counter() - time_i))


pd.DataFrame({'seq': seqsAll, 'y': fstar}).to_csv(
    outdir + '/map.txt', index=False, header=['sequence', 'phenotype'])

print("Done!")
