import argparse
parser = argparse.ArgumentParser()
parser.add_argument("a", help="alphabet size", type=int)
parser.add_argument("l", help="sequence length", type=int)
parser.add_argument("-name", help="name of output folder")
parser.add_argument("-data", help="path to input data",
                    type=str, required=True)
parser.add_argument("-lambdas", help="path to lambdas",
                    type=str, required=True)
parser.add_argument("-MAP", help="path to MAP estimate", dest="MAP")
parser.add_argument("-step_size", help="initial leapfrog stepsize",
                    dest="step_size", type=float, default=1e-6)
parser.add_argument("-n_steps", help="number of leapfrog steps per iteration",
                    dest="n_steps", type=int, required=True)
parser.add_argument("-n_samples", help="number of samples to draw from the posterior",
                    dest="n_samples", type=int, required=True)
parser.add_argument("-n_tunes", help="number of HMC steps used for tuning the step_size parameter",
                    dest="n_tunes", type=int, default=100)
parser.add_argument("-starting_position",  help="starting position",
                    choices=['mode', 'random', 'custom'], dest="qstart", default='mode')
parser.add_argument("-starting_position_path",
                    help="path to starting position vector if using custom starting position", dest="qpath")
parser.add_argument("-intermediate_output", help="output intemediate samples and variance",
                    default=False)
parser.add_argument("-sample_name", help="name of the hmc sample")


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
from scipy.sparse import csr_matrix, dia_matrix, save_npz, load_npz
from scipy.optimize import minimize
from scipy.special import comb
from scipy.spatial.distance import hamming
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import cg
import arrow


import vc_regression as vc


def str2bool(v):
    return v.lower() in ("True", "true", "t", "1")

start_time = arrow.now().format('YYYY-MM-DD-HH-mm')
vc.start_time = start_time


args = parser.parse_args()
a = args.a
l = args.l


args.intermediate_output = str2bool(args.intermediate_output)

if args.name == None:
    args.name = "my_project"

name = args.name
outdir = name

outpath = outdir
if not os.path.exists(outpath):
    os.makedirs(outpath)

vc.outpath = outpath


if args.sample_name == None:
    args.sample_name = start_time
vc.sample_name = args.sample_name


if args.qstart == 'mode':
    q_start = np.zeros(a**l)
else:
    if args.qstart == 'random':
        q_start = np.random.normal(loc=0, scale=1, size=a**l)
    else:
        q_start = np.array(pd.read_csv(args.qpath, sep='\t', header=None)).T
        q_start = np.ndarray.flatten(q_start)


if a**l > 5000000:
    print("sequence space is too big!")
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
vc.seqsAll = seqsAll

if np.shape(seqs)[1] != l:
    print("seqs file dimension incompatible!")
    exit()


ys = np.array(data[1])
sig2s = np.array(data[2])

lda_star = pd.DataFrame(np.array(pd.read_csv(
    args.lambdas, header=None, index_col=0)))
lda_star = np.array(lda_star).flatten()
vc.lda_star = lda_star

MAP = np.array(pd.read_csv(args.MAP)['phenotype'])


vc.set_data_as_global_parameters(seqs, np.zeros(len(seqs)), sig2s)


##

vc.prepare_pos_sampling()


samples, varsample = vc.hamiltonian_monte_carlo(
    n_samples=args.n_samples,
    m=1,
    initial_position=q_start,
    tune=args.n_tunes,
    initial_step_size=args.step_size,
    num_steps=args.n_steps,
    max_energy_change=1000.0,
    intermediate_output=args.intermediate_output

)

samples = MAP + samples

finish_time = arrow.now().format('YYYY-MM-DD-HH-mm')

pd.DataFrame(samples.T, index=seqsAll).to_csv(outpath + "/" +
                                              args.sample_name + "_hmc_samples.txt", index=seqsAll, header=False)

pd.DataFrame(varsample.T, index=seqsAll).to_csv(outpath + "/" +
                                                args.sample_name + "_hmc_variance.txt", index=seqsAll, header=False)

# pd.DataFrame(lda_star).to_csv(outpath + "/" + args.sample_name + "_hmc_input_lambdas.txt", index = None, header = False)


summary = pd.DataFrame([args.sample_name, args.step_size, vc.step_size, vc.acceptance_ratio, args.n_tunes,
                        args.n_steps, args.n_samples, start_time, finish_time,
                        arrow.get(finish_time, 'YYYY-MM-DD-HH-mm') -
                        arrow.get(start_time, 'YYYY-MM-DD-HH-mm'), lda_star], index=["sample_name", "initial_step_size", "final_step_size", "final_acceptance_ratio",
                                                                                     "ntunes", "n_steps", "n_samples", "start_time",
                                                                                     "finish_time", "total_time", "input_lambdas"])

print(summary)

pd.DataFrame(summary).to_csv(outpath + "/" + args.sample_name +
                             "_hmc_summary.txt", index=True, header=False)
