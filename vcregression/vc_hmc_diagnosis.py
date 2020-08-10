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
import matplotlib.pyplot as plt
import vc_regression as vc
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-samples", help="tuple of paths to HMC samples",
                    dest="samplepaths", nargs="+", type=str, required=True)
parser.add_argument("-name", help="project name")

args = parser.parse_args()
samplepaths = args.samplepaths


paths = samplepaths[0].split(',')
n = len(paths)

if args.name == None:
    args.name = "my_project"

name = args.name
outdir = name

vc.outpath = outdir

allsamples = []

for i in range(n):
    dat = np.array(pd.read_csv(paths[i], header=None))
    seqsAll = dat[:, 0]
    dat = dat[:, 1:]
    allsamples.append(dat)

vc.seqsAll = seqsAll

G, chainlength = allsamples[0].shape

allsamples = np.array(allsamples)


initial_phis = np.zeros([G, n])
for i in range(n):
    initial_phis[:, i] = allsamples[i, :, 0]

pot_scale_red_facs = vc.compute_pot_scale_red_facs(
    G, allsamples, [0, chainlength])

pd.DataFrame(pot_scale_red_facs, index=seqsAll).to_csv(
    outdir + "/R_hat.txt", header=None)


print('potential scale reduction factor:')
print('%.3f (%s) ~ %.3f (%s)' %
      (pot_scale_red_facs.min(), np.array(seqsAll)[pot_scale_red_facs.argmin()], pot_scale_red_facs.max(), np.array(seqsAll)[pot_scale_red_facs.argmax()]))


from_sample, to_sample = 1, chainlength
colors = ['blue', 'orange', 'green', 'red', 'purple',
          'pink', 'brown', 'violet', 'maroon', 'gold']
xlimits, ylimits = None, None

# # ----------

j = pot_scale_red_facs.argmax()
vc.plot_trajectory([j, j], initial_phis, allsamples, [
                   from_sample, to_sample], colors, xlimits, ylimits)
