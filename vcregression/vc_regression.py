import numpy as np
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
import scipy.stats as st
from tqdm import tqdm

# basic functions


def inv(v):
    """return a vector a reciprocals for v"""

    return([1 / v[i] for i in range(len(v))])


def gd(x, y):
    """calculate reciprocal of variable, if variable=0, return 0"""
    if y == 0:
        return 0
    else:
        return x / y

# a,l - dependent functions


def seq2pos(seq):
    """return lexicographical position of a numerical sequence"""
    nums = []
    for i in range(len(seq)):
        nums.append(seq[i] * (a**(len(seq) - i - 1)))
    return(sum(nums))


def sequence_to_position(seq):
    """return lexicographical position of a numerical sequence using the
     seq_to_pos_converter"""

    return sp.sum(seq * seq_to_pos_converter)


def w(k, d):
    """return value of the Krwatchouk polynomial for k, d"""
    s = 0
    for q in range(l + 1):
        s += (-1)**q * (a - 1)**(k - q) * comb(d, q) * comb(l - d, k - q)
    return 1 / a**l * s


def W_kd_mat():
    """return full matrix l+1 by l+1 Krawtchouk matrix"""

    # Construct W_kd
    global W_kd
    W_kd = np.zeros([l + 1, l + 1])
    for k in range(l + 1):
        for d in range(l + 1):
            W_kd[k, d] = w(k, d)
    return W_kd


def construct_L_sparse():
    """construct graph Laplacian for using the global a and l"""

    # Generate bases and sequences
    bases = sp.array(range(a))
    seqs = list(itertools.product(bases, repeat=l))

    # Find indices of L at which L = -1
    row_ids, col_ids, values = [], [], []
    for i in range(G):
        row_ids.append(i)
        col_ids.append(i)
        values.append(l * (a - 1))
        for site in range(l):
            for base in bases:
                seq_i = sp.array(seqs[i])
                if base != seq_i[site]:
                    seq_i[site] = base
                    j = sequence_to_position(seq_i)
                    row_ids.append(i)
                    col_ids.append(j)
                    values.append(-1)

    # Save L_sparse as sparse matrix
    L_sparse = csr_matrix((values, (row_ids, col_ids)), shape=(G, G))

    # Return
    return L_sparse


def prepare_L_sparse(from_dir='sparse_matrix/L/'):
    """If L given the a and l is present in the sparse_matrix folder
    ,load the matrix. Otherwise, the matrix is constructed"""

    # Set global parameters for later use
    global L_sparse

    # Get list of current sparse matrices
    spm_list = os.listdir(from_dir)

    # If the matrix desired has been made already, load it. Otherwise,
    # construct and save it
    file_name = 'L_sparse_a' + str(a) + '_' + 'l' + str(l) + '.npz'

    if file_name in spm_list:

        print('Loading L_sparse ...')
        start_time = time.perf_counter()
        L_sparse = load_npz(from_dir + file_name)
        print('%.2f sec' % (time.perf_counter() - start_time))

    else:

        print('Constructing L_sparse ...')
        start_time = time.perf_counter()
        L_sparse = construct_L_sparse()
        print('%.2f sec' % (time.perf_counter() - start_time))

        save_npz(from_dir + file_name, L_sparse)


def construct_MAT():
    """Construct entries of powers of L. 
    Column: powers of L. 
    Row: Hamming distance"""

    # Set global parameters for later use
    global MAT, MAT_inv

    # Construct C
    C = np.zeros([l + 1, l + 1])
    for i in range(l + 1):
        for j in range(l + 1):
            if i == j:
                C[i, j] = i * (a - 2)
            if i == j + 1:
                C[i, j] = i
            if i == j - 1:
                C[i, j] = (l - j + 1) * (a - 1)

    # Construct D
    D = sp.array(np.diag(l * (a - 1) * np.ones(l + 1), 0))

    # Construct B
    B = D - C

    # Construct u
    u = np.zeros(l + 1)
    u[0], u[1] = l * (a - 1), -1

    # Construct MAT column by column
    MAT = np.zeros([l + 1, l + 1])
    MAT[0, 0] = 1
    for j in range(1, l + 1):
        MAT[:, j] = sp.array(sp.mat(B)**(j - 1) * sp.mat(u).T).ravel()

    # Invert MAT
    MAT_inv = np.linalg.inv(MAT)

# coefficients of d-adajcent matrix in L**k


def solve_c_k(d):
    """return the coefficients of the adjacency matrix at Hamming
    distance d as a polynomial of L"""

    # Set global parameters for later use
    global c_k

    # Solve for c_k: coefficients in L**k for the matrix with (i,j)=1 if
    # d(i,j)=d
    c_k = MAT_inv[:, d]


# data dependent functions

def construct_A_sparse():
    """construct a sparse matrix that maps the m-length vector to a
    a**l vector"""

    # Set global parameters for later use
    global A_sparse

    # Convert each sequence to its corresponding position
    poss = sp.array(sp.mat(seqs) * sp.mat(seq_to_pos_converter).T).ravel()

    # Construct A with positions of the input sequences
    row_ids, col_ids, values = poss, np.arange(num_data), np.ones(num_data)

    # Save A as sparse matrix
    A_sparse = csr_matrix((values, (row_ids, col_ids)), shape=(G, num_data))


def construct_E_sparse():
    """construct the diagonal matrix containing noise variance"""

    # Set global parameters for later use
    global E_sparse

    # Construct E and save it as sparse matrix
    E_sparse = dia_matrix((sig2s, sp.array([0])), shape=(num_data, num_data))


# Sparse matrix operators

# multiply by L
def L_opt(phi):
    """multiply L by the vector phi"""
    return L_sparse.dot(phi)


def M_opt(v, lambdas):
    """
    return the value M.v

    Keyward arguments:

    v -- a vector of length a**l

    lambdas -- l + 1 vector of nonnegative values that specify the 
    eigenvalues of a matrix M
    """

    b_k = sp.array(sp.mat(
        MAT_inv) * sp.mat(sp.array(sp.mat(W_kd).T * sp.mat(lambdas).T).ravel()).T).ravel()
    max_power = len(b_k) - 1
    Lsv = np.zeros([G, len(b_k)])
    Lsv[:, 0] = b_k[0] * v
    power = 1
    while power <= max_power:
        v = L_opt(v)
        Lsv[:, power] = b_k[power] * v
        power += 1
    Kv = Lsv.sum(axis=1)
    return Kv


def R_d_opt(v):
    """multiply by d-adjacenty matrix"""
    max_power = len(c_k) - 1
    Lsv = np.zeros([G, len(c_k)])
    Lsv[:, 0] = c_k[0] * v
    power = 1
    while power <= max_power:
        v = L_opt(v)
        Lsv[:, power] = c_k[power] * v
        power += 1
    Rdv = Lsv.sum(axis=1)
    return Rdv


def K_opt(v):
    """multiply by the whole covariance matrix. b_k is set by lda_star"""

    max_power = len(b_k) - 1
    Lsv = np.zeros([G, len(b_k)])
    Lsv[:, 0] = b_k[0] * v
    power = 1
    while power <= max_power:
        v = L_opt(v)
        Lsv[:, power] = b_k[power] * v
        power += 1
    Kv = Lsv.sum(axis=1)
    return Kv


def K_BB_E(v):
    """ multiply by the m by m matrix K_BB + E"""
    return A_sparse.T.dot(K_opt(A_sparse.dot(v))) + E_sparse.dot(v)


# basic inference of rhod and lambda


def compute_rhod_and_Nd():
    """compute autocorrelation function rho_d and number of sequence 
    pairs for all distance classes, N_d"""

    Ays = A_sparse.dot(ys)
    A1s = A_sparse.dot(np.ones(len(ys)))

    # Construct MAT and its inverse
    construct_MAT()

    # Compute rho_d and N_d
    rho_d, N_d = np.zeros(l + 1), np.zeros(l + 1)
    for d in range(l + 1):
        solve_c_k(d)
        rho_d[d] = sp.sum(Ays * R_d_opt(Ays))
        N_d[d] = sp.sum(A1s * R_d_opt(A1s))

    rho_d = [gd(rho_d[i], N_d[i]) for i in range(len(rho_d))]

    # Return
    return rho_d, N_d


def Frob(lda, M, a):
    """calculate the cost function given lambdas and a"""
    Frob1 = (sp.mat(lda) * sp.mat(M) * sp.mat(lda).T)[0, 0]
    Frob2 = 2 * sp.sum(lda * a)
    return Frob1 - Frob2


def grad_Frob(lda, M, a):
    """gradient of the function Frob(lda, M, a)"""
    grad_Frob1 = 2 * sp.array(sp.mat(M) * sp.mat(lda).T).ravel()
    grad_Frob2 = 2 * a
    return grad_Frob1 - grad_Frob2


def solve_for_lambda(rho_d, N_d):
    """solve for lambdas using least squares with the given rho_d 
    and N_d"""

    # Construct rho_d_prime
    rho_d_prime = rho_d.copy()
    # rho_d_prime[0] -= sp.array(sig2s).mean()

    # Construct M
    M = np.zeros([l + 1, l + 1])
    for i in range(l + 1):
        for j in range(l + 1):
            for d in range(l + 1):
                M[i, j] += N_d[d] * w(i, d) * w(j, d)

    # Construct a
    a = np.zeros(l + 1)
    for i in range(l + 1):
        for d in range(l + 1):
            a[i] += N_d[d] * w(i, d) * rho_d_prime[d]

    # Minimize the objective function Frob with lda > 0
    lda_initial = sp.array(np.linalg.inv(M) * sp.mat(a).T).ravel()
    res = minimize(fun=Frob, jac=grad_Frob, args=(
        M, a), x0=lda_initial, method='L-BFGS-B', bounds=[(0, None)] * (l + 1))
    if not res.success:
        print(res.message)
        print()
    lda_star = res.x

    # Return
    return lda_star


# pass variables to vc module

def set_global_parameters(a0, l0):
    """Set global parameters for later use"""
    global a, l, G, seq_to_pos_converter

    a = a0
    l = l0
    G = a**l
    seq_to_pos_converter = np.flip(a**sp.array(range(l)), axis=0)


def preliminary_preparation(a, l):
    """preparing all data-indepdent objects given a and l"""

    # Pass a,l
    set_global_parameters(a, l)

    # Construct L_sparse
    if not os.path.exists('sparse_matrix/L'):
        os.makedirs('sparse_matrix/L')

    time_i = time.perf_counter()
    prepare_L_sparse()
    # print('%.2f sec' % (time.perf_counter() - time_i))

    # Construct entries for L^k
    construct_MAT()

    # Construct matrix W_kd containing krawtchouk polynomials
    W_kd_mat()

    # Construct second order difference matrix for regularization
    # Diff2 = np.zeros((l - 1, l + 1))

    # for i in range(Diff2.shape[0]):
    #     Diff2[i, i:i + 3] = [-1, 2, -1]

    Diff2 = np.zeros((l - 2, l))

    for i in range(Diff2.shape[0]):
        Diff2[i, i:i + 3] = [-1, 2, -1]

    global Rmat

    Rmat = Diff2.T.dot(Diff2)


def set_data_as_global_parameters(seqs0, ys0, sig2s0):
    """Set global parameters for later use"""
    global num_data, seqs, ys, sig2s

    num_data = len(seqs0)
    seqs = seqs0
    ys = ys0
    sig2s = sig2s0


def initialize_computation(seqs, ys, sig2s):
    """calculate all data-dependent objects"""

    # Set global parameters for later use
    global lda_star, rho_d, N_d

    # Set data as global parameters
    set_data_as_global_parameters(seqs, ys, sig2s)

    # Construct A_sparse
    construct_A_sparse()

    # Construct E_sparse
    construct_E_sparse()

    # Compute rho_d and N_d
    print('Computing rho_d & N_d ...')
    time_i = time.perf_counter()
    rho_d, N_d = compute_rhod_and_Nd()
    print('%.2f sec' % (time.perf_counter() - time_i))

    # Solve for lambda
    print('Solving for lambda ...')
    time_i = time.perf_counter()
    lda_star = solve_for_lambda(rho_d, N_d)
    print('%.2f sec' % (time.perf_counter() - time_i))


# Crossvalidation lambda utilities

def construct_A_sparse_custom(poss):
    """construct A_sparse using the list of positions poss"""

    num_data = len(poss)

    # Construct A with positions of the input sequences
    row_ids, col_ids, values = poss, np.arange(num_data), np.ones(num_data)

    # Save A as sparse matrix
    A_sparse = csr_matrix((values, (row_ids, col_ids)), shape=(G, num_data))

    return A_sparse


def compute_rhod_and_Nd_custom(val, list):
    """calculate rho_d and N_d with specified list of values 
    val and list"""

    A_sparse_custom = construct_A_sparse_custom(list)

    # Compute rho_d and N_d
    rho_d, N_d = np.zeros(l + 1), np.zeros(l + 1)
    for d in range(l + 1):
        solve_c_k(d)
        Ays = A_sparse_custom.dot(val)
        A1s = A_sparse_custom.dot(np.ones(len(val)))
        rho_d[d] = sp.sum(Ays * R_d_opt(Ays))
        N_d[d] = sp.sum(A1s * R_d_opt(A1s))
    rho_d = [gd(rho_d[i], N_d[i]) for i in range(len(rho_d))]

    # Return
    return rho_d, N_d


def Frob_reg(theta, M, a, beta):
    """cost function for regularized least square method for inferring 
    lambdas"""
    Frob1 = sp.exp(theta).dot(M).dot(sp.exp(theta))
    Frob2 = 2 * sp.exp(theta).dot(a)
    # return Frob1 - Frob2 + beta * theta. dot(Rmat).dot(theta)
    return Frob1 - Frob2 + beta * theta[1:]. dot(Rmat).dot(theta[1:])


def solve_lambda(ys, list, vars, betas):
    """solve for lambdas using regularized least squares for a list of 
    regularization parameters betas

    Keyword arguments:
    ys -- vector of values
    list -- list of positions for values in ys
    vars -- noise variance for values in ys
    betas -- list of regulariztion parameter
    """

    # Construct rho_d_prime
    rho_d_prime, N_d = compute_rhod_and_Nd_custom(ys, list)
    rho_d_prime[0] -= sp.array(vars).mean()

    M = np.zeros([l + 1, l + 1])
    for i in range(l + 1):
        for j in range(l + 1):
            for d in range(l + 1):
                M[i, j] += N_d[d] * w(i, d) * w(j, d)

    # Construct a
    a = np.zeros(l + 1)
    for i in range(l + 1):
        for d in range(l + 1):
            a[i] += N_d[d] * w(i, d) * rho_d_prime[d]

    ldas = []
    thetas = []

    for i in range(len(betas)):

        beta = betas[i]

        res = minimize(fun=Frob_reg, args=(M, a, beta), x0=np.zeros(
            l + 1), method='Powell', options={'xtol': 1e-8, 'ftol': 1e-8})

        thetas.append(np.array(res.x))
        ldas.append(np.exp(res.x))

    return thetas, ldas


def solve_lambda_single_beta(ys, list, vars, beta):
    """solve for lambdas using regularized least squares with 
    regularization parameter provided by beta"""

    # Construct rho_d_prime
    rho_d_prime, N_d = compute_rhod_and_Nd_custom(ys, list)
    rho_d_prime[0] -= sp.array(vars).mean()

    M = np.zeros([l + 1, l + 1])
    for i in range(l + 1):
        for j in range(l + 1):
            for d in range(l + 1):
                M[i, j] += N_d[d] * w(i, d) * w(j, d)

    # Construct a
    a = np.zeros(l + 1)
    for i in range(l + 1):
        for d in range(l + 1):
            a[i] += N_d[d] * w(i, d) * rho_d_prime[d]

    res = minimize(fun=Frob_reg, args=(M, a, beta), x0=np.zeros(
        l + 1), method='Powell', options={'xtol': 1e-8, 'ftol': 1e-8})

    # Return
    return np.exp(res.x)


def lambdaCV(val, tr, var, betas, nfolds=10):
    """
    Solve for optimal lambdas using regularized least squares with 
    regularization parameter chosen using crossvalidation

    """
    m = len(betas)
    order = list(range(len(tr)))
    rd.shuffle(order)
    folds = np.array_split(np.array(order), nfolds)

    rhods = np.zeros((10, l + 1))
    Nds = np.zeros((10, l + 1))
    for i in range(nfolds):
        rhods[i], Nds[i] = compute_rhod_and_Nd_custom(
            val[folds[i]], tr[folds[i]])

    mses = np.zeros((nfolds, m))

    for k in tqdm(range(nfolds)):
        sub = list(set(range(len(tr))).difference(set(folds[k])))
        ldas = solve_lambda(val[sub], tr[sub], var[sub], betas)[1]
        rhods[k][0] -= np.mean(var[sub])
        mses[k] = [np.sum((W_kd.T.dot(ldas[i]) - rhods[k])**2 *
                          (1 / np.sum(Nds[k])) * Nds[k]) for i in range(len(ldas))]

    betaopt = betas[np.argmin(np.mean(mses, axis=0))]

    return mses, betaopt


# Estimate MAP
def compute_posterior_mean(target_seqs=None):
    """compute the MAP"""

    # Set b_k = coeffs of K in L**k
    global b_k

    # Solve b_k
    W_kd = W_kd_mat()
    rho = sp.array(sp.mat(W_kd).T * sp.mat(lda_star).T).ravel()
    b_k = sp.array(sp.mat(MAT_inv) * sp.mat(rho).T).ravel()

    # Solve 'a' with ys

    Kop = LinearOperator((num_data, num_data), matvec=K_BB_E)
    a_star = sp.sparse.linalg.minres(Kop, ys, tol=1e-9)

    # Compute posterior mean

    post_means = K_opt(A_sparse.dot(a_star[0]))

    # If target_seqs is specified, pick up the corresponding components
    # if target_seqs is not None:
    #     poss = []
    #     for target_seq in target_seqs:
    #         poss.append(sequence_to_position(target_seq))
    #     post_means = np.take(post_means, poss)
    #
    # # Return
    return post_means

# function for calculating MAP externally


def compute_posterior_mean_ex(seqs, ys, sig2s):
    """compute MAP for external use"""

    global lda_star

    set_data_as_global_parameters(seqs, ys, sig2s)
    construct_A_sparse()
    construct_E_sparse()

    betas = 10 ** np.arange(-2, 6, .5)

    tr = np.array([seq2pos(seqs[i]) for i in range(len(seqs))])

    out = lambdaCV(ys, tr, sig2s, betas, 10)
    beta_cv = out[1]
    lda_cv = solve_lambda_single_beta(ys, tr, sig2s, beta_cv)
    lda_star = lda_cv

    fstar = compute_posterior_mean()

    return lda_cv, fstar

# posterior variance


def compute_posterior_variance(target_seqs=None):
    """compute posterior variances for a list of sequences
    Keyword arguments:
    target_seqs -- list of sequences to compute for"""

    # Set global parameters for later use
    global b_k, seqlist

    seqlist = [sequence_to_position(seq) for seq in seqs]

    # Construct rho & b_k
    W_kd = W_kd_mat()
    rho = W_kd.T.dot(lda_star)
    b_k = np.array(np.mat(MAT_inv) * np.mat(rho).T).ravel()
    Kop = LinearOperator((num_data, num_data), matvec=K_BB_E)

    K_ii = rho[0]

    # If target_seqs is not specified, set it to all sequences
    if target_seqs is None:
        bases = np.array(range(a))
        target_seqs = list(itertools.product(bases, repeat=l))

    K_Bi = np.zeros([len(seqlist), len(target_seqs)])
    for i in range(len(target_seqs)):
        vec = np.zeros(G)
        vec[sequence_to_position(target_seqs[i])] = 1
        K_Bi[:, i] = np.array(K_opt(vec))[seqlist]

    # Compute posterior variance for each target sequence
    post_vars = []
    print("computing posterior variance")
    for i in tqdm(range(len(target_seqs))):

        alph = cg(Kop, K_Bi[:, i])[0]
        post_vars.append(K_ii - np.sum(K_Bi[:, i] * alph))

    # Return
    return np.array(post_vars)


# posterior sampling

def prepare_pos_sampling():
    """construct all objects needed for posterior sampling using HMC"""
    global E0
    construct_A_sparse()
    E0 = dia_matrix((A_sparse.dot(
        [1 / sig2s[i] for i in range(len(sig2s))]), np.array([0])), shape=(G, G))


def grad_U(q):
    """returns the gradient of the potential energy at position q"""
    return(M_opt(q, inv(lda_star)) + E0.dot(q))


def S(q):
    """returns the potential energy for position q"""
    return (1 / 2) * q.dot(grad_U(q))


def leapfrog(q, p, m, potential, path_len, step_size):
    """perform leapfrog integration

    Keyword arguments:
    q -- starting position vector
    p -- starting momentum vector
    m -- scaling factor for kinetic energy
    potential -- function for calculating the potential energy
    path_len --- leapfrog integration length
    step_size -- size of one leapfrog step
    """

    q, p = np.copy(q), np.copy(p)

    _, dVdq = potential(q)

    p -= step_size * dVdq / 2  # half step
    for _ in tqdm(np.arange(np.round(path_len / step_size) - 1)):
        q += step_size / m * p  # whole step
        V, dVdq = potential(q)
        p -= step_size * dVdq  # whole step
    q += step_size / m * p  # whole step
    V, dVdq = potential(q)
    p -= step_size * dVdq / 2  # half step

    # momentum flip at end
    return q, -p, V, dVdq


class DualAveragingStepSize:
    """update stepsize for the leapfrog function during tuning steps"""

    def __init__(self, initial_step_size, target_accept=0.5, gamma=0.05, t0=10.0, kappa=0.75):
        # proposals are biased upwards to stay away from 0.
        self.mu = np.log(10 * initial_step_size)
        self.target_accept = target_accept
        self.gamma = gamma
        self.t = t0
        self.kappa = kappa
        self.error_sum = 0
        self.log_averaged_step = 0

    def update(self, p_accept):
        # Running tally of absolute error. Can be positive or negative. Want to
        # be 0.
        self.error_sum += self.target_accept - p_accept

        # This is the next proposed (log) step size. Note it is biased towards
        # mu.
        log_step = self.mu - self.error_sum / (np.sqrt(self.t) * self.gamma)

        # Forgetting rate. As `t` gets bigger, `eta` gets smaller.
        eta = self.t ** -self.kappa

        # Smoothed average step size
        self.log_averaged_step = eta * log_step + \
            (1 - eta) * self.log_averaged_step

        # This is a stateful update, so t keeps updating
        self.t += 1

        # Return both the noisy step size, and the smoothed step size
        return np.exp(log_step), np.exp(self.log_averaged_step)


def potential(q):
    """return both potential energy and its gradient at position q"""
    return [S(q), grad_U(q)]


def hamiltonian_monte_carlo(
        n_samples,
        m,
        initial_position,
        tune=500,
        initial_step_size=0.1,
        num_steps=100,
        max_energy_change=1000.0,
        intermediate_output=True):
    """perform HMC sampling

    Keyword arguments:
    n_samples: number of samples to draw from the posterior
    m: scaling factor for the kinetic energy
    initial_position: starting position
    tune: number of tuning steps
    initial_step_size: initial leapfrog step size
    num_steps: number of leapfrog steps per iteraction
    max_energy_change: set maximum difference between starting ann and 
    final total energy: for numerical stability
    intermediate_output: output intermediate samples and variances

    """

    global step_size, acceptance_ratio

    num_acceptance = 0

    initial_position = np.array(initial_position)

    initial_potential, initial_potential_grad = potential(initial_position)

    frac = 0.02

    batch_size = 100

    # collect all our samples in a list
    samples = [initial_position]

    # Keep a single object for momentum resampling
    momentum = st.norm(0, m)

    step_size = initial_step_size
    step_size_tuning = DualAveragingStepSize(step_size)
    # If initial_position is a 10d vector and n_samples is 100, we want 100 x 10 momentum draws
    # we can do this in one call to np.random.normal, and iterate over rows
    size = (n_samples + tune,) + initial_position.shape[:1]
    for idx, p0 in enumerate(momentum.rvs(size=size)):
        print("currently working on step " + str(idx + 1))

        # num_steps_r = np.random.randint(int((1 - frac)*num_steps), int((1 + frac)*num_steps) + 1)
        # step_size_r =  np.random.uniform((1 - frac)*step_size, (1 + frac)*step_size)

        numerical_check = False

        while numerical_check == False:

            num_steps_r = np.random.randint(
                int((1 - frac) * num_steps), int((1 + frac) * num_steps) + 1)
            step_size_r = np.random.uniform(
                (1 - frac) * step_size, (1 + frac) * step_size)

            q_new, p_new, final_potential, final_dVdq = leapfrog(
                samples[-1],
                p0,
                m,
                potential,
                step_size_r * num_steps_r,
                step_size_r)

            numerical_check = isinstance(final_potential, (float))

        print("initial_potential = ", initial_potential)
        print("final_potential = ", final_potential)

        # print("log_df(p0) = ", np.sum(momentum.logpdf(p0)))
        # print("log_df(p_new) = ", np.sum(momentum.logpdf(p_new)))

        start_energy = -np.sum(momentum.logpdf(p0)) + initial_potential
        new_energy = -np.sum(momentum.logpdf(p_new)) + final_potential

        # print("new_log_p = ", str(new_log_p))

        if new_energy == float("inf"):
            new_energy = 1e+10

        energy_change = new_energy - start_energy

        if energy_change < -10:
            energy_change = -10

        print("energy_change = ", str(energy_change))

        # Check Metropolis acceptance criterion
        p_accept = min(1, np.exp(-energy_change))
        print("p accept = ", str(np.exp(-energy_change)))
        if np.random.rand() < p_accept:
            samples.append(q_new)
            initial_potential = final_potential
            initial_potential_grad = final_dVdq
            num_acceptance += 1
        else:
            samples.append(np.copy(samples[-1]))

        # output batch samples&variance

        if intermediate_output is True:
            if idx + 1 > tune and (idx + 1 - tune) % batch_size == 0:
                varsample = np.var(samples[1 + tune:], axis=0)

                batch = int((idx + 1 - tune) / batch_size)

                pd.DataFrame(varsample.T, index=seqsAll).to_csv(outpath + "/" + str(sample_name) +
                                                                "_hmc_vars_batch_" + str(batch) + ".txt", index=seqsAll, header=None)

                pd.DataFrame(np.array(samples[1 + tune + (batch - 1) * batch_size: 1 + tune + batch * batch_size]).T, index=seqsAll).to_csv(
                    outpath + "/" + str(sample_name) + "_hmc_samples_batch_" + str(batch) + ".txt", index=seqsAll, header=None)

        if idx < tune - 1:
            step_size, _ = step_size_tuning.update(p_accept)
        elif idx == tune - 1:
            _, step_size = step_size_tuning.update(p_accept)

        acceptance_ratio = num_acceptance / (idx + 1)
        print("%s%s" % ("current step size = ", step_size))
        print("%s%s" % ("current acceptance ratio = ", acceptance_ratio))

    return np.array(samples[1 + tune:]), np.var(samples[1 + tune:], axis=0)


def compute_pot_scale_red_facs(components, multi_phi_samples0, from_to):
    """
    calculate the potential scale reduction factors for all sequences
    """

    # Crop the original multi_phi_samples
    multi_phi_samples = multi_phi_samples0[:, :, from_to[0]:from_to[1]]

    num_simulations, G, num_samples = \
        multi_phi_samples.shape[0], multi_phi_samples.shape[
            1], multi_phi_samples.shape[2]

    num_chains, len_chain = 2 * num_simulations, int(num_samples / 2)

    # Re-shape the multi phi samples into a shape of (num_chains, G, len_chain)
    a = []
    for k in range(num_simulations):
        a.append(multi_phi_samples[k, :, :len_chain])
        a.append(multi_phi_samples[k, :, len_chain:])

    multi_phi_samples_reshaped = np.array(a)

    # Compute potential scale reduction factor for the specified components of
    # phi
    if type(components) == int:
        I = range(components)
    elif type(components) == list:
        I = components
    else:
        print('"components" not recognized.')
        sys.exit()

    pot_scale_red_facs = []
    for i in I:

        # Collect the chains of samples of phi_i
        i_collector = np.zeros([len_chain, num_chains])
        for j in range(num_chains):
            i_collector[:, j] = multi_phi_samples_reshaped[j, i, :]

        # Compute the between-chain variance
        mean_0 = i_collector.mean(axis=0)
        mean_01 = mean_0.mean()
        B = len_chain / (num_chains - 1) * np.sum((mean_0 - mean_01)**2)

        # Compute the within-chain variance
        s2 = np.zeros(num_chains)
        for j in range(num_chains):
            s2[j] = 1 / (len_chain - 1) * \
                np.sum((i_collector[:, j] - mean_0[j])**2)
        W = s2.mean()

        # Estimate the marginal posterior variance
        var = (len_chain - 1) / len_chain * W + 1 / len_chain * B

        # Compute the potential scale reduction factor
        pot_scale_red_fac = sp.sqrt(var / W)

        # Save
        pot_scale_red_facs.append(pot_scale_red_fac)

    # Return
    return np.array(pot_scale_red_facs)


def plot_trajectory(ij, starting_phis, multi_phi_samples0, from_to, colors, xlimits, ylimits):
    """plot HMC trajectory for specify samples i and j"""

    # Crop the original multi_phi_samples
    multi_phi_samples = multi_phi_samples0[:, :, from_to[0]:from_to[1]]

    num_simulations, G, num_samples = \
        multi_phi_samples.shape[0], multi_phi_samples.shape[
            1], multi_phi_samples.shape[2]

    # Plot trajectory of i-j components of phi
    i, j = ij[0], ij[1]

    if i == j:

        plt.figure(0, figsize=(6, 5))

        for k in range(num_simulations):

            phi_start = starting_phis[:, k]
            phi_samples = multi_phi_samples[k, :, :]

            plt.figure(0)
            plt.scatter(sp.array([0]), phi_start[i],
                        marker='x', s=150, color=colors[k], zorder=2)
            plt.scatter(sp.array([0]), phi_samples[i, 0],
                        marker='+', s=200, color=colors[k], zorder=2)
            plt.plot(sp.array(range(num_samples)), phi_samples[
                     i, :], color=colors[k], a=0.4, zorder=1)

        plt.xlabel('Iteration', fontsize=14)
        plt.ylabel(r'$\phi_{%s}$' % np.array(seqsAll)[i], fontsize=16)
        if (xlimits is None) and (ylimits is None):
            plt.xlim(0,)
        else:
            plt.xlim(xlimits)
            plt.ylim(ylimits)

        plt.savefig(outpath + '/R_hat.png')
        plt.show()

    else:

        plt.figure(0, figsize=(6, 6))

        for k in range(num_simulations):

            phi_start = starting_phis[:, k]
            phi_samples = multi_phi_samples[k, :, :]

            plt.figure(0)
            plt.scatter(phi_start[i], phi_start[j],
                        marker='x', s=150, color=colors[k], zorder=2)
            plt.scatter(phi_samples[i, 0], phi_samples[
                        j, 0], marker='+', s=200, color=colors[k], zorder=2)
            plt.plot(phi_samples[i, :], phi_samples[j, :],
                     color=colors[k], a=0.4, zorder=1)

        plt.xlabel(r'$\phi_{%d}$' % i, fontsize=16)
        plt.ylabel(r'$\phi_{%d}$' % j, fontsize=16)
        plt.xlim(xlimits)
        plt.ylim(ylimits)

        plt.savefig(outpath + '/trajectories.png')
        plt.show()


def combine_samples(multi_phi_samples0, from_to):
    """combine multiple hmc samples for compute_pot_scale_red_facs()"""

    # Crop the original multi_phi_samples
    multi_phi_samples = multi_phi_samples0[:, :, from_to[0]:from_to[1]]

    num_simulations, G, num_samples = \
        multi_phi_samples.shape[0], multi_phi_samples.shape[
            1], multi_phi_samples.shape[2]

    # Combine the selected phi samples
    phi_samples = np.zeros([G, num_simulations * num_samples])
    for k in range(num_simulations):
        phi_samples[
            :, k * num_samples:(k + 1) * num_samples] = multi_phi_samples[k, :, :]

    # Convert phi samples to Q samples
    Q_samples = np.zeros([G, num_simulations * num_samples])
    for k in range(num_simulations * num_samples):
        Q_samples[:, k] = sp.exp(-phi_samples[:, k]) / \
            sp.sum(sp.exp(-phi_samples[:, k]))

    # Return
    return phi_samples, Q_samples


def plot_distribution(components, map_estimate, samples, num_bins, colors, xlimits):

    """Plot posterior distribution of specified components of phi or Q"""
    if type(components) == int:
        I = range(components)
    elif type(components) == list:
        I = components
    else:
        print('"components" not recognized.')
        sys.exit()

    for i in I:

        print('Component # %d' % i)
        if map_estimate is not None:
            print('MAP estimate = %.3f' % map_estimate[i])

        plt.figure(i, figsize=(6, 4))

        for j in range(len(samples)):

            ci = np.percentile(samples[j][i, :], [2.5, 97.5])
            print('95%% CI = [%.3f, %.3f]' % (ci[0], ci[1]))

            hist, bin_edges = np.histogram(
                samples[j][i, :], bins=num_bins[j], density=True)
            bin_width = bin_edges[1] - bin_edges[0]
            bin_centers = np.linspace(
                bin_edges[0] + bin_width / 2, bin_edges[-1] - bin_width / 2, len(bin_edges) - 1)

            plt.figure(i)
            plt.bar(bin_centers, hist, width=bin_width,
                    color=colors[j], a=0.5, edgecolor=colors[j])

        plt.xlim(xlimits)
        plt.show()
