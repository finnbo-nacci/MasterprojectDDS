"""
Optimal Transport Transfer Operator approximation (OTTO)
approximate the singular or eigenvalues of the transfer operator from a entropic optimal transport plan
inspired on https://arxiv.org/abs/2006.16085
"""

import numpy as np
import matplotlib.pyplot as plt

def l22(mu, nu):
    """computes the squared Euclidean distance l_2^2"""
    return np.sum(mu ** 2, 0)[:, None] + np.sum(nu ** 2, 0)[None, :] - 2 * mu.transpose().dot(nu)


def l22_1d(mu, nu):
    return (mu ** 2)[:, None] + (nu ** 2)[None, :] - 2 * mu.transpose().dot(nu)


def Sinkhorn(a, b, C, epsilon, threshold=1e-12, max_iter=5000):
    # a,b: 2d n x m array, for m datapoints of n dimensions
    # C  : m x m cost matrix
    # eps: regularisation parameter

    K = np.exp(-C / epsilon)  # Gibbs kernel K=e^{-C/epsilon}
    v = np.ones(a.size)  # initialize
    u = a

    i = 0
    residuals_a = np.empty(max_iter + 1)  # make arrays to store residuals
    residuals_b = np.empty(max_iter + 1)
    residuals_a[i] = np.sum(np.abs(np.multiply(u, np.dot(K, v)) - a))  # P1 = u . Kv, so ||P1-a|| = ||u.Kv-a||
    residuals_b[i] = np.sum(np.abs(np.multiply(v, np.dot(np.transpose(K), u)) - b))  # v . K^T u

    # Repeat the iteration till convergence or max. iterations
    while i < max_iter and (residuals_a[i] > threshold or residuals_b[i] > threshold):
        i += 1
        u = a / (np.dot(K, v))
        residuals_b[i] = np.sum(np.abs(np.multiply(v, np.dot(np.transpose(K), u)) - b))
        v = b / (np.dot(np.transpose(K), u))
        residuals_a[i] = np.sum(np.abs(np.multiply(u, np.dot(K, v)) - a))
    P = np.dot(np.dot(np.diag(u), K), np.diag(v)) # compute the resulting optimal transport plan
    print("Sinkhorn iterations:", i, "out of", max_iter)
    return P

def otto(mu, nu, eps=.1):
    """computes the OTTO
    inputs:
    - mu: initial positions, N,d array
    - nu: final positions, N,d array
    - eps: regularisation strength epsilon, float
    outputs:
    - model: N*N array
    - evals : dominant eigenvalues, k array
    - evecs : dominant eigenvectors
    """

    # Initialize
    n = mu.shape[1]
    if n != nu.shape[1]: print("Mu and nu differ in number of data points. \n mu: %d, nu: %d" % (n, nu.shape[1]))
    a = np.full(n, 1/n) #np.ones(n) / n

    # Compute distance/cost matrix
    cost = l22(mu, nu)

    # Solve EOT problem
    eot_plan = Sinkhorn(a, a, cost, eps, threshold=1e-8, max_iter=5000)
    model = eot_plan * n  # *n == /a
    return model


def traj_to_distr(data, tau):
    # data is np array containing coordinates of trajectory
    samples_t0 = np.arange(0, data.shape[0] - tau, tau, dtype=int)
    samples_t1 = samples_t0 + tau

    mu = np.transpose(data[samples_t0])
    nu = np.transpose(data[samples_t1])
    return mu, nu



def estimate_model(traj, tau=100, eps=.1):
    mu, nu = traj_to_distr(traj, tau)
    return otto(mu, nu, eps)



def display_plan_traj(P, traj, ts):
    f, (ax1, ax2) = plt.subplots(2, figsize=(5, 7.5),
                                 gridspec_kw={'height_ratios': [1, 2]})
    ax1.plot(ts, traj, ts[::1000], traj[::1000], '.')
    ax1.set_title('Trajectory')
    ax1.margins(x=0, y=.2)
    ax1.set(xticklabels=[])
    ax2.imshow(np.log(P + 1e-5))
    ax2.axis('off')
    ax2.set_title('Entropic optimal transport plan')
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    f.subplots_adjust(hspace=0)
    #plt.tight_layout()
    plt.show()
