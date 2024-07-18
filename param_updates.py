import numpy as np
from scipy.stats import norm
from scipy.stats import gamma

def logit(x):
    return np.log(x / (1 - x))

def expit(x):
    return 1 / (1 + np.exp(-x))

def sample_beta(tau2, theta, Z, sig2, num_island_region, isl_reg, num_group, num_island, A, m0):
    var_beta = (np.tile(tau2, num_island) / np.repeat(num_island_region, num_group)).reshape(num_island, num_group)
    tmZ = theta - Z
    mean_beta = [(tmZ[isl_reg[isl]]).mean(0) for isl in range(num_island)]
    var_t = tau2 + (tau2 + sig2) / m0
    beta_thres = ((1 - A) + np.sqrt((A - 1) ** 2 + 4 * (A - 1 / var_t))) / 2
    beta_thres = np.array([logit(np.maximum(0, beta_thres))] * num_island)
    beta_max = norm.cdf(beta_thres, mean_beta, np.sqrt(var_beta))
    if (beta_max == 0).any():
        return beta
    u = np.random.uniform(0, beta_max)
    beta = norm.ppf(u, mean_beta, np.sqrt(var_beta))
    return beta

def sample_Z(Z, tau2, sig2, theta, beta, num_adj, island_id, adj, num_region):
    tmb = theta - beta[island_id]
    var_Z = [1 / (1 / tau2 + i / sig2) for i in range(max(num_adj) + 1)]
    for i in range(num_region):
        mean_Z = var_Z[num_adj[i]] * (tmb[i] / tau2 + Z[adj[i]].sum(0) / sig2)
        Z[i] = np.random.normal(mean_Z, np.sqrt(var_Z[num_adj[i]]))
    Z -= Z.mean(0)
    return Z

def sample_sig2(beta, Z, tau2, num_island_region, adj, num_adj, num_region, num_group, num_island, sigma_a, sigma_b, m0, A):
    a_sig = (num_region - num_island) / 2 + sigma_a
    sum_adj = np.array([Z[i] * Z[adj[i]].sum(0) for i in range(num_region)]).sum(0)
    b_sig = 1 / (((num_adj * Z.T ** 2).sum(1) - sum_adj) / 2 + sigma_b)
    pi = (beta * np.repeat(num_island_region, num_group).reshape(num_island, num_group) / num_region).sum(0)
    pi = np.exp(pi) / (1 + np.exp(pi))
    sig_thres = (1 / ((A + pi) * (1 - pi)) - tau2 * (1 + 1 / m0)) * m0
    sig_thres = np.maximum(0, sig_thres)
    with np.errstate(divide = "ignore"):
        sig_max = gamma.cdf(1 / sig_thres, a_sig, scale = b_sig)
    u = np.random.uniform(0, sig_max)
    sig2 = 1 / gamma.ppf(u, a_sig, scale = b_sig)
    return sig2

def sample_tau2(theta, beta, Z, sig2, island_id, num_island_region, num_region, num_group, num_island, tau_a, tau_b, A, m0):
    a_tau = num_region / 2 + tau_a
    b_tau = 1 / (((theta - beta[island_id] - Z) ** 2).sum(0) / 2 + tau_b)
    pi = (beta * np.repeat(num_island_region, num_group).reshape(num_island, num_group) / num_region).sum(0)
    pi = np.exp(pi) / (1 + np.exp(pi))
    tau_thres = (1 / ((A + pi) * (1 - pi)) - sig2 / m0) / (1 + 1 / m0)
    tau_thres = np.maximum(0, tau_thres)
    with np.errstate(divide = "ignore"):
        tau_max = gamma.cdf(1 / tau_thres, a_tau, scale = b_tau)
    u = np.random.uniform(0, tau_max, num_group)
    tau2 = 1 / gamma.ppf(u, a_tau, scale = b_tau)
    return tau2

def sample_theta(theta, beta, Z, tau2, Y, n, theta_sd, theta_acpt, island_id, num_region, num_group):
    theta_star = np.random.normal(theta, theta_sd)
    rk1 = Y * (theta_star - theta)
    rk2 = n * (np.log(1 + np.exp(theta_star)) - np.log(1 + np.exp(theta)))
    rk3a = (theta_star - beta[island_id] - Z) ** 2
    rk3b = (theta      - beta[island_id] - Z) ** 2
    rk = np.exp(rk1 - rk2 - 1 / (2 * tau2) * (rk3a - rk3b))
    rcand = rk >= np.random.uniform(0, 1, num_region * num_group).reshape(num_region, num_group)
    if rcand.any():
        theta[rcand] = theta_star[rcand]
        theta_acpt[rcand] += 1e-2
    return [theta, theta_acpt]
