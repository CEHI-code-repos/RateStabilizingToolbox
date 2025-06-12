try:
    import arcpy
except ImportError:
    arcpy_available = False
else:
    arcpy_available = True
import numpy as np
import pandas as pd
import json
from . import param_updates

# Calculate sets of contiguous adjacency structures
def get_islands(adj):
    f = set(range(len(adj)))
    isl_reg = []
    while f:
        active_list = {next(iter(f))}
        inactive_list = set()
        while active_list:
            Na = adj[active_list.pop()]
            active_list |= set(Na) - inactive_list
            inactive_list |= active_list
        isl_reg.append(sorted(inactive_list))
        f -= inactive_list
    return isl_reg

# Create initial values
def get_inits(Y, n, isl_reg, island_id):
    with np.errstate(divide = "ignore"):
        beta = np.array([param_updates.logit(Y[i].sum(0) / n[i].sum(0)) for i in isl_reg])
        theta = param_updates.logit(Y / n)
    beta[~np.isfinite(beta)] = param_updates.logit(Y.sum() / n.sum())
    theta[~np.isfinite(theta)] = beta[island_id][~np.isfinite(theta)]
    tau2, sig2 = [theta.var(0) * 2] * 2
    Z = theta - beta[island_id]
    return [theta, beta, Z, tau2, sig2]

# Create spatial data
def get_spdat(adj):
    if (min(min(adj)) == 1):
        adj = [list(np.array(adj[i]) - 1) for i in range(len(adj))]
    num_adj = [len(adj[i]) for i in range(len(adj))]
    isl_reg = get_islands(adj)
    num_island = len(isl_reg)
    num_island_region = [len(isl_reg[isl]) for isl in range(num_island)]
    island_id = [0] * len(adj)
    for isl in list(range(0, num_island)):
        island_id = [isl if i in isl_reg[isl] else x for i, x in enumerate(island_id)]
    return [adj, num_adj, isl_reg, num_island, num_island_region, island_id]

# Restricted UCAR Gibbs sampler
def gibbs_rucar(Y, n, adj, std_pop):
    np.random.seed(1) # For replicability
    adj, num_adj, isl_reg, num_island, num_island_region, island_id = get_spdat(adj) # Spatial data
    theta, beta, Z, tau2, sig2 = get_inits(Y, n, isl_reg, island_id) # Inits
    tau_a, tau_b, sigma_a, sigma_b = [1e-3] * 4 # Hyperparameters
    theta_sd = theta * 0 + 0.025
    theta_acpt = theta * 0
    num_region, num_group = Y.shape # Auxiliary variables
    m0 = 3
    A = Y.sum(0) / n.sum(0) * std_pop
    A = 6 * A / sum(A)
    theta_out = np.zeros([num_region, num_group, 400])
    if arcpy_available: arcpy.SetProgressor("step", "Generating estimates...", 0, 6000, 1)
    for s in range(6000):
        sig2 = param_updates.sample_sig2(beta, Z, tau2, num_island_region, adj, num_adj, num_region, num_group, num_island, sigma_a, sigma_b, m0, A)
        tau2 = param_updates.sample_tau2(tau2, theta, beta, Z, sig2, island_id, num_island_region, num_region, num_group, num_island, tau_a, tau_b, A, m0)
        beta = param_updates.sample_beta(beta, tau2, theta, Z, sig2, num_island_region, isl_reg, num_group, num_island, A, m0)
        Z = param_updates.sample_Z(Z, tau2, sig2, theta, beta, num_adj, island_id, isl_reg, adj, num_region, num_island)
        theta, theta_acpt = param_updates.sample_theta(theta, beta, Z, tau2, Y, n, theta_sd, theta_acpt, island_id, num_region, num_group)
        if (s + 1) % 100 == 0:
            theta_acpt = np.clip(theta_acpt, 0.20, 0.75)
            theta_sd  *= theta_acpt / 0.43
            theta_acpt *= 0
        if (s > 2000) & ((s + 1) % 10 == 0):
            theta_out[:, :, (s - 2000) // 10] = theta
        if arcpy_available: arcpy.SetProgressorPosition()
        print(str(round(s / 6000 * 100, 1)) + "%", end = "\r")
    print()
    return theta_out

# Age-standardization algorithm
def age_std(output, ages, std_pop, wtages):
    wt = std_pop[np.isin(ages, wtages)]
    wt = (wt / sum(wt)).reshape(1, len(wtages), 1)
    output_new = (output[:, np.where(np.isin(ages, wtages))[0], :] * wt).sum(1, keepdims = True)
    return np.concatenate((output, output_new), 1)

# ***MODIFIED*** Now outputs minimum CI for reliable estimates in each region, along with true/false value for reliability. Does not suppress estimates.
# Calculate medians
def get_medians(output, regions, ages, ci_pct):
    ci_values = [0.5, 0.75, 0.9, 0.95, 0.99]
    medians = pd.DataFrame(np.median(output, 2), regions, ages)
    alpha = (1 - ci_pct) / 2
    ci_lo = pd.DataFrame(np.quantile(output, alpha, 2), regions, ages)
    ci_hi = pd.DataFrame(np.quantile(output, 1 - alpha, 2), regions, ages)
    ci_chart = np.zeros(medians.shape)
    for ci in ci_values:
        alpha = (1 - ci) / 2
        lo = np.quantile(output, alpha, 2)
        hi = np.quantile(output, 1 - alpha, 2)
        rp = medians / (hi - lo)
        ci_chart[rp >= 1] = ci
    ci_chart = pd.DataFrame(ci_chart, regions, ages)
    if ages == [""]:
        medians = medians.rename(columns={"": "median"})
        ci_lo = ci_lo.rename(columns={"": "ci_low"})
        ci_hi = ci_hi.rename(columns={"": "ci_high"})
        ci_chart = ci_chart.rename(columns={"": "max_reliable_ci"})
    return [medians, ci_lo, ci_hi, ci_chart]

# Model testing
def get_estimates(data_url, adj_url, region_id, event_id, pop_id, estimates_dir, estimates_name, rates_per = 1e5, group_id = None, std_pop_yr = None, age_std_groups = None, age_std_group_names = None):
    print("Preparing data...")
    # Import data
    data = pd.read_csv(data_url)
    with open(adj_url) as f:
        adj = json.load(f)

    ### ***ADDED*** contingencies for single-group data
    age_groups = [""]
    num_group = 1
    data = data.sort_values(by = [region_id])
    regions = data[region_id].unique().tolist()
    num_region = data[region_id].nunique()
    if group_id:
        data = data.sort_values(by = [region_id, group_id])
        age_groups = data[group_id].unique().tolist()
        num_group = data[group_id].nunique()
    

    Y = np.array(data[event_id]).reshape([num_region, num_group])
    n = np.array(data[pop_id]).reshape([num_region, num_group])

    # Set up standard population
    # ***ADDED*** Contingencies for data that does not age-standardize
    std_pop = 1
    # Klein, et al. 2000 standard population
    if (std_pop_yr == "2000"): 
        std_pop = np.array([18987, 39977, 38077, 37233, 44659, 37030, 23961, 18136, 12315, 4259])
    # 2010 standard population
    if (std_pop_yr == "2010"):
        std_pop = np.array([20203362, 41025851, 43626342, 41063948, 41070606, 45006716, 36482729, 21713429, 13061122, 5493433])
    if std_pop_yr:
        pop_ages = ["0-4", "5-14", "15-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85up"]
        std_pop = std_pop[np.isin(pop_ages, age_groups)]

    # Generate estimates
    print("Generating estimates...")
    theta_out = gibbs_rucar(Y, n, adj, std_pop)
    output = param_updates.expit(theta_out) * rates_per
    np.save(estimates_dir + "/" + "output.npy", output)
    if age_std_groups:
        for ages in age_std_groups:
            output = age_std(output, age_groups, std_pop, ages)
        age_groups.extend(age_std_group_names)
    medians, ci_chart, reliable = get_medians(output, regions, age_groups)

    # ***ADDED*** Saves CI chart and reliability logical into same folder as medians. Also, the region_id column is now named for all three csv's
    medians.to_csv(estimates_dir + "/" + estimates_name, index_label = region_id)
    ci_chart.to_csv(estimates_dir + "/" + "ci_chart.csv", index_label = region_id)
    reliable.to_csv(estimates_dir + "/" + "reliable.csv", index_label = region_id)
    print("Model finished!")
