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