import geopandas as gpd
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict, List, Union
import model

AGE_GROUPS = ["0-4", "5-14", "15-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85up"]

def rst(
    input_table_path: str, 
    input_table_fields: Dict[str, str], 
    input_feature_path: str, 
    input_feature_fields: Dict[str, str], 
    additional_options: Dict[str, Union[float, int]], 
    age_group_field: Optional[str] = None,
    std_pop_year: Optional[int] = None, 
    std_age_groups: Optional[Dict[str, List[str]]] = None
) -> pd.DataFrame:
    """Get RST estimates

    Parameters
    ----------
    input_table_path : str
        The file location of the table containing event and population data
        at the region or region-age group level. Geodatabase or geopackage
        features can be specified as "db.gdb/table" or "pkg.gpkg/table"
    input_table_fields : Dict[str, str]
        Dictionary containing "region_id", "event_count", and 
        "population_count" keys which specify their respective field names 
        within the input table.
    input_feature_path : str
        The file location of the table containing event and population data
        at the region or region-age group level. Geodatabase or geopackage
        features can be specified as "db.gdb/table" or "pkg.gpkg/table"
    input_table_fields : Dict[str, str]
        Dictionary containing "region_id" key which specifies the region id
        field name within the input feature.
    additional_options : Dict[str, Union[float, int]]
        Dictionary containing "ci" key which specifies level of confidence
        used to generate confidence intervals (e.g. 0.95), "rate_per" key 
        which specifies the denomerator of output rates (e.g. 100,000), and
        "n_years" key which specifies the total number of years of data
        represented by the input table events.
    age_group_field : str, optional
        Name of age group field within the input table.
    std_pop_year : int, optional
        Standard population year used to age adjust output rates. 2000 or 
        2010.
    std_age_groups : Dict[str, List[str]], optional
        Dictionary with keys representing the names of groups to age
        standardize for (e.g. "65up") and values being a list of 
        constituent age groups (e.g. ["65-74", "75-84", "85up"])

    Returns
    -------
    pd.DataFrame
        A dataframe containing the output rates.

    Notes
    _______
    This is a RST implementation without the arcgis dependency, should allow for easier testing
    there will not be as many checks on the data, so maybe run it in ArcGIS Pro first.

    Examples
    --------
    >>> rst(
    >>>     input_table_path = "data.gdb/mi_joined_event_pop_strat",
    >>>     input_table_fields = {
    >>>         "region_id": "GEOID",
    >>>         "event_count": "EventCount",
    >>>         "population_count": "PopulationCount"
    >>>     },
    >>>     input_feature_path = "data.gdb/mi_carto",
    >>>     input_feature_fields = {
    >>>         "region_id": "GEOID"
    >>>     },
    >>>     additional_options = {
    >>>         "ci": 0.95,
    >>>         "rate_per": 100_000,
    >>>         "n_years": 1
    >>>     },
    >>>     age_group_field = "AgeGroup",
    >>>     std_pop_year = 2000,
    >>>     std_age_groups = {
    >>>         "35-64": ["35-44", "45-54", "55-64"],
    >>>         "35up": ["35-44", "45-54", "55-64", "65-74", "75-84", "85up"],
    >>>         "65up": ["65-74", "75-84", "85up"]
    >>>     }
    >>> )
    """
    input_table_id_name = input_table_fields["region_id"]
    input_table_event_name = input_table_fields["event_count"]
    input_table_pop_name = input_table_fields["population_count"]
    input_table = read_geom(input_table_path).sort_values(by = [input_table_id_name])

    input_feature_id_name = input_feature_fields["region_id"]
    input_feature = read_geom(input_feature_path).sort_values(by = [input_feature_id_name])
    input_feature_adj = calculate_adj(input_feature)

    region_ids = input_table[input_table_id_name].unique().tolist()
    num_region = input_table[input_table_id_name].nunique()

    ci = additional_options["ci"]
    rate_per = additional_options["rate_per"]
    n_years = additional_options["n_years"]

    n_age_groups = 1
    if std_pop_year:
        if std_pop_year == 2000: 
            std_pop = np.array([18987, 39977, 38077, 37233, 44659, 37030, 23961, 18136, 12315, 4259])
        elif std_pop_year == 2010:
            std_pop = np.array([20203362, 41025851, 43626342, 41063948, 41070606, 45006716, 36482729, 21713429, 13061122, 5493433])
        else:
            raise ValueError(f"No standard population found for year: {std_pop_year}")

        input_table_age_name = age_group_field
        input_table = input_table.sort_values(by = [input_table_id_name, input_table_age_name])
        age_groups = input_table[input_table_age_name].unique().tolist()
        n_age_groups = input_table[input_table_age_name].nunique()
        std_pop = std_pop[np.isin(AGE_GROUPS, age_groups)]

    Y = np.array(input_table[input_table_event_name]).reshape([num_region, n_age_groups])
    n = np.array(input_table[input_table_pop_name]).reshape([num_region, n_age_groups])

    theta_out = model.runner.gibbs_rucar(Y, n, input_feature_adj, std_pop)
    output = model.param_updates.expit(theta_out) * rate_per / n_years
    
    if std_pop_year:
        for _, constituient_age_group in std_age_groups.items():
            output = model.runner.age_std(output, age_groups, std_pop, constituient_age_group)
        age_groups.extend(std_age_groups.keys())

    medians, ci_lo, ci_hi, ci_chart = model.runner.get_medians(output, region_ids, age_groups, ci)

    return(medians)

def read_geom(path_str: str) -> gpd.GeoDataFrame:
    path = Path(path_str)

    if path.suffix:
        return gpd.read_file(path_str)

    parent = path.parent
    container = parent.suffix.lower()

    if container not in [".gdb", ".gpkg"]:
        raise ValueError(f"Could not read path: {path_str}")

    return gpd.read_file(str(parent), layer=path.name)


def calculate_adj(feature: gpd.GeoDataFrame) -> list[list[int]]:
    adjacency: list[list[int]] = [[] for _ in range(len(feature))]
    sindex = feature.sindex

    for i, geom in enumerate(feature.geometry):
        possible_idx = list(sindex.intersection(geom.bounds))

        neighbors = []
        for j in possible_idx:
            if j == i:
                continue
            other = feature.geometry.iloc[j]
            if geom.touches(other):
                neighbors.append(j)

        adjacency[i] = neighbors

    return adjacency
