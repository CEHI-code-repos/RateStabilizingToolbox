import requests
import pandas as pd
from . import constants

def get_census(survey: str, year: str, geography: str, state: str, age_stratified: bool):
    geography = geography.lower()
    fips = constants.state_to_fips[state]

    if survey == "5-year ACS":
        survey = "acs"
        file = "acs5"
        year = year.split("-")[1]
    elif survey == "Decennial":
        survey = "dec"
        if int(year) == 2020:
            file = "dp"
        elif int(year) != 2020:
            file = "sf1"

    age_vars_dict = getattr(constants, file + "_age_vars")
    if age_stratified:
        var_cols = [var for var_list in age_vars_dict.values() for var in var_list]
    else:
        var_cols = [getattr(constants, file + "_tot_var")]
    var_str = ",".join(var_cols)

    req_url = f"https://api.census.gov/data/{year}/{survey}/{file}?get={var_str},GEO_ID,NAME&for={geography}:*&in=state:{fips}"

    resp = requests.get(req_url)
    resp_df = pd.DataFrame.from_dict(resp.json())
    resp_df.columns = resp_df.iloc[0].to_list()
    resp_df = resp_df.drop(0)

    geo_cols = [col for col in resp_df.columns if col not in var_cols]
    resp_df = resp_df[geo_cols + var_cols]
    resp_df["GEO_ID"] = resp_df["GEO_ID"].map(lambda geoid: geoid.split("US")[1])
    resp_df[var_cols] = resp_df[var_cols].apply(pd.to_numeric)
    if age_stratified:
        for age_grp in age_vars_dict:
            resp_df[age_grp] = resp_df[age_vars_dict[age_grp]].sum(axis=1)
        resp_df = resp_df.drop(var_cols, axis = 1)
        resp_df = resp_df.melt(
            id_vars = geo_cols, 
            value_vars = age_vars_dict.keys(),
            var_name = 'age_group', 
            value_name = 'pop_count')
    else:
        resp_df = resp_df.rename(columns = {getattr(constants, file + "_tot_var"): "pop_count"})

    resp_df = resp_df.rename(columns = {"GEO_ID": "geoid", "NAME": "name"})

    return resp_df