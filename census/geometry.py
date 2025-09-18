from . import constants
import arcpy
import tempfile
import requests
import os
from zipfile import ZipFile

def get_geometry(geometry: str, geom_type: str, year: str, state: str, output_path: str):
    state_geoid = constants.STATE_TO_FIPS[state]
    if year in constants.ACS_YEARS:
        year = year.split("-")[1]

    with tempfile.TemporaryDirectory() as tempdir:
        if geometry == "County":
            req_url, where_exp = get_county_url(geom_type, year, state_geoid)
        elif geometry == "Tract":
            req_url, where_exp = get_tract_url(geom_type, year, state_geoid)

        resp = requests.get(req_url)
        temp_zip_path = os.path.join(tempdir, "temp.zip")

        with open(temp_zip_path, "wb") as f:
            f.write(resp.content)
        
        with ZipFile(temp_zip_path, 'r') as zip_ref:
            zfiles = zip_ref.namelist()
            zip_ref.extractall(tempdir)

        shape_file = [file for file in zfiles if file.endswith('.shp')][0]
        shape_file_path = os.path.join(tempdir, shape_file)

        arcpy.conversion.ExportFeatures(shape_file_path, output_path, where_exp)

def get_county_url(geom_type: str, year: str, state_geoid: str) -> tuple[str, str]:
    year_int = int(year)
    year_suffix = year[2:]

    if geom_type == "Cartographic":
        if year_int >= 2014:
            req_url = f"https://www2.census.gov/geo/tiger/GENZ{year}/shp/cb_{year}_us_county_500k.zip"
            where_exp = f"STATEFP = '{state_geoid}'"
        elif year_int == 2013:
            req_url = f"https://www2.census.gov/geo/tiger/GENZ{year}/cb_{year}_us_county_500k.zip"
            where_exp = f"STATEFP = '{state_geoid}'"
        elif year_int == 2010:
            req_url = f"https://www2.census.gov/geo/tiger/GENZ2010/gz_2010_us_050_00_500k.zip"
            where_exp = f"STATE = '{state_geoid}'"
        elif year_int == 2000:
            req_url = f"https://www2.census.gov/geo/tiger/PREVGENZ/co/co{year_suffix}shp/co99_d{year_suffix}_shp.zip"
            where_exp = f"ST = '{state_geoid}'"
    else:
        if year_int >= 2011:
            req_url = f"https://www2.census.gov/geo/tiger/TIGER{year}/COUNTY/tl_{year}_us_county.zip"
            where_exp = f"STATEFP = '{state_geoid}'"
        elif year_int in [2000, 2010]:
            req_url = f"https://www2.census.gov/geo/tiger/TIGER2010/COUNTY/{year}/tl_2010_us_county{year_suffix}.zip"
            where_exp = f"STATEFP{year_suffix} = '{state_geoid}'"
    
    return (req_url, where_exp)

def get_tract_url(geom_type: str, year: str, state_geoid: str) -> tuple[str, str]:
    year_int = int(year)
    year_suffix = year[2:]

    if geom_type == "Cartographic":
        if year_int >= 2014:
            req_url = f"https://www2.census.gov/geo/tiger/GENZ{year}/shp/cb_{year}_{state_geoid}_tract_500k.zip"
        elif year_int == 2013:
            req_url = f"https://www2.census.gov/geo/tiger/GENZ{year}/cb_{year}_{state_geoid}_tract_500k.zip"
        elif year_int == 2010:
            req_url = f"https://www2.census.gov/geo/tiger/GENZ2010/gz_2010_{state_geoid}_140_00_500k.zip"
        elif year_int == 2000:
            req_url = f"https://www2.census.gov/geo/tiger/PREVGENZ/tr/tr{year_suffix}shp/tr{state_geoid}_d{year_suffix}_shp.zip"
    else:
        if year_int >= 2011:
            req_url = f"https://www2.census.gov/geo/tiger/TIGER{year}/TRACT/tl_{year}_{state_geoid}_tract.zip"
        elif year_int in [2000, 2010]:
            req_url = f"https://www2.census.gov/geo/tiger/TIGER2010/TRACT/{year}/tl_2010_{state_geoid}_tract{year_suffix}.zip"

    return (req_url, None)
    