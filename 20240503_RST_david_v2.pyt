# -*- coding: utf-8 -*-
import arcpy
import pandas as pd
import numpy as np
import helpers

class Toolbox:
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Rate Stabilizing Toolbox"
        self.alias = "RST"

        # List of tool classes associated with this toolbox
        self.tools = [RST]


class RST:
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Rate Stabilizing Tool"
        self.description = "Accurate and precise estimates of local-level epidemiologic measures are critical to informing policy and program decisions, but they often require advanced statistical knowledge, programming/coding skills, and extensive computing power. In response, we developed the Rate Stabilizing Tool (RST), an ArcGIS-based tool that enables users to input their own record-level data to generate more reliable age-standardized measures of chronic disease (eg, prevalence rates, mortality rates) or other population health outcomes at the county or census tract levels. The RST uses 2 forms of empirical Bayesian modeling (nonspatial and spatial) to estimate age-standardized rates and 95% credible intervals for user-specified geographic units. The RST also provides indicators of the reliability of point estimates. In addition to reviewing the RST’s statistical techniques, we present results from a simulation study that illustrates the key benefit of smoothing. We demonstrate the dramatic reduction in root mean-squared error (rMSE), indicating a better compromise between accuracy and stability for both smoothing approaches relative to the unsmoothed estimates. Finally, we provide an example of the RST’s use. This example uses heart disease mortality data for North Carolina census tracts to map the RST output, including reliability of estimates, and demonstrates a subsequent statistical test."

    def getParameterInfo(self):
        """Define the tool parameters."""

        param_data_table = arcpy.Parameter(
            displayName="Input Table",
            name="InputTable",
            datatype="GPTableView",
            parameterType="Required",
            direction="Input"
        )

        param_data_fields = arcpy.Parameter(
            displayName="Input Table Fields",
            name="InputTableFields",
            datatype="GPValueTable",
            parameterType="Required",
            direction="Input"
        )
        param_data_fields.parameterDependencies = [param_data_table.name]
        param_data_fields.columns = [['Field', 'Region ID'], ['Field', 'Age Group'], ['Field', 'Event'], ['Field', 'Population']]
        param_data_fields.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'

        param_feature = arcpy.Parameter(
            displayName="Input Feature",
            name="InputFeature",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input"
        )
        param_feature.filter.list = ["Polygon"]

        param_feature_field = arcpy.Parameter(
            displayName="Input Feature Field",
            name="InputFeatureField",
            datatype="GPValueTable",
            parameterType="Required",
            direction="Input"
        )
        param_feature_field.parameterDependencies = [param_feature.name]
        param_feature_field.columns = [['Field', 'Region ID']]
        param_feature_field.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'

        param_std_pop_yr = arcpy.Parameter(
            displayName="Standard Population Year",
            name="StandardPopulationYear",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )
        param_std_pop_yr.filter.list = ["2000", "2010"]

        param_out_table = arcpy.Parameter(
            displayName="Output Table",
            name="OutputTable",
            datatype="GPTableView",
            parameterType="Required",
            direction="Output"
        )

        params = [
            param_data_table,
            param_data_fields,
            param_feature,
            param_feature_field,
            param_std_pop_yr,
            param_out_table
        ]
        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""
        data_url = parameters[0].valueAsText
        data_exists = False
        data_fields = parameters[1].values
        data_region_id = data_event_id = data_pop_id = data_region_id_type = ""
        if data_fields is not None:
            data_region_id = data_fields[0][0].value
            data_event_id = data_fields[0][2].value
            data_pop_id = data_fields[0][3].value

        feature_url = parameters[2].valueAsText
        feature_exists = False
        feature_fields = parameters[3].values
        feature_region_id = feature_region_id_type = ""
        if feature_fields is not None:
            feature_region_id = feature_fields[0][0].value

        estimates_out = parameters[5].valueAsText
        estimates_out_exists = False

        if data_url is not None and arcpy.Exists(data_url):
            data_exists = True
        if feature_url is not None and arcpy.Exists(feature_url):
            feature_exists = True
        if estimates_out is not None and arcpy.Exists(estimates_out):
            estimates_out_exists = True

        # Check if Region ID types are the same
        if data_exists and feature_exists:
            data_region_id_type = [f.type for f in arcpy.ListFields(data_url) if f.name == data_region_id]
            feature_region_id_type = [f.type for f in arcpy.ListFields(feature_url) if f.name == feature_region_id]
            if len(data_region_id_type) + len(feature_region_id_type) != 2:
                pass
            elif data_region_id_type[0] != feature_region_id_type[0]:
                parameters[3].setErrorMessage("Input Feature Region ID field type does not match Input Table Region ID field type")

        # Check if Population ID is an integer
        if data_exists:
            data_pop_id_type = [f.type for f in arcpy.ListFields(data_url) if f.name == data_pop_id]
            if len(data_pop_id_type) != 1:
                pass
            elif data_pop_id_type[0] not in ["SmallInteger", "Integer", "BigInteger"]:
                parameters[1].setErrorMessage("Population field type is not a integer")

        # Check if Event ID is an integer
        if data_exists:
            data_event_id_type = [f.type for f in arcpy.ListFields(data_url) if f.name == data_event_id]
            if len(data_event_id_type) != 1:
                pass
            elif data_event_id_type[0] not in ["SmallInteger", "Integer", "BigInteger"]:
                parameters[1].setErrorMessage("Event field type is not a integer")

        # Check if Output Table exists
        if estimates_out_exists:
            parameters[5].setErrorMessage("Output Table already Exists")

        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        # Process parameters
        data_url = parameters[0].valueAsText
        data_fields = parameters[1].values
        data_region_id = data_group_id = data_event_id = data_pop_id = ""
        if data_fields is not None:
            data_region_id = data_fields[0][0].value
            data_group_id = data_fields[0][1].value
            data_event_id = data_fields[0][2].value
            data_pop_id = data_fields[0][3].value
        data_fields = [data_region_id, data_group_id, data_event_id, data_pop_id]

        feature_url = parameters[2].valueAsText
        feature_fields = parameters[3].values
        feature_region_id = ""
        if feature_fields is not None:
            feature_region_id = feature_fields[0][0].value
        
        std_pop_yr = parameters[4].valueAsText
        estimates_out = parameters[5].valueAsText
        age_std_groups = [
            ["35-44", "45-54", "55-64"],
            ["65-74", "75-84", "85up"],
            ["35-44", "45-54", "55-64", "65-74", "75-84", "85up"]
        ]
        age_std_group_names = ["35-64", "65up", "35up"]        
        ci = 0.95

        # Read in data
        data = pd.DataFrame(data = arcpy.da.SearchCursor(data_url, data_fields), columns = data_fields)
        data = data.sort_values(by = [data_region_id, data_group_id])
        regions = data[data_region_id].unique().tolist()
        num_region = data[data_region_id].nunique()
        age_groups = data[data_group_id].unique().tolist()
        num_group = data[data_group_id].nunique()

        Y = np.array(data[data_event_id]).reshape([num_region, num_group])
        n = np.array(data[data_pop_id]).reshape([num_region, num_group])
        
        # Check if data and feature regions are identical
        feature_regions = pd.DataFrame(data = arcpy.da.SearchCursor(feature_url, [feature_region_id]), columns = [feature_region_id])
        feature_regions = feature_regions[feature_region_id].unique().tolist()
        feature_regions.sort()
        if len(feature_regions) != len(regions) or feature_regions != regions:
            raise Exception("Input Table and Input Feature do not have an identical regions")

        # Create adjacency matrix
        arcpy.management.Sort(in_dataset = feature_url, out_dataset = r"in_memory\feature", sort_field = feature_region_id + " ASCENDING")
        arcpy.management.CalculateField(
            in_table = r"in_memory\feature",
            field = "NUM_REGION",
            expression = "autoIncrement()",
            expression_type = "PYTHON3",
            code_block = """rec = 0
def autoIncrement(start=1, interval=1):
    global rec
    if rec == 0:
        rec = start
    else:
        rec += interval
    return rec""",
            field_type = "LONG",
            enforce_domains = "NO_ENFORCE_DOMAINS"
        )
        arcpy.analysis.PolygonNeighbors(in_features=r"in_memory\feature", out_table=r"in_memory\adjTable", in_fields="NUM_REGION")
        adj_table = arcpy.da.TableToNumPyArray(r"in_memory\adjTable", ("src_NUM_REGION", "nbr_NUM_REGION"))
        adj_table = np.array([*adj_table.astype(object)]) 
        arcpy.management.Delete(r"in_memory\feature")
        arcpy.management.Delete(r"in_memory\adjTable")
        adj_dict = {}
        for source, neighbor in adj_table:
            if source in adj_dict:
                adj_dict[source].append(neighbor)
            else:
                adj_dict[source] = [neighbor]
        adj = list(adj_dict.values())

        # Set up standard population
        # Klein, et al. 2000 standard population
        if (std_pop_yr == "2000"): 
            std_pop = np.array([18987, 39977, 38077, 37233, 44659, 37030, 23961, 18136, 12315, 4259])
        # 2010 standard population
        if (std_pop_yr == "2010"):
            std_pop = np.array([20203362, 41025851, 43626342, 41063948, 41070606, 45006716, 36482729, 21713429, 13061122, 5493433])
        pop_ages = ["0-4", "5-14", "15-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85up"]
        std_pop = std_pop[np.isin(pop_ages, age_groups)]

        # Generate estimates
        messages.AddMessage("Generating estimates...")
        theta_out = helpers.gibbs_rucar(Y, n, adj, std_pop)
        output = helpers.expit(theta_out) * 1e5

        for ages in age_std_groups:
            output = helpers.age_std(output, age_groups, std_pop, ages)
        age_groups.extend(age_std_group_names)

        medians = helpers.get_medians(output, regions, age_groups, ci).reset_index()
        messages.AddMessage("Model finished!")

        # Output feature with joined estimates
        median_fields = medians.columns.tolist()
        median_fields[0] = data_region_id
        medians_np = np.rec.fromrecords(medians.values, names = median_fields)
        arcpy.da.NumPyArrayToTable(medians_np, estimates_out)

        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return
