# -*- coding: utf-8 -*-
import arcpy
import pandas as pd
import numpy as np
import helpers

import importlib
importlib.reload(helpers)

class Toolbox:
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Rate Stabilizing Toolbox"
        self.alias = "RST"

        # List of tool classes associated with this toolbox
        self.tools = [RST, IDP]


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
        param_data_fields.columns = [['Field', 'Region ID'], ['Field', 'Event Count'], ['Field', 'Population Count']]
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

        param_out_table = arcpy.Parameter(
            displayName="Output Table",
            name="OutputTable",
            datatype="GPTableView",
            parameterType="Required",
            direction="Output"
        )

        param_age_grp_field = arcpy.Parameter(
            displayName="Age Group Field",
            name="AgeGroupField",
            datatype="Field",
            parameterType="Optional",
            direction="Input",
            category="Age Standardization (optional)"
        )
        param_age_grp_field.parameterDependencies = [param_data_table.name]

        param_std_pop_yr = arcpy.Parameter(
            displayName="Standard Population Year",
            name="StandardPopulationYear",
            datatype="GPString",
            parameterType="Optional",
            direction="Input",
            category="Age Standardization (optional)"
        )
        param_std_pop_yr.filter.list = ["2000", "2010"]

        param_age_std_groups = arcpy.Parameter(
            displayName="Standardized Age Groups",
            name="StandardizedAgeGroups",
            datatype="GPValueTable",
            parameterType="Optional",
            direction="Input",
            category="Age Standardization (optional)"
        )
        param_age_std_groups.columns = [['String', 'Lower age value'], ['String', 'Upper age value']]
        param_age_std_groups.filters[0].type = "ValueList"
        param_age_std_groups.filters[0].list = ["0", "5", "15", "25", "35", "45", "55", "65", "75"]
        param_age_std_groups.filters[1].type = "ValueList"
        param_age_std_groups.filters[1].list = ["14", "24", "34", "44", "54", "64", "74", "84", "up"]

        param_rates_per = arcpy.Parameter(
            displayName="Rate",
            name="Rate",
            datatype="GPValueTable",
            parameterType="Required",
            direction="Input"
        )
        param_rates_per.columns = [['Long', 'Per']]
        param_rates_per.values = [[100000]]
        param_rates_per.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'

        params = [
            param_data_table,
            param_data_fields,
            param_feature,
            param_feature_field,
            param_rates_per,
            param_out_table,
            param_age_grp_field,
            param_std_pop_yr,
            param_age_std_groups
        ]
        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        data_url = parameters[0]
        data_fields = parameters[1]
        feature_url = parameters[2]
        feature_fields = parameters[3]
        rates_per = parameters[4]
        estimates_out = parameters[5]
        data_ageGrp_id = parameters[6]
        std_pop_yr = parameters[7]
        age_std_groups = parameters[8]

        # Restrict age groups to only those present within data
        age_groups_column = helpers.get_fieldAsList(data_url.valueAsText, data_ageGrp_id.valueAsText)
        if age_groups_column:
            age_groups = list(set(age_groups_column))

            grp_not_in_consts = [group for group in age_groups if group not in helpers.const_age_grps]

            if not grp_not_in_consts:
                lvs = sorted([int(group.split("-")[0]) for group in age_groups if group != "85up"])
                uvs = sorted([int(group.split("-")[1]) if group != "85up" else 85 for group in age_groups])
                age_std_groups.filters[0].list = [str(lv) for lv in lvs[:-1]]
                age_std_groups.filters[1].list = ["up" if uv == 85 else str(uv) for uv in uvs[1:]]

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""

        data_url = parameters[0]
        data_fields = parameters[1]
        feature_url = parameters[2]
        feature_fields = parameters[3]
        rates_per = parameters[4]
        estimates_out = parameters[5]
        data_ageGrp_id = parameters[6]
        std_pop_yr = parameters[7]
        age_std_groups = parameters[8]

        data_region_id_str = data_event_id_str = data_pop_id_str = None
        
        if data_fields.values is not None:
            data_fields_str = [str(field) for field in data_fields.values[0]]
            data_region_id_str = data_fields_str[0]
            data_event_id_str = data_fields_str[1]
            data_pop_id_str = data_fields_str[2]


        feature_region_id_str = None
        if feature_fields.values is not None:
            feature_region_id_str = str(feature_fields.values[0][0])

        # Check if all fields are filled in for age standardization
        if data_ageGrp_id.valueAsText is None and (std_pop_yr.valueAsText is not None or age_std_groups.valueAsText is not None):
            data_ageGrp_id.setErrorMessage("Age Group Field necessary for age standardization")
        if std_pop_yr.valueAsText is None and (data_ageGrp_id.valueAsText is not None or age_std_groups.valueAsText is not None):
            std_pop_yr.setErrorMessage("Standard Population Year necessary for age standardization")
        if age_std_groups.valueAsText is None and (data_ageGrp_id.valueAsText is not None or std_pop_yr.valueAsText is not None):
            age_std_groups.setErrorMessage("Age Groups necessary for age standardization")

        data_region_id_type = helpers.get_fieldType(data_url.valueAsText, data_region_id_str)
        feature_region_id_type = helpers.get_fieldType(data_url.valueAsText, feature_region_id_str)
        data_regions = helpers.get_fieldAsList(data_url.valueAsText, data_region_id_str)
        feature_regions = helpers.get_fieldAsList(feature_url.valueAsText, feature_region_id_str)
        # Check if Region ID types are the same
        if data_region_id_type and feature_region_id_type and data_region_id_type != feature_region_id_type:
            feature_fields.setErrorMessage("Input Feature Region ID field type does not match Input Table Region ID field type")
        elif feature_regions and data_regions:
            data_regions_set = set(data_regions)
            feature_regions_set = set(feature_regions)
            regions_only_feature = feature_regions_set - data_regions_set
            regions_only_data = data_regions_set - feature_regions_set
        # Check if data contains regions not present in feature
            if regions_only_data:
                data_fields.setErrorMessage(data_region_id_str + ' "' + str(list(regions_only_data)[0]) + '" not present in Input Feature ' + feature_region_id_str)
        # Check if feature contains regions not present in data
            elif regions_only_feature:
                feature_fields.setErrorMessage(feature_region_id_str + ' "' + str(list(regions_only_feature)[0]) + '" not present in Input Table ' + data_region_id_str)

        # Check if Event Count is an integer
        data_event_id_type = helpers.get_fieldType(data_url.valueAsText, data_event_id_str)
        if data_event_id_type and data_event_id_type not in ["SmallInteger", "Integer", "BigInteger"]:
            data_fields.setErrorMessage("Event Count field type is not an Integer")
 
        # Check if Population Count is an integer
        data_pop_id_type = helpers.get_fieldType(data_url.valueAsText, data_pop_id_str)
        if data_pop_id_type and data_pop_id_type not in ["SmallInteger", "Integer", "BigInteger"]:
            data_fields.setErrorMessage("Population Count field type is not an Integer")

        # Check if Output Table exists
        if helpers.exists(estimates_out.valueAsText):
            estimates_out.setErrorMessage("Output Table already exists")

        # Check for empty standardized age group bounds
        if age_std_groups.valueAsText is not None:
            for lv, uv, in age_std_groups.values:
                if lv == "" or uv == "":
                    age_std_groups.setErrorMessage("At least one lower or upper age value is empty")

        data_ageGrp_type = helpers.get_fieldType(data_url.valueAsText, data_ageGrp_id.valueAsText)
        ageGrp_column = helpers.get_fieldAsList(data_url.valueAsText, data_ageGrp_id.valueAsText)
        # Check if age group is a String
        if data_ageGrp_type and data_ageGrp_type != "String":
            data_ageGrp_id.setErrorMessage("Age Group field type is not a String")
        # Check if invalid age groups are present
        elif ageGrp_column:
            age_groups = list(set(ageGrp_column))
            grp_not_in_consts = [group for group in age_groups if group not in helpers.const_age_grps]
            if grp_not_in_consts:
                data_ageGrp_id.setErrorMessage('Age Group "' + grp_not_in_consts[0] + '" is not a valid age group')
            elif data_regions:
                regions_grp_dict = {}
                for i, region in enumerate(data_regions):
                    if region in regions_grp_dict:
                        regions_grp_dict[region].append(ageGrp_column[i])
                    else:
                        regions_grp_dict[region] = [ageGrp_column[i]]
                unique_grps_num_list = [len(set(regions_grp_dict[region])) for region in regions_grp_dict]
                grps_num_list = [len(regions_grp_dict[region]) for region in regions_grp_dict]
        # Check if there duplicate age groups
                if unique_grps_num_list != grps_num_list:
                    data_ageGrp_id.setErrorMessage("At least one Age Group is duplicated within a Region")
        # Check if there missing age groups
                elif len(set(grps_num_list)) != 1:
                    data_ageGrp_id.setErrorMessage("At least one Region is missing an Age Group")

        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        messages.AddMessage("Preparing data...")

        data_url = parameters[0]
        data_fields = parameters[1]
        feature_url = parameters[2]
        feature_fields = parameters[3]
        rates_per = parameters[4]
        estimates_out = parameters[5]
        data_ageGrp_id = parameters[6]
        std_pop_yr = parameters[7]
        age_std_groups = parameters[8]

        data_fields_str = [str(field) for field in data_fields.values[0]]
        if data_ageGrp_id.valueAsText is not None:
            data_fields_str.append(data_ageGrp_id.valueAsText)
        data_region_id_str = data_fields_str[0]
        data_event_id_str = data_fields_str[1]
        data_pop_id_str = data_fields_str[2]
        feature_region_id_str = str(feature_fields.values[0][0])

        rates_per = int(rates_per.valueAsText)

        # Get the age group distribution
        age_std_groups_arr = []
        age_std_groups_names = []
        if age_std_groups.values is not None:
            for lv, uv, in age_std_groups.values:
                lv_index = [i for i, grp in enumerate(helpers.const_age_grps) if grp.startswith(lv)][0]
                uv_index = [i for i, grp in enumerate(helpers.const_age_grps) if grp.endswith(uv)][0]
                age_std_groups_arr.append(helpers.const_age_grps[lv_index:(uv_index + 1)])
                age_std_groups_names.append( lv + ("to" + uv if uv != "85up" else "up") )

        # Read in data
        data = helpers.get_pandas(data_url.valueAsText, data_fields_str)
        age_groups = [""]
        num_group = 1
        data = data.sort_values(by = [data_region_id_str])
        regions = data[data_region_id_str].unique().tolist()
        num_region = data[data_region_id_str].nunique()
        if data_ageGrp_id.valueAsText is not None:
            data = data.sort_values(by = [data_region_id_str, data_ageGrp_id.valueAsText])
            age_groups = data[data_ageGrp_id.valueAsText].unique().tolist()
            num_group = data[data_ageGrp_id.valueAsText].nunique()
        elif num_region != len(data.index):
            data = data.groupby(data_region_id_str).agg({data_event_id_str : 'sum', data_pop_id_str : "sum"})
            messages.addWarningMessage("Repeated Region IDs were detected. Population and Event Counts were aggregated to totals.")
        
        Y = np.array(data[data_event_id_str]).reshape([num_region, num_group])
        n = np.array(data[data_pop_id_str]).reshape([num_region, num_group])

        std_pop = 1
        # Klein, et al. 2000 standard population
        if (std_pop_yr.valueAsText == "2000"): 
            std_pop = np.array([18987, 39977, 38077, 37233, 44659, 37030, 23961, 18136, 12315, 4259])
        # 2010 standard population
        if (std_pop_yr.valueAsText == "2010"):
            std_pop = np.array([20203362, 41025851, 43626342, 41063948, 41070606, 45006716, 36482729, 21713429, 13061122, 5493433])
        if std_pop_yr.valueAsText is not None:
            std_pop = std_pop[np.isin(helpers.const_age_grps, age_groups)]

        # Create adjacency matrix
        arcpy.AddMessage("Creating adjacency matrix ...")
        arcpy.management.Sort(in_dataset = feature_url.valueAsText, out_dataset = r"in_memory\feature", sort_field = feature_region_id_str + " ASCENDING")
        arcpy.management.AddField(in_table=r"in_memory\feature", field_name="NUM_REGION", field_type="LONG")
        with arcpy.da.UpdateCursor(r"in_memory\feature", ["NUM_REGION"]) as cursor:
            for index, row in enumerate(cursor):
                row[0] = index + 1
                cursor.updateRow(row)
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

        # Generate estimates
        messages.AddMessage("Generating estimates...")
        theta_out = helpers.gibbs_rucar(Y, n, adj, std_pop)
        output = helpers.expit(theta_out) * rates_per

        # If age standardized, generate age_groups
        if age_std_groups.values is not None:
            for ages in age_std_groups_arr:
                output = helpers.age_std(output, age_groups, std_pop, ages)
            age_groups.extend(age_std_groups_names)
        
        medians, ci_chart, reliable = helpers.get_medians(output, regions, age_groups)

        # If age standardized, combine output with prefixes, else just rename the output cols
        output = output_cols = []
        if age_std_groups.values is not None:
            medians = medians.add_prefix("median_")
            ci_chart = ci_chart.add_prefix("maxCI_")
            reliable = reliable.add_prefix("reliable_")
            output = pd.concat([medians, ci_chart, reliable], axis = 1)
            for i, med in enumerate(medians.columns):
                output_cols.extend([ medians.columns[i], ci_chart.columns[i], reliable.columns[i] ])
            output = output[output_cols].reset_index()
        else:
            output = pd.concat([medians, ci_chart, reliable], axis = 1).reset_index()
            output_cols = ["median", "maxCI", "reliable"]

        # Write out the final table
        output_cols = [data_region_id_str] + output_cols
        output_np = np.rec.fromrecords(output, names = output_cols)
        arcpy.da.NumPyArrayToTable(output_np, estimates_out.valueAsText)
    
        messages.AddMessage("Model finished!")

        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

class IDP:
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Individual Data Processing"
        self.description = ""

    def getParameterInfo(self):
        """Define the tool parameters."""

        param_byAge = arcpy.Parameter(
            displayName="By Age",
            name="ByAge",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input"
        )
        param_byAge.value = True

        param_idv_data = arcpy.Parameter(
            displayName="Input Individual Data",
            name="InputIndividualData",
            datatype="GPTableView",
            parameterType="Required",
            direction="Input"
        )

        param_idvWAge_data_fields = arcpy.Parameter(
            displayName="Input Individual Data Fields",
            name="InputIndividualDataFieldsWAge",
            datatype="GPValueTable",
            parameterType="Optional",
            direction="Input"
        )
        param_idvWAge_data_fields.parameterDependencies = [param_idv_data.name]
        param_idvWAge_data_fields.columns = [['Field', 'Region ID'], ['Field', 'Age']]
        param_idvWAge_data_fields.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'
        
        param_idvWOAge_data_fields = arcpy.Parameter(
            displayName="Input Individual Data Fields",
            name="InputIndividualDataFieldsWOAge",
            datatype="GPValueTable",
            parameterType="Optional",
            direction="Input"
        )
        param_idvWOAge_data_fields.parameterDependencies = [param_idv_data.name]
        param_idvWOAge_data_fields.columns = [['Field', 'Region ID']]
        param_idvWOAge_data_fields.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'
        param_idvWOAge_data_fields.enabled = False

        param_pop_data = arcpy.Parameter(
            displayName="Input Population Data",
            name="InputPopulationData",
            datatype="GPTableView",
            parameterType="Required",
            direction="Input"
        )

        param_popWAge_data_fields = arcpy.Parameter(
            displayName="Input Population Data Fields",
            name="InputPopulationDataFieldsWAge",
            datatype="GPValueTable",
            parameterType="Optional",
            direction="Input"
        )
        param_popWAge_data_fields.parameterDependencies = [param_pop_data.name]
        param_popWAge_data_fields.columns = [['Field', 'Region ID'], ['Field', 'Population Count'], ['Field', 'Age Group']]
        param_popWAge_data_fields.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'
        
        param_popWOAge_data_fields = arcpy.Parameter(
            displayName="Input Population Data Fields",
            name="InputPopulationDataFieldsWOAge",
            datatype="GPValueTable",
            parameterType="Optional",
            direction="Input"
        )
        param_popWOAge_data_fields.parameterDependencies = [param_pop_data.name]
        param_popWOAge_data_fields.columns = [['Field', 'Region ID'], ['Field', 'Population Count']]
        param_popWOAge_data_fields.controlCLSID = '{1A1CA7EC-A47A-4187-A15C-6EDBA4FE0CF7}'
        param_popWOAge_data_fields.enabled = False

        param_out_table = arcpy.Parameter(
            displayName="Output Table",
            name="OutputTable",
            datatype="GPTableView",
            parameterType="Required",
            direction="Output"
        )

        params = [
            param_byAge,
            param_idv_data,
            param_idvWAge_data_fields,
            param_idvWOAge_data_fields,
            param_pop_data,
            param_popWAge_data_fields,
            param_popWOAge_data_fields,
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
        
        byAge = parameters[0]
        idv_data_url = parameters[1]
        idvWAge_data_fields = parameters[2]
        idvWOAge_data_fields = parameters[3]
        pop_data_url = parameters[4]
        popWAge_data_fields = parameters[5]
        popWOAge_data_fields = parameters[6]
        out_table = parameters[7]

        idvWAge_data_fields.enabled = byAge.value
        idvWOAge_data_fields.enabled = not byAge.value
        popWAge_data_fields.enabled = byAge.value
        popWOAge_data_fields.enabled = not byAge.value

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""

        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        byAge = parameters[0]
        idv_data_url = parameters[1]
        idvWAge_data_fields = parameters[2]
        idvWOAge_data_fields = parameters[3]
        pop_data_url = parameters[4]
        popWAge_data_fields = parameters[5]
        popWOAge_data_fields = parameters[6]
        out_table = parameters[7]

        # Process individual data fields
        idv_data_fields_str = None
        idv_data_region_field = None
        idv_data_age_field = None
        if byAge.value:
            idv_data_fields_str = [str(field) for field in idvWAge_data_fields.values[0]]
            idv_data_age_field = idv_data_fields_str[1]
        else:
            idv_data_fields_str = [str(field) for field in idvWOAge_data_fields.values[0]]
        idv_data_region_field = idv_data_fields_str[0]

        # Process population data
        pop_data_fields_str = None
        pop_data_region_field = None
        pop_data_pop_field = None
        pop_data_age_field = None
        if byAge.value:
            pop_data_fields_str = [str(field) for field in popWAge_data_fields.values[0]]
            pop_data_age_field = pop_data_fields_str[2]
        else:
            pop_data_fields_str = [str(field) for field in popWOAge_data_fields.values[0]]
        pop_data_region_field = pop_data_fields_str[0]
        pop_data_pop_field = pop_data_fields_str[1]

        idv_data = helpers.get_pandas(idv_data_url.valueAsText, idv_data_fields_str)
        pop_data = helpers.get_pandas(pop_data_url.valueAsText, pop_data_fields_str)

        idv_data_groups = None
        pop_data_groups = None
        output_cols = None
        if byAge.value:
            idv_data["AgeGroup"] = idv_data[idv_data_age_field].apply(helpers.categorize_age)
            idv_data_groups = [idv_data_region_field, "AgeGroup"]
            pop_data_groups = [pop_data_region_field, pop_data_age_field]
            output_cols = ["RegionID", "AgeGroup", "EventCount", "PopCount"]
        else:
            idv_data_groups = [idv_data_region_field]
            pop_data_groups = [pop_data_region_field]
            output_cols = ["RegionID", "EventCount", "PopCount"]
        
        event_data = idv_data.groupby(idv_data_groups, as_index= False).size()
        event_data = event_data.merge(pop_data, 
            right_on= pop_data_groups, 
            left_on= idv_data_groups, 
            how = "right")
        event_data["EventCount"] = event_data["size"].fillna(value=0).astype(int)
        event_data = event_data[pop_data_groups + ["EventCount", pop_data_pop_field]]
        
        output_np = np.rec.fromrecords(event_data, names = output_cols)
        arcpy.da.NumPyArrayToTable(output_np, out_table.valueAsText)

        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""

        return